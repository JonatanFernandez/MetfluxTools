// Validation of js/core/reginfer.js (Phase 1):
// Tests BOTH Mode A (metabolomics-only) and Mode B (metabolomics + paired RNA-seq) against real published
// experiments with known ground-truth transcriptional regulators, plus whole-sample permutation controls.
//
// Datasets:
//   1. GSE153906 (E. coli, Pal et al.): single regulator KOs (Δfnr, ΔarcA, Δihf vs WT) — RNA-seq ULM benchmark
//   2. ST002076 + GSE195954 (E. coli, Pal et al.): paired metabolomics + RNA-seq for double KOs
//      (FA = ΔfnrΔarcA, FI = ΔfnrΔihf, AI = ΔarcAΔihf vs WT) — Mode A vs Mode B head-to-head
//   3. ST004848 (mouse MPC11, methionine restriction 120 h) — Mode A (metabolomics-only) + permutation null
//   4. Dietmair 2012 (human HEK293, producer vs non-producer) — Mode A (metabolomics-only) + permutation null
//
// Usage: node tools/validation/reginfer_validate.mjs
import fs from "node:fs";
import zlib from "node:zlib";
import { contrast, normaliseByMedian, metaboliteNetwork, bottlenecks, permutationNull } from "../../js/core/metcontrast.js";
import { buildMetIndex, resolveMetabolite } from "../../js/core/metnames.js";
import { ulm, inferRegulatorsMetab, inferRegulatorsMultiOmic } from "../../js/core/reginfer.js";

const R = new URL("../../", import.meta.url).pathname;
const loadJson = (p) => JSON.parse(fs.readFileSync(p, "utf8"));
const cache = {};
const model = (s) => (cache[s] ??= loadJson(`${R}data/${s}/model.json`));
const layout = (s) => (cache["L" + s] ??= loadJson(`${R}data/${s}/layout.json`));
const regulons = (s) => (cache["R" + s] ??= loadJson(`${R}data/${s}/regulons.json`));
const thermoOf = (s) => { try { return (cache["T" + s] ??= loadJson(`${R}data/${s}/thermo.json`)).rxns; } catch { return null; } };

/** Read gzipped raw-count RNA-seq matrix and compute CPM log2FC + empirical-Bayes moderated t-stat. */
function rnaDiffFromCounts(gzPath, colsA, colsB) {
  const txt = zlib.gunzipSync(fs.readFileSync(gzPath)).toString("utf8");
  const lines = txt.split(/\r?\n/).filter((l) => l.trim());
  const head = lines[0].split("\t").map((x) => x.trim());
  const idxA = colsA.map((c) => head.indexOf(c)).filter((i) => i >= 0);
  const idxB = colsB.map((c) => head.indexOf(c)).filter((i) => i >= 0);
  const allIdx = [...idxA, ...idxB];

  // Parse raw counts & compute library sizes
  const raw = [];
  const libSize = new Map(allIdx.map((i) => [i, 0]));
  for (let i = 1; i < lines.length; i++) {
    const c = lines[i].split("\t");
    const gene = c[0]?.trim();
    if (!gene) continue;
    const vals = allIdx.map((k) => Number(c[k]) || 0);
    allIdx.forEach((k, p) => libSize.set(k, libSize.get(k) + vals[p]));
    raw.push({ gene, vals });
  }

  // log2(CPM + 1), filter low-expression genes (mean CPM < 1)
  const rows = [];
  const vars = [];
  const nA = idxA.length, nB = idxB.length;
  for (const r of raw) {
    const cpms = allIdx.map((k, p) => (1e6 * r.vals[p]) / Math.max(1, libSize.get(k)));
    const meanCpm = cpms.reduce((s, x) => s + x, 0) / cpms.length;
    if (meanCpm < 1) continue;
    const lg = cpms.map((x) => Math.log2(x + 1));
    const a = lg.slice(0, nA), b = lg.slice(nA);
    const ma = a.reduce((s, x) => s + x, 0) / nA;
    const mb = b.reduce((s, x) => s + x, 0) / nB;
    const va = nA > 1 ? a.reduce((s, x) => s + (x - ma) ** 2, 0) / (nA - 1) : 0;
    const vb = nB > 1 ? b.reduce((s, x) => s + (x - mb) ** 2, 0) / (nB - 1) : 0;
    const s2 = ((nA - 1) * va + (nB - 1) * vb) / Math.max(1, nA + nB - 2);
    vars.push(s2);
    rows.push({ gene: r.gene, log2fc: ma - mb, s2 });
  }
  vars.sort((a, b) => a - b);
  const s02 = vars[vars.length >> 1] || 0.05; // empirical-Bayes variance prior
  const seScale = Math.sqrt(1 / nA + 1 / nB);
  const map = new Map();
  for (const r of rows) {
    // moderated t-stat shrinks noisy low-df variances toward the genome-wide median s0^2
    const stat = r.log2fc / (Math.sqrt(0.5 * r.s2 + 0.5 * s02) * seScale);
    map.set(r.gene, { gene: r.gene, log2fc: r.log2fc, stat, q: NaN });
  }
  return map;
}

/** Parse Metabolomics Workbench mwTab table. */
function mwtab(path) {
  const txt = fs.readFileSync(path, "utf8");
  const factors = new Map();
  for (const line of txt.split(/\r?\n/)) {
    const m = /^SUBJECT_SAMPLE_FACTORS\s+(.*)$/.exec(line);
    if (!m) continue;
    const f = m[1].split("\t").map((x) => x.trim());
    if (f.length >= 3) factors.set(f[1], f[2]);
  }
  const block = txt.split("MS_METABOLITE_DATA_START")[1];
  const lines = block.split(/\r?\n/).filter((l) => l.trim());
  const samples = lines[0].split("\t").slice(1).map((x) => x.trim());
  const rows = [];
  for (const line of lines.slice(1)) {
    if (/^MS_METABOLITE_DATA_END/.test(line)) break;
    const cells = line.split("\t");
    const label = cells[0]?.trim();
    if (!label || /^(Samples|Factors|#)/i.test(label)) continue;
    rows.push({ label, values: cells.slice(1).map((x) => Number(x)) });
  }
  return { samples, factors, rows };
}

function splitMwtab(table, isA, isB, index) {
  const ia = table.samples.map((s, i) => [s, i]).filter(([s]) => isA(table.factors.get(s) || s)).map(([, i]) => i);
  const ib = table.samples.map((s, i) => [s, i]).filter(([s]) => isB(table.factors.get(s) || s)).map(([, i]) => i);
  const rows = table.rows.map((r) => {
    const { base } = resolveMetabolite(r.label, index);
    return { label: r.label, bases: base === undefined ? [] : [base], a: ia.map((i) => r.values[i]), b: ib.map((i) => r.values[i]) };
  });
  return { rows, nA: ia.length, nB: ib.length };
}

function prepMetContrast(group) {
  const raw = contrast(group.rows);
  const needed = Math.abs(raw.panelMedian) > 0.3;
  const rows = needed ? normaliseByMedian(group.rows).rows : group.rows;
  const c = needed ? contrast(rows) : raw;
  const fcByBase = new Map(), qByBase = new Map();
  for (const r of c.rows) if (r.bases.length) { fcByBase.set(r.bases[0], r.log2fc); qByBase.set(r.bases[0], r.q); }
  return { rows, c, fcByBase, qByBase, needed };
}

function rankOf(list, targetRegex) {
  const i = list.findIndex((x) => targetRegex.test(x.tf));
  return i >= 0 ? { rank: i + 1, of: list.length, entry: list[i] } : { rank: -1, of: list.length, entry: null };
}

// ============================================================================
// 1. GSE153906: E. coli single-KO RNA-seq ULM validation (RegulonDB)
// ============================================================================
console.log("=".repeat(78));
console.log("1. E. coli single TF knockouts — RNA-seq ULM (GSE153906 + RegulonDB)");
console.log("=".repeat(78));
const ecReg = regulons("ecoli-gem");
const gse153 = `${R}tools/cache/GSE153906_Raw_counts_msi.txt.gz`;
for (const [label, colsA, truthRe, truthName] of [
  ["Δfnr vs WT", ["Δfnr_1", "Δfnr_2"], /^FNR$/i, "FNR"],
  ["ΔarcA vs WT", ["ΔarcA_1", "ΔarcA_2"], /^ArcA$/i, "ArcA"],
  ["Δihf vs WT", ["Δihf_1", "Δihf_2"], /^IHF$/i, "IHF"],
]) {
  const rna = rnaDiffFromCounts(gse153, colsA, ["WT_1", "WT_2"]);
  const statMap = new Map([...rna].map(([g, x]) => [g, x.stat]));
  const res = ulm(ecReg, statMap, { minTargets: 5, metabolicOnly: false });
  const r = rankOf(res, truthRe);
  console.log(`\n  ${label} (${rna.size} expressed genes):`);
  console.log(`    GROUND TRUTH ${truthName}: rank ${r.rank} of ${r.of}`
    + (r.entry ? ` | score = ${r.entry.score.toFixed(2)} (expected < 0), q = ${r.entry.q.toExponential(2)}, targets ${r.entry.nObs}/${r.entry.nTotal}` : " (not found)"));
  console.log("    Top 5 TFs:", res.slice(0, 5).map((x) => `${x.tf}(${x.score > 0 ? "+" : ""}${x.score.toFixed(1)}, q=${x.q.toExponential(1)})`).join(", "));
}

// ============================================================================
// 2. ST002076 + GSE195954: Paired metabolomics + RNA-seq (double KOs)
// ============================================================================
console.log("\n" + "=".repeat(78));
console.log("2. E. coli double TF knockouts — PAIRED Metabolomics (ST002076) + RNA-seq (GSE195954)");
console.log("=".repeat(78));
const ecM = model("ecoli-gem"), ecNet = metaboliteNetwork(ecM, layout("ecoli-gem"));
const ecIdx = buildMetIndex(ecM), ecThermo = thermoOf("ecoli-gem");
const st2076 = mwtab(`${R}tools/validation/metabolomics/ST002076.mwtab`);
const gse195 = `${R}tools/cache/GSE195954_Raw_counts_msi.txt.gz`;

for (const [code, rnaCols, truths] of [
  ["FA", ["FA_1", "FA_2"], [["FNR", /^FNR$/i], ["ArcA", /^ArcA$/i]]],
  ["FI", ["FI_1", "FI_2"], [["FNR", /^FNR$/i], ["IHF", /^IHF$/i]]],
  ["AI", ["AI_1", "AI_2"], [["ArcA", /^ArcA$/i], ["IHF", /^IHF$/i]]],
]) {
  const grp = splitMwtab(st2076, (f) => new RegExp(`KO Type:${code}\\b`).test(f), (f) => /Wildtype/.test(f), ecIdx);
  const mc = prepMetContrast(grp);
  const steps = bottlenecks(ecNet, mc.fcByBase, mc.qByBase, { thermo: ecThermo });
  const rna = rnaDiffFromCounts(gse195, rnaCols, ["WT_1", "WT_2"]);

  const modeA = inferRegulatorsMetab(ecM, ecNet, ecReg, mc.fcByBase, mc.qByBase, { thermo: ecThermo, minTargets: 2 });
  const modeB = inferRegulatorsMultiOmic(ecM, ecNet, ecReg, mc.fcByBase, mc.qByBase, rna, { thermo: ecThermo, bottleneckSteps: steps });

  console.log(`\n  ${code} (${truths.map((x) => x[0]).join(" + ")} KO) vs WT:`);
  console.log(`    Metabolomics panel observability: ${modeA.nObservableRxns} rxns, ${modeA.nObservableGenes} GPR genes`);
  console.log(`    [Mode A: metabolomics-only] ${modeA.tfs.length} TFs testable (minTargets=2)`
    + (modeA.tfs.length ? " | Top 3: " + modeA.tfs.slice(0, 3).map((x) => `${x.tf}(${x.score.toFixed(2)}, ${x.nObs}g)`).join(", ") : " — panel too small to test TFs"));
  console.log(`    [Mode B: RNA-seq + metabolomics] ${modeB.tfs.length} TFs tested:`);
  for (const [tName, tRe] of truths) {
    const rk = rankOf(modeB.tfs, tRe);
    const e = rk.entry;
    console.log(`      GROUND TRUTH ${tName}: rank ${rk.rank}/${rk.of}`
      + (e ? ` | comb=${e.combinedScore.toFixed(2)}, zRna=${e.zRna.toFixed(2)} (q=${e.qRna.toExponential(1)}), concordant met chains=${e.concordant.length}` : ""));
  }
  console.log("      Top 5:", modeB.tfs.slice(0, 5).map((x) => `${x.tf}(${x.combinedScore > 0 ? "+" : ""}${x.combinedScore.toFixed(1)}, chains=${x.concordant.length})`).join(", "));
}

// ============================================================================
// 3 & 4. Mode A (metabolomics-only) + permutation null on Mouse ST004848 & Human Dietmair 2012
// ============================================================================
console.log("\n" + "=".repeat(78));
console.log("3. Mode A (metabolomics-only) + permutation null on larger eukaryotic panels");
console.log("=".repeat(78));

function testModeA(title, slug, group) {
  const m = model(slug), net = metaboliteNetwork(m, layout(slug)), reg = regulons(slug), th = thermoOf(slug);
  const mc = prepMetContrast(group);
  const res = inferRegulatorsMetab(m, net, reg, mc.fcByBase, mc.qByBase, { thermo: th, minTargets: 3 });
  const perm = res.tfs.length
    ? permutationNull(mc.rows, (rr) => {
        const cc = contrast(rr);
        const fc = new Map(), qq = new Map();
        for (const r of cc.rows) if (r.bases.length) { fc.set(r.bases[0], r.log2fc); qq.set(r.bases[0], r.q); }
        const r2 = inferRegulatorsMetab(m, net, reg, fc, qq, { thermo: th, minTargets: 3 });
        return r2.tfs.length ? Math.abs(r2.tfs[0].score) : 0;
      }, { rounds: 60 })
    : null;
  console.log(`\n  ${title} (${slug}):`);
  console.log(`    Observable: ${res.nObservableRxns} rxns, ${res.nObservableGenes} GPR genes -> ${res.tfs.length} TFs with >=3 observable targets`);
  for (const t of res.tfs.slice(0, 6)) {
    console.log(`      ${t.tf.padEnd(10)} t=${(t.score > 0 ? "+" : "") + t.score.toFixed(2)}  p=${t.p.toExponential(1)}  q=${t.q.toFixed(3)}  obs=${t.nObs}/${t.nTotal}  top=${t.targets.slice(0, 3).map((g) => `${g.gene}(${g.sign > 0 ? "+" : "-"})`).join(",")}`);
  }
  if (perm) console.log(`    Permutation null (60 rounds, whole samples): real |t|=${perm.real.toFixed(2)}, noise >= real in ${perm.better}/${perm.rounds} -> p = ${perm.p.toFixed(3)}`);
}

// Mouse ST004848 (methionine restriction 120 h)
{
  const mIdx = buildMetIndex(model("mouse-gem"));
  const tbl = mwtab(`${R}tools/validation/metabolomics/ST004848.mwtab`);
  const grp = splitMwtab(tbl, (f) => /MetR/.test(f) && /120h/.test(f), (f) => /Control/.test(f) && /120h/.test(f), mIdx);
  testModeA("ST004848 Mouse MPC11: Methionine restriction 120 h", "mouse-gem", grp);
}

// Human Dietmair 2012 (HEK293 producer vs non-producer)
{
  const hIdx = buildMetIndex(model("human-gem"));
  const d = loadJson(`${R}tools/validation/dietmair2012_intracellular.json`);
  const rows = d.metabolites.map((m) => {
    const { base } = resolveMetabolite(m.metabolite || m.label, hIdx);
    return { label: m.label, bases: base === undefined ? [] : [base], a: m.producer.map((x) => x.mean), b: m.nonProducer.map((x) => x.mean) };
  });
  testModeA("Dietmair 2012 HEK293: Producer vs Non-producer", "human-gem", { rows, nA: 4, nB: 4 });
}
