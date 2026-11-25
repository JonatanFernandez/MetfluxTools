// Validation Suite 1: Human CollecTRI TF activity inference (ULM) on published human DESeq2 datasets
// Dataset: GSE211244 (Baniulyte et al., Mol. Cell. Biol. 2023) — HCT116 TP53+/+ and TP53-/- cells
// treated with p53 activators (Nutlin-3A, Etoposide) or ISR/ATF4 activators (Histidinol, Tunicamycin).
// Includes built-in genetic negative controls: Nutlin-3A and Etoposide in TP53-/- cells!
//
// Run: node tools/validation/reginfer_human_validate.mjs
import fs from "node:fs";
import zlib from "node:zlib";
const R = new URL("../../", import.meta.url).pathname;
const { parseRnaDiffTable, inferRegulatorsMultiOmic } = await import(R + "js/core/reginfer.js");

const model = JSON.parse(fs.readFileSync(R + "data/human-gem/model.json", "utf8"));
const regulons = JSON.parse(fs.readFileSync(R + "data/human-gem/regulons.json", "utf8"));

function loadDeseq2(gzPath, useStat = true) {
  const txt = zlib.gunzipSync(fs.readFileSync(gzPath)).toString("utf8");
  const lines = txt.split(/\r?\n/);
  const hdr = lines[0].split("\t").map((s) => s.replace(/^"|"$/g, ""));
  const iSym = hdr.indexOf("hgnc_symbol");
  const iLfc = hdr.indexOf("log2FoldChange");
  const iStat = hdr.indexOf("stat");
  const iBase = hdr.indexOf("baseMean");
  const out = ["gene\tlog2FC\tstat"];
  const seen = new Map(); // keep highest baseMean per symbol
  for (let i = 1; i < lines.length; i++) {
    if (!lines[i]) continue;
    const c = lines[i].split("\t");
    const sym = (c[iSym] || "").replace(/^"|"$/g, "").trim();
    if (!sym || sym === "NA") continue;
    const lfc = Number(c[iLfc]), stat = Number(c[iStat]), bm = Number(c[iBase]);
    if (!Number.isFinite(lfc) || !Number.isFinite(stat)) continue;
    const prev = seen.get(sym);
    if (!prev || bm > prev.bm) seen.set(sym, { sym, lfc, stat, bm });
  }
  for (const { sym, lfc, stat } of seen.values()) {
    out.push(`${sym}\t${lfc}\t${useStat ? stat : lfc}`);
  }
  return out.join("\n");
}

const CASES = [
  { id: "Nutlin_WT",     file: "GSE211244_deseq2_DMSO_vs_Nutlin_HCT116.tsv.gz",           expectedTf: "TP53", dir: +1, role: "positive (Nutlin-3A in TP53+/+)" },
  { id: "Nutlin_p53KO",  file: "GSE211244_deseq2_DMSO_vs_Nutlin_HCT116_p53KO.tsv.gz",     expectedTf: "TP53", dir:  0, role: "negative control (Nutlin-3A in TP53-/-)" },
  { id: "Etoposide_WT",  file: "GSE211244_deseq2_DMSO_vs_Etoposide_HCT116.tsv.gz",        expectedTf: "TP53", dir: +1, role: "positive (Etoposide in TP53+/+)" },
  { id: "Etoposide_KO",  file: "GSE211244_deseq2_DMSO_vs_Etoposide_HCT116_p53KO.tsv.gz",  expectedTf: "TP53", dir:  0, role: "negative control (Etoposide in TP53-/-)" },
  { id: "Histidinol_WT", file: "GSE211244_deseq2_DMSO_vs_Histidinol_HCT116.tsv.gz",       expectedTf: "ATF4", dir: +1, role: "positive (Histidinol AA-starvation ISR)" },
  { id: "Tunicamycin_WT",file: "GSE211244_deseq2_DMSO_vs_Tunicamycin_HCT116.tsv.gz",      expectedTf: "ATF4", dir: +1, role: "positive (Tunicamycin ER-stress ISR)" },
];

const results = [];
console.log("== Human CollecTRI ULM validation on GSE211244 (Baniulyte et al. 2023) ==");
for (const c of CASES) {
  const path = R + "tools/cache/geo_human/" + c.file;
  const t0 = Date.now();
  const tsv = loadDeseq2(path, true);
  const parsed = parseRnaDiffTable(tsv);
  const res = inferRegulatorsMultiOmic(model, null, regulons, null, null, parsed.map, { minTargets: 5 });
  const ms = Date.now() - t0;
  const idx = res.tfs.findIndex((x) => x.tf === c.expectedTf);
  const row = idx >= 0 ? res.tfs[idx] : null;
  const top5 = res.tfs.slice(0, 5).map((x) => `${x.tf}(t=${x.zRna.toFixed(1)},q=${x.qRna.toExponential(1)},met=${x.metGenes.length})`);
  console.log(`\n[${c.id}] ${c.role} — ${parsed.map.size} genes, ${res.tfs.length} TFs scored in ${ms} ms`);
  console.log(`  Expected ${c.expectedTf}: rank=${idx >= 0 ? idx + 1 : "NA"}/${res.tfs.length}, t=${row ? row.zRna.toFixed(2) : "NA"}, q=${row ? row.qRna.toExponential(2) : "NA"}, nTargets=${row?.nObsRna ?? 0}, metGenes=${row?.metGenes.length ?? 0}, metRxns=${row?.metRxns.length ?? 0}`);
  console.log(`  Top 5: ${top5.join(" | ")}`);
  results.push({
    id: c.id, role: c.role, expectedTf: c.expectedTf, expectedDir: c.dir,
    nGenes: parsed.map.size, nTfsScored: res.tfs.length, ms,
    rank: idx >= 0 ? idx + 1 : null,
    t: row ? +row.zRna.toFixed(3) : null,
    q: row ? row.qRna : null,
    nTargets: row?.nObsRna ?? 0,
    nMetGenes: row?.metGenes.length ?? 0,
    nMetRxns: row?.metRxns.length ?? 0,
    top5: res.tfs.slice(0, 5).map((x) => ({ tf: x.tf, t: +x.zRna.toFixed(2), q: x.qRna, metGenes: x.metGenes.length })),
  });
}

// -------------------------------------------------------------------------
// Part 2: Human Hypoxia vs Normoxia (GSE101526, Galbraith et al. 2017, HCT116 WT)
// -------------------------------------------------------------------------
console.log("\n== Human CollecTRI ULM validation on GSE101526 (HCT116 Hypoxia vs Normoxia) ==");
{
  const txt = zlib.gunzipSync(fs.readFileSync(R + "tools/cache/geo_human/GSE101526_counts.txt.gz")).toString("utf8");
  const lines = txt.split(/\r?\n/).filter(Boolean);
  const hdr = lines[0].split("\t");
  const idxNorm = hdr.map((h, i) => (/HCT116wt_norm/i.test(h) ? i : -1)).filter((i) => i > 0);
  const idxHyp = hdr.map((h, i) => (/HCT116wt_hypox/i.test(h) ? i : -1)).filter((i) => i > 0);
  const libSize = hdr.map(() => 0);
  const rawRows = [];
  for (let i = 1; i < lines.length; i++) {
    const c = lines[i].split("\t");
    if (!c[0] || c[0].startsWith("__")) continue;
    const vals = c.slice(1).map(Number);
    vals.forEach((v, k) => { if (Number.isFinite(v)) libSize[k + 1] += v; });
    rawRows.push([c[0], vals]);
  }
  const rows = [];
  const vars = [];
  for (const [sym, vals] of rawRows) {
    const gN = idxNorm.map((i) => Math.log2((vals[i - 1] / libSize[i]) * 1e6 + 0.5));
    const gH = idxHyp.map((i) => Math.log2((vals[i - 1] / libSize[i]) * 1e6 + 0.5));
    const mean = (a) => a.reduce((s, x) => s + x, 0) / a.length;
    const mN = mean(gN), mH = mean(gH);
    if (mN < 0 && mH < 0) continue; // skip unexpressed (< ~0.5 CPM)
    const vN = gN.reduce((s, x) => s + (x - mN) ** 2, 0);
    const vH = gH.reduce((s, x) => s + (x - mH) ** 2, 0);
    const s2 = (vN + vH) / (gN.length + gH.length - 2);
    rows.push({ sym, lfc: mH - mN, s2, nN: gN.length, nH: gH.length });
    vars.push(s2);
  }
  vars.sort((a, b) => a - b);
  const s02 = vars[Math.floor(vars.length / 2)] || 0.1;
  const d0 = 4;
  const tsv = ["gene\tlog2FC\tstat"];
  for (const r of rows) {
    const df = r.nN + r.nH - 2;
    const modVar = (d0 * s02 + df * r.s2) / (d0 + df);
    const t = r.lfc / Math.sqrt(modVar * (1 / r.nN + 1 / r.nH));
    tsv.push(`${r.sym}\t${r.lfc.toFixed(4)}\t${t.toFixed(4)}`);
  }
  const t0 = Date.now();
  const parsed = parseRnaDiffTable(tsv.join("\n"));
  const res = inferRegulatorsMultiOmic(model, null, regulons, null, null, parsed.map, { minTargets: 5 });
  const ms = Date.now() - t0;
  const idx = res.tfs.findIndex((x) => x.tf === "HIF1A");
  const row = idx >= 0 ? res.tfs[idx] : null;
  const idxE = res.tfs.findIndex((x) => x.tf === "EPAS1");
  const rowE = idxE >= 0 ? res.tfs[idxE] : null;
  const top5 = res.tfs.slice(0, 5).map((x) => `${x.tf}(t=${x.zRna.toFixed(1)},q=${x.qRna.toExponential(1)},met=${x.metGenes.length})`);
  console.log(`[Hypoxia_WT] HCT116 Hypoxia vs Normoxia (n=${idxHyp.length} vs ${idxNorm.length}) — ${parsed.map.size} expressed genes, ${res.tfs.length} TFs scored in ${ms} ms`);
  console.log(`  Expected HIF1A: rank=${idx >= 0 ? idx + 1 : "NA"}/${res.tfs.length}, t=${row?.zRna.toFixed(2)}, q=${row?.qRna.toExponential(2)}, nTargets=${row?.nObsRna}, metGenes=${row?.metGenes.length}`);
  console.log(`  Related EPAS1 (HIF2A): rank=${idxE >= 0 ? idxE + 1 : "NA"}/${res.tfs.length}, t=${rowE?.zRna.toFixed(2)}, q=${rowE?.qRna.toExponential(2)}`);
  console.log(`  Top 5: ${top5.join(" | ")}`);
  results.push({
    id: "Hypoxia_WT", role: "positive (HCT116 Hypoxia vs Normoxia, GSE101526)", expectedTf: "HIF1A", expectedDir: +1,
    nGenes: parsed.map.size, nTfsScored: res.tfs.length, ms,
    rank: idx >= 0 ? idx + 1 : null,
    t: row ? +row.zRna.toFixed(3) : null,
    q: row ? row.qRna : null,
    nTargets: row?.nObsRna ?? 0,
    nMetGenes: row?.metGenes.length ?? 0,
    nMetRxns: row?.metRxns.length ?? 0,
    epas1Rank: idxE >= 0 ? idxE + 1 : null,
    epas1T: rowE ? +rowE.zRna.toFixed(3) : null,
    top5: res.tfs.slice(0, 5).map((x) => ({ tf: x.tf, t: +x.zRna.toFixed(2), q: x.qRna, metGenes: x.metGenes.length })),
  });
}

fs.writeFileSync(R + "tools/validation/reginfer_human_results.json", JSON.stringify({
  datasets: ["GSE211244 (Baniulyte et al., Mol. Cell. Biol. 2023)", "GSE101526 (Galbraith et al., Cell Rep. 2017)"],
  regulons: regulons.meta,
  cases: results,
}, null, 2));
