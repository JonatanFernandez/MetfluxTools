// Validation 7: does the interpretation engine find the cause that is already known?
//
// Each dataset here is a published experiment whose cause the authors established by other means. The
// engine (js/core/metcontrast.js) is given only the metabolite table and asked for its ranked hypotheses.
// The question is narrow and falsifiable: does the known cause appear in the top three?
//
// A detector that cannot do that on an experiment designed around one enzyme or one missing nutrient has no
// business being shown to a user as a hypothesis. The permutation null answers the opposite question — can
// the engine say "nothing here" when the labels are shuffled.
//
//   node tools/validation/metcontrast_validate.mjs
// Needs tools/validation/fetch_metabolomics.sh and prepare_dietmair.py to have been run.
import fs from "node:fs";
import { contrast, metaboliteNetwork, observability, bottlenecks, pathwayShifts, ratioChanges, permutationNull, normaliseByMedian } from "../../js/core/metcontrast.js";
import { buildMetIndex, resolveMetabolite } from "../../js/core/metnames.js";

const R = new URL("../../", import.meta.url).pathname;
const cache = {};
const model = (slug) => (cache[slug] ??= JSON.parse(fs.readFileSync(`${R}data/${slug}/model.json`)));
const layout = (slug) => (cache["L" + slug] ??= JSON.parse(fs.readFileSync(`${R}data/${slug}/layout.json`)));
const thermoOf = (slug) => { try { return (cache["T" + slug] ??= JSON.parse(fs.readFileSync(`${R}data/${slug}/thermo.json`)).rxns); } catch { return null; } };

/** mwTab: the metabolite table plus the factor of every sample. */
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
  if (!block) return null;
  const lines = block.split(/\r?\n/).filter((l) => l.trim());
  const samples = lines[0].split("\t").slice(1).map((x) => x.trim());
  const rows = [];
  for (const line of lines.slice(1)) {
    if (/^(MS_METABOLITE_DATA_END|Factors)/.test(line)) continue;
    if (/^MS_METABOLITE_DATA_END/.test(line)) break;
    const cells = line.split("\t");
    const label = cells[0]?.trim();
    if (!label || /^(Samples|Factors|#)/i.test(label)) continue;
    rows.push({ label, values: cells.slice(1).map((x) => Number(x)) });
  }
  return { samples, factors, rows };
}

/** Rows for `contrast`, built from the samples whose factor matches each side. */
function twoGroups(table, isA, isB, index) {
  const ia = table.samples.map((s, i) => [s, i]).filter(([s]) => isA(table.factors.get(s) || s)).map(([, i]) => i);
  const ib = table.samples.map((s, i) => [s, i]).filter(([s]) => isB(table.factors.get(s) || s)).map(([, i]) => i);
  const rows = [];
  for (const r of table.rows) {
    const { base } = resolveMetabolite(r.label, index);
    rows.push({ label: r.label, bases: base === undefined ? [] : [base],
      a: ia.map((i) => r.values[i]), b: ib.map((i) => r.values[i]) });
  }
  return { rows, nA: ia.length, nB: ib.length };
}

function report(name, slug, group, truth) {
  const m = model(slug), net = metaboliteNetwork(m, layout(slug));
  // normalise the samples before any test: a panel that moved as a whole makes every metabolite
  // significant and flattens every pathway (see normaliseByMedian)
  // Normalising is a decision with a cost, and the panel median is what decides it. Above the threshold the
  // whole panel has moved and the per-metabolite statistics are worthless without it; below, normalising
  // removes real signal — on Dietmair, whose panel median is only -0.14, it takes the significant
  // metabolites from 18 to 2 and every candidate step with them.
  const SHIFT = 0.3; // log2, about 23 %
  const raw = contrast(group.rows);
  const needed = Math.abs(raw.panelMedian) > SHIFT;
  const norm = needed ? normaliseByMedian(group.rows) : { rows: group.rows, factors: [] };
  group = { ...group, rows: norm.rows };
  const c = needed ? contrast(group.rows) : raw;
  const spread = norm.factors.length ? `${Math.min(...norm.factors.map((f) => f.factor)).toFixed(2)}–${Math.max(...norm.factors.map((f) => f.factor)).toFixed(2)}` : "none";
  console.log(`\n${"=".repeat(78)}\n== ${name}  (${slug}, ${group.nA} vs ${group.nB} samples)`);
  console.log(needed
    ? `   panel median ${raw.panelMedian.toFixed(3)} is beyond ${SHIFT}: samples normalised to their median (factors ${spread}) -> panel median now ${c.panelMedian.toFixed(3)}.`
      + " If the total metabolite content really fell, part of a real effect went with the artefact."
    : `   panel median ${raw.panelMedian.toFixed(3)} is within ${SHIFT}: left as measured, no normalisation applied.`);
  const withBase = c.rows.filter((r) => r.bases.length);
  const fcByBase = new Map(), qByBase = new Map();
  for (const r of withBase) for (const b of r.bases) { fcByBase.set(b, r.log2fc); qByBase.set(b, r.q); }
  const obs = observability(net, fcByBase.keys());
  console.log(`   ${withBase.length} of ${c.rows.length} metabolites placed on the network; ${c.nSignificant} differ (q < 0.05)`);
  console.log(`   panel median log2 ${c.panelMedian.toFixed(3)} — a class has to beat this to be biology rather than normalisation`);
  console.log(`   the panel can argue about ${obs.nBothSides} reactions (substrate and product measured) and touches ${obs.nOneSide} more on one side`);

  const bn = bottlenecks(net, fcByBase, qByBase, { thermo: thermoOf(slug), groups: [] });
  console.log(`\n   bottlenecks, best 5 of ${bn.length}:`);
  for (const b of bn.slice(0, 5)) console.log(`     ${b.score.toFixed(2)}  ${b.id} ${b.name.slice(0, 46).padEnd(46)} ${b.ss.slice(0, 26)}${b.couldBeReversal ? "  [could be a reversal]" : ""}`);

  const ps = pathwayShifts(net, withBase, { panelMedian: c.panelMedian });
  console.log(`\n   pathways moving as a whole, best 5 of ${ps.length}:`);
  for (const p of ps.slice(0, 5)) console.log(`     q = ${p.q.toFixed(4)}  median ${p.median > 0 ? "+" : ""}${p.median.toFixed(2)} (${p.vsPanel > 0 ? "+" : ""}${p.vsPanel.toFixed(2)} vs panel, ${p.n} measured)  ${p.ss}`);

  // did it find what the authors found?
  const verdict = [];
  if (truth.pathway) {
    const rank = ps.findIndex((p) => truth.pathway.test(p.ss));
    verdict.push(`pathway ${truth.pathway} -> ${rank < 0 ? "NOT FOUND" : `rank ${rank + 1} of ${ps.length}` + (rank < 3 ? "  TOP 3" : "")}`);
  }
  if (truth.reaction) {
    const rank = bn.findIndex((b) => truth.reaction.test(b.id) || truth.reaction.test(b.name));
    verdict.push(`reaction ${truth.reaction} -> ${rank < 0 ? "NOT FOUND" : `rank ${rank + 1} of ${bn.length}` + (rank < 3 ? "  TOP 3" : "")}`);
  }
  if (truth.metabolites) {
    const hits = truth.metabolites.map((rx) => {
      const row = c.rows.find((r) => rx.test(r.label));
      return row ? `${row.label} ${row.log2fc > 0 ? "+" : ""}${row.log2fc.toFixed(2)} (q ${row.q.toExponential(1)})` : `${rx} not measured`;
    });
    verdict.push("expected metabolites: " + hits.join("; "));
  }
  console.log(`\n   KNOWN CAUSE: ${truth.what}`);
  for (const v of verdict) console.log(`     ${v}`);

  if (truth.ratios) {
    const byLabel = new Map(c.rows.map((r) => [r.label, r.log2fc]));
    const rc = ratioChanges(truth.ratios, byLabel);
    if (rc.length) { console.log("\n   ratio changes (valid with relative data):");
      for (const x of rc) console.log(`     ${x.a}/${x.b}: ${x.dlog2 > 0 ? "+" : ""}${x.dlog2.toFixed(2)} log2 (x${x.fold.toFixed(2)})  ${x.what}`); }
  }

  // can it say "nothing here"? Shuffle which replicate belongs to which condition.
  const topScore = (rows) => {
    const cc = contrast(rows);
    const fc = new Map(), qq = new Map();
    for (const r of cc.rows) for (const b of r.bases) { fc.set(b, r.log2fc); qq.set(b, r.q); }
    const b = bottlenecks(net, fc, qq, { thermo: thermoOf(slug) });
    return b.length ? b[0].score : 0;
  };
  const perm = permutationNull(group.rows, topScore, { rounds: 60 });
  console.log(`\n   permutation null (60 rounds, sample labels shuffled): best real score ${perm.real.toFixed(2)}, reached by noise in ${perm.better} rounds, p = ${perm.p.toFixed(3)}`);
  return { name, nSignificant: c.nSignificant, bottlenecks: bn.length, perm: perm.p };
}

// ---------------------------------------------------------------- the datasets ----
const results = [];
const have = (p) => fs.existsSync(p);

// 1. methionine restriction: the cause is a nutrient removed, and the route it feeds should move
{
  const p = `${R}tools/validation/metabolomics/ST004848.mwtab`;
  if (have(p)) {
    const table = mwtab(p);
    const index = buildMetIndex(model("mouse-gem"));
    // the 120 h time point, where the effect has had time to appear
    const g = twoGroups(table, (f) => /MetR/.test(f) && /120h/.test(f), (f) => /Control/.test(f) && /120h/.test(f), index);
    results.push(report("ST004848, mouse MPC11: methionine restriction against control, 120 h", "mouse-gem", g, {
      what: "methionine was removed from the medium; the transsulfuration route and the polyamines should follow",
      pathway: /methionine|cysteine|methane|sulfur/i,
      reaction: /methionine|homocyst|cystathion|adenosylmethionine/i,
      metabolites: [/^Methionine$/, /^Cystathionine$/, /^Homocysteine$/, /^Putrescine$/, /^Spermidine$/],
    }));
  } else console.log("ST004848 not downloaded, skipped");
}

// 2. E. coli regulator knockouts: absolute units, four genotypes; the machinery test
{
  const p = `${R}tools/validation/metabolomics/ST002076.mwtab`;
  if (have(p)) {
    const table = mwtab(p);
    const index = buildMetIndex(model("ecoli-gem"));
    const g = twoGroups(table, (f) => /^FA_|KO Type:(?!Wildtype)/.test(f) || /FA/.test(f), (f) => /Wildtype/.test(f), index);
    results.push(report("ST002076, E. coli: a regulator knockout against the wild type", "ecoli-gem", g, {
      what: "a transcriptional regulator was deleted; no single reaction is expected, this tests the machinery on absolute amounts",
      pathway: /amino acid|glutamate|arginine/i,
    }));
  } else console.log("ST002076 not downloaded, skipped");
}

// 3. HEK293 producer against non-producer: the cause the authors found with transcriptomics
{
  const p = `${R}tools/validation/dietmair2012_intracellular.json`;
  if (have(p)) {
    const d = JSON.parse(fs.readFileSync(p));
    const index = buildMetIndex(model("human-gem"));
    const rows = d.metabolites.map((m) => {
      const { base } = resolveMetabolite(m.metabolite || m.label, index);
      return { label: m.label, bases: base === undefined ? [] : [base],
        a: m.producer.map((x) => x.mean), b: m.nonProducer.map((x) => x.mean) };
    });
    results.push(report("Dietmair 2012, HEK293: producer against non-producer", "human-gem", { rows, nA: 4, nB: 4 }, {
      what: "the producer line has oxidative phosphorylation genes down (their transcriptomics) and consumes less glucose and glutamine",
      reaction: /ATP|adenylate|nucleoside.?diphosphate/i,
      metabolites: [/^ATP$/, /^AMP$/, /^LAC$/, /^GLN$/],
      ratios: [["ATP", "AMP", "energy charge"], ["LAC", "PYR", "cytosolic NADH/NAD+"]],
    }));
  } else console.log("Dietmair not prepared, skipped");
}

console.log(`\n${"=".repeat(78)}\nsummary`);
for (const r of results) console.log(`  ${r.name.slice(0, 56).padEnd(56)} ${String(r.nSignificant).padStart(3)} changed, ${String(r.bottlenecks).padStart(4)} candidate steps, permutation p ${r.perm.toFixed(3)}`);
