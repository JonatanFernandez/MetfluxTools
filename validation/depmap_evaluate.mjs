// Scores the predictions of depmap_essentiality.mjs against DepMap Public 26Q3 (tools/validation/depmap/depmap26q3_metabolic.json,
// made by depmap_prepare.py). Three definitions of "essential in a line": Chronos gene effect <= -0.5, <= -1, and
// probability of dependency >= 0.5 (the DepMap portal's call). Genes: the metabolic genes of Human-GEM screened by DepMap.
// Besides all genes, the scores are given without DepMap's common essentials (genes needed by nearly every line), and for
// the line-specific genes (essential in 10-90 % of the lines scored), where only a cell type can make a difference.
//   node tools/validation/depmap_evaluate.mjs  -> prints the tables, writes tools/validation/depmap/evaluation.json
import fs from "node:fs";
const R = new URL("../../", import.meta.url).pathname, DIR = R + "tools/validation/depmap/";
const D = JSON.parse(fs.readFileSync(DIR + "depmap26q3_metabolic.json"));
// TAG: the corrections of the run (e.g. TAG=nocof+serum reads pred_<variant>_nocof+serum_<k>of<n>.json)
const TAG = process.env.TAG ? "_" + process.env.TAG : "";
const load = (v) => Object.assign({}, ...fs.readdirSync(DIR).filter((f) => new RegExp(`^pred_${v}${TAG.replace("+", "\\+")}_\\d+of\\d+\\.json$`).test(f) && !f.includes("0of1")).map((f) => JSON.parse(fs.readFileSync(DIR + f))));
const generic = JSON.parse(fs.readFileSync(DIR + (process.env.GENERIC || `pred_generic${TAG}_0of1.json`))).generic, block = load("block"), complex = load("complex");
// a variant still running (no pred files yet) is left out
const lines = Object.keys(block).filter((a) => (!Object.keys(complex).length || complex[a]) && D.effect[a]).sort();
const common = new Set(D.commonEssential);
const mcc = ({ tp, fp, fn, tn }) => { const d = Math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)); return d ? (tp * tn - fp * fn) / d : 0; };
const zero = () => ({ tp: 0, fp: 0, fn: 0, tn: 0 });
const add = (a, b) => { for (const k of ["tp", "fp", "fn", "tn"]) a[k] += b[k]; };
const fmt = (c) => ({ ...c, precision: +(c.tp / (c.tp + c.fp || 1)).toFixed(3), recall: +(c.tp / (c.tp + c.fn || 1)).toFixed(3), mcc: +mcc(c).toFixed(3) });
const PRED = { generic: () => new Set(generic.essential), block: (a) => new Set(block[a].essential), complex: (a) => new Set(complex[a].essential) };
if (!Object.keys(complex).length) { delete PRED.complex; PRED.complex = PRED.block; }
const TRUTH = { "Chronos <= -0.5": (a, i) => (D.effect[a][i] === null ? null : D.effect[a][i] <= -0.5), "Chronos <= -1": (a, i) => (D.effect[a][i] === null ? null : D.effect[a][i] <= -1),
  "P(dependency) >= 0.5": (a, i) => (D.prob[a][i] === null ? null : D.prob[a][i] >= 0.5) };
const out = { lines: lines.length, genes: D.genes.length, commonEssential: common.size, results: {} };
console.log(`${lines.length} lines x ${D.genes.length} metabolic genes (${common.size} DepMap common essentials)`);
for (const [tname, truth] of Object.entries(TRUTH)) {
  const res = {};
  for (const [pname, pred] of Object.entries(PRED)) {
    const all = zero(), noCommon = zero();
    const byLineage = {};
    for (const a of lines) {
      const p = pred(a), c = zero(), cn = zero();
      D.genes.forEach((g, i) => { const t = truth(a, i); if (t === null) return; const k = p.has(g) ? (t ? "tp" : "fp") : t ? "fn" : "tn"; c[k]++; if (!common.has(g)) cn[k]++; });
      add(all, c); add(noCommon, cn);
      const L = D.lines[a].lineage || "?"; (byLineage[L] ||= zero()); add(byLineage[L], c);
    }
    // line-specific genes: essential in 10-90 % of these lines
    const sel = D.genes.map((g, i) => { const t = lines.map((a) => truth(a, i)).filter((x) => x !== null); const e = t.filter(Boolean).length / (t.length || 1); return e >= 0.1 && e <= 0.9 ? i : -1; }).filter((i) => i >= 0);
    const spec = zero(); let varying = 0, right = 0;
    for (const i of sel) {
      const g = D.genes[i], cc = zero();
      for (const a of lines) { const t = truth(a, i); if (t === null) continue; cc[pred(a).has(g) ? (t ? "tp" : "fp") : t ? "fn" : "tn"]++; }
      add(spec, cc);
      if (cc.tp + cc.fp > 0 && cc.fn + cc.tn > 0) { varying++; if (mcc(cc) > 0) right++; }
    }
    res[pname] = { all: fmt(all), withoutCommon: fmt(noCommon), lineSpecific: { genes: sel.length, ...fmt(spec), predictedToVary: varying, varyInRightDirection: right },
      byLineage: Object.fromEntries(Object.entries(byLineage).map(([k, v]) => [k, +mcc(v).toFixed(3)])) };
  }
  out.results[tname] = res;
  console.log(`\n== essential = ${tname}`);
  for (const [p, r] of Object.entries(res)) console.log(`  ${p.padEnd(8)} all: precision ${r.all.precision} recall ${r.all.recall} MCC ${r.all.mcc} | without common essentials: MCC ${r.withoutCommon.mcc} (TP ${r.withoutCommon.tp}, FP ${r.withoutCommon.fp}) | line-specific genes (${r.lineSpecific.genes}): MCC ${r.lineSpecific.mcc}, predicted to vary ${r.lineSpecific.predictedToVary}, of which in the right direction ${r.lineSpecific.varyInRightDirection}`);
}
// most frequent errors (complex variant, Chronos <= -0.5)
const fp = new Map(), fn = new Map();
for (const a of lines) { const p = PRED.complex(a); D.genes.forEach((g, i) => { const x = D.effect[a][i]; if (x === null) return; if (p.has(g) && x > -0.5) fp.set(g, (fp.get(g) || 0) + 1); if (!p.has(g) && x <= -0.5) fn.set(g, (fn.get(g) || 0) + 1); }); }
out.topFalseEssential = [...fp].sort((a, b) => b[1] - a[1]).slice(0, 30);
out.topMissedEssential = [...fn].sort((a, b) => b[1] - a[1]).slice(0, 30);
console.log(`\npredicted essential but not in DepMap (number of lines): ${out.topFalseEssential.slice(0, 20).map(([g, k]) => `${g} ${k}`).join(", ")}`);
console.log(`essential in DepMap but missed: ${out.topMissedEssential.slice(0, 20).map(([g, k]) => `${g} ${k}`).join(", ")}`);
fs.writeFileSync(DIR + `evaluation${TAG}.json`, JSON.stringify(out, null, 1));
