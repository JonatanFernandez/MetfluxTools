// Validation Suite 4: E. coli iML1515 gene essentiality in aerobic M9 glucose
// against the published Keio Collection / Monk et al. 2017 (Nature Biotechnology 35:904-908,
// Supplementary Data set 11, Table 1, n = 1,503 genes) experimental growth screen AND
// gene-for-gene comparison against Monk et al.'s COBRApy computational predictions.
//
// Run: node tools/validation/ecoli_essentiality_validate.mjs
import fs from "node:fs";
const R = new URL("../../", import.meta.url).pathname;
const { default: loadHighs } = await import(R + "vendor/highs/highs.mjs");
const { FbaProblem } = await import(R + "js/core/fba.js");
const { presetMedium } = await import(R + "js/core/medium.js");
const { Simulator } = await import(R + "js/sim.js");
const { gprActivity } = await import(R + "js/core/gpr.js");

const model = JSON.parse(fs.readFileSync(R + "data/ecoli-gem/model.json", "utf8"));
const st = JSON.parse(fs.readFileSync(R + "data/ecoli-gem/structures.json", "utf8"));
const am = JSON.parse(fs.readFileSync(R + "data/ecoli-gem/atommaps.json", "utf8"));
const monk = JSON.parse(fs.readFileSync(R + "tools/cache/ecoli_ess/monk2017_glucose.json", "utf8"));

const highs = await loadHighs({ print: () => {}, printErr: () => {} });
let pb;
const backend = async (m) => {
  if (m.type === "init") { pb = new FbaProblem(highs, m.net); return true; }
  const r = pb.solve(m.params); r.ms = 0; return r;
};

const sim = new Simulator(model, st, am, backend);
const sup = presetMedium(model, "m9glc");
const objIdx = model.rxns.findIndex((r) => r.obj);
const cfg = { objective: objIdx, maximize: true, parsimonious: true };

const geneRxns = new Map();
model.rxns.forEach((r, j) => {
  if (!r.g) return;
  for (const tok of new Set(r.g.match(/[^\s()]+/g).filter((x) => !/^(and|or)$/i.test(x)))) {
    if (!geneRxns.has(tok)) geneRxns.set(tok, []);
    geneRxns.get(tok).push(j);
  }
});

// model.genes entries are [locus_id, symbol, name] tuples (e.g. ["b2551", "glyA", ""])
const sampleGene = model.genes[0];
const gprUsesBnum = geneRxns.has("b0002");
const bnumToGpr = new Map();
for (const g of model.genes) {
  if (Array.isArray(g)) {
    const [id, sym] = g;
    if (id && sym) bnumToGpr.set(id, geneRxns.has(sym) ? sym : id);
    else if (id) bnumToGpr.set(id, id);
  } else if (typeof g === "string") {
    bnumToGpr.set(g, g);
  }
}

const t0 = Date.now();
const wtFrame = (await sim.runSteady(cfg, sup, st)).frames[0];
const refObj = (await sim.screenKnockouts(cfg, sup, [new Int32Array(0)]))[0].objective;
console.log(`WT M9 aerobic glucose growth: ${refObj.toFixed(4)} h^-1 (status=${wtFrame.status}), sampleGene=${JSON.stringify(sampleGene)}, gprUsesBnum=${gprUsesBnum}, mappedBnums=${bnumToGpr.size}`);

// Group single-gene knockouts by identical lost reaction set (exact same algorithm as app.js geneScreen)
const groups = new Map();
const geneToGroup = new Map();
for (const [sym, rx] of geneRxns) {
  const levels = new Map([[sym, 0]]);
  const lost = rx.filter((j) => gprActivity(model.rxns[j].g, levels) === 0).sort((a, b) => a - b);
  if (!lost.length) continue;
  const key = lost.join(",");
  if (!groups.has(key)) groups.set(key, { rxns: Int32Array.from(lost), genes: [] });
  groups.get(key).genes.push(sym);
  geneToGroup.set(sym, key);
}
const list = [...groups.values()];
const res = await sim.screenKnockouts(cfg, sup, list.map((g) => g.rxns));
const ratioByGene = new Map();
res.forEach((r, i) => {
  const ratio = r.status === "optimal" ? r.objective / refObj : 0;
  for (const g of list[i].genes) ratioByGene.set(g, ratio);
});
const elapsedMs = Date.now() - t0;

// If bnumToGpr is empty (e.g., model.genes is array of symbols and we need b-number -> gene symbol map),
// load b-number -> symbol from RegulonDB or model metadata if present.
if (!gprUsesBnum && bnumToGpr.size === 0) {
  // Check model.geneIds or similar
  console.log("Model keys:", Object.keys(model));
  if (model.geneIds) model.geneIds.forEach((id, i) => bnumToGpr.set(id, model.genes[i]));
}

// Evaluate at both Monk's 50% growth cutoff (Computational_growth < 0.5 * ref) and MetfluxLive's 1% cutoff
function evaluate(predCut) {
  // Convention 1 (Essential = Positive, standard in DepMap/gene essentiality):
  //   TP = exp essential & pred essential
  //   FP = exp viable & pred essential
  //   TN = exp viable & pred viable
  //   FN = exp essential & pred viable
  let tp = 0, fp = 0, tn = 0, fn = 0, matched = 0, sameAsMonk = 0;
  const diffs = [];
  for (const row of monk) {
    const sym = bnumToGpr.get(row.gene);
    if (!sym && !geneRxns.has(row.gene)) continue;
    matched++;
    const key = sym || row.gene;
    const ratio = ratioByGene.has(key) ? ratioByGene.get(key) : 1.0;
    const expEss = row.measuredGrowth < row.cutoff; // measured < 0.5 => inviable/essential
    const predEss = ratio < predCut;
    const monkPredEss = row.monkCompGrowth < row.cutoff;
    if (predEss === monkPredEss) sameAsMonk++;
    else diffs.push({ gene: row.gene, sym: key, ratio: +ratio.toFixed(4), monkComp: +row.monkCompGrowth.toFixed(4), expEss });
    if (expEss && predEss) tp++;
    else if (!expEss && predEss) fp++;
    else if (!expEss && !predEss) tn++;
    else fn++;
  }
  const acc = (tp + tn) / (tp + fp + tn + fn);
  const prec = tp / (tp + fp);
  const rec = tp / (tp + fn);
  const denom = Math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
  const mcc = denom ? (tp * tn - fp * fn) / denom : 0;
  return { matched, tp, fp, tn, fn, acc, prec, rec, mcc, sameAsMonk, diffs };
}

const e50 = evaluate(0.5);
const e01 = evaluate(0.01);
console.log(`\nScreened ${geneRxns.size} genes (${list.length} unique LPs) in ${elapsedMs} ms`);
console.log(`Matched to Monk 2017 Glucose table: ${e50.matched}/${monk.length}`);
console.log(`Cutoff <50% WT (Monk 2017): TP=${e50.tp} FP=${e50.fp} TN=${e50.tn} FN=${e50.fn} | Acc=${(100*e50.acc).toFixed(2)}% Prec=${(100*e50.prec).toFixed(1)}% Rec=${(100*e50.rec).toFixed(1)}% MCC=${e50.mcc.toFixed(4)}`);
console.log(`Agreement with Monk 2017 COBRApy predictions: ${e50.sameAsMonk}/${e50.matched} (${(100*e50.sameAsMonk/e50.matched).toFixed(2)}%)`);
if (e50.diffs.length) console.log("Differences vs Monk COBRApy (first 15):", e50.diffs.slice(0, 15));
console.log(`Cutoff <1% WT (MetfluxLive UI): TP=${e01.tp} FP=${e01.fp} TN=${e01.tn} FN=${e01.fn} | Acc=${(100*e01.acc).toFixed(2)}% Prec=${(100*e01.prec).toFixed(1)}% Rec=${(100*e01.rec).toFixed(1)}% MCC=${e01.mcc.toFixed(4)}`);

fs.writeFileSync(R + "tools/validation/ecoli_essentiality_results.json", JSON.stringify({
  model: "ecoli-gem (iML1515)",
  benchmark: "Monk et al. 2017 Nature Biotechnology 35:904-908, Supplementary Data 11 Table 1 (Glucose)",
  wtGrowth: +refObj.toFixed(5),
  screenMs: elapsedMs,
  uniqueLps: list.length,
  cutoff50: { ...e50, diffs: e50.diffs.slice(0, 25) },
  cutoff01: { ...e01, diffs: undefined },
}, null, 2));
