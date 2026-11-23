// Validation 5: cancer cell lines. Essential metabolic genes predicted by the app's own knockout screen
// (Simulator.screenKnockouts, with the same grouping as geneScreen() in app.js) against the CRISPR screens of DepMap
// Public 26Q1 (Chronos gene effect; tools/validation/depmap/depmap26q1_metabolic.json, made from gene_effect.csv and
// Cellosaurus DepMap cross-references).
//
// Variants of the cell type ("context"):
//   generic  no cell type: the same prediction for every line
//   block    the line's HPA cell type, reactions without a detected enzyme blocked (the app's default mode)
//   complex  the same, large complexes kept by majority and mitochondrial DNA genes never off (docs section 6)
// ("Only make costlier" cannot change essentiality: it acts after the growth optimum is fixed, so it equals generic.)
// Setting: Human-GEM, standard medium, NCI-60 measured exchange fluxes (the app's defaults), growth objective.
// A gene is predicted essential when its knockout leaves < 1 % of the growth (the app's "essential" class); genes
// switched off in the cell type are predicted non-essential (the cell already lives without them).
//
//   [STEP=6] node tools/validation/depmap_essentiality.mjs <variant> <shard k> <of n> [max lines]
//   -> tools/validation/depmap/pred_<variant>_<k>of<n>.json  { line: { growth, essential: [genes], lowered: [genes] } }
import fs from "node:fs";
const R = new URL("../../", import.meta.url).pathname;
const { default: loadHighs } = await import(R + "vendor/highs/highs.mjs");
const { FbaProblem } = await import(R + "js/core/fba.js");
const { presetMedium } = await import(R + "js/core/medium.js");
const { measuredPreset, measuredNgam } = await import(R + "js/core/measured.js");
const { Simulator } = await import(R + "js/sim.js");
const { gprActivity } = await import(R + "js/core/gpr.js");
const { offGenes, blockedReactions, isMitoGene } = await import(R + "js/core/contexts.js");

const [variant = "generic", k = "0", n = "1", maxLines = "100000"] = process.argv.slice(2);
const rd = (f) => JSON.parse(fs.readFileSync(R + "data/human-gem/" + f));
const model = rd("model.json"), st = rd("structures.json"), am = rd("atommaps.json"), ctxs = rd("contexts.json");
// line list: DEPMAP=depmap26q3_metabolic.json for the full run (the 122-line subset was drawn from 26Q1)
const D = JSON.parse(fs.readFileSync(R + "tools/validation/depmap/" + (process.env.DEPMAP || "depmap26q1_metabolic.json")));
const highs = await loadHighs({ print: () => {}, printErr: () => {} });
let pb; const backend = async (m) => { if (m.type === "init") { pb = new FbaProblem(highs, m.net); return true; } const r = pb.solve(m.params); r.ms = 0; return r; };
// Corrections under test (CORR=a,b,...; the output file name gets the suffix):
//   pool   cofactor pool of the biomass without the catabolites of retinoids, vitamin D and vitamin E
//   nocof  objective "growth without the cofactor pool" (MFL_biomass_nocof; cofactors are catalytic, their stocks outlast a screen)
//   serum  medium "standard + serum lipids" (defaultSerum: carnitine, cholesterol, palmitate and oleate of 10 % FBS)
const CORR = (process.env.CORR || "").split(",").filter(Boolean), tag = CORR.length ? "_" + CORR.join("+") : "";
if (CORR.includes("pool")) {
  const r = model.rxns.find((x) => x.id === "MAR10065"), keep = r.m.map((i) => !/retin|vitamin|chromanol|CEHC/i.test(model.bases[model.mets[i].b].n));
  r.m = r.m.filter((_, q) => keep[q]); r.c = r.c.filter((_, q) => keep[q]);
}
const sim = new Simulator(model, st, am, backend);
const sup = presetMedium(model, CORR.includes("serum") ? "defaultSerum" : "default");
const meas = { items: measuredPreset(model, sim.exchanges, "nci60").map((x) => ({ met: x.met, rate: x.rate })), tol: 0.2, ngam: measuredNgam("nci60") };
const cfg0 = { objective: model.rxns.findIndex((r) => (CORR.includes("nocof") ? r.id === "MFL_biomass_nocof" : r.obj)), maximize: true, parsimonious: true, measured: meas };
const geneRxns = new Map();
model.rxns.forEach((r, j) => { if (!r.g) return; for (const g of new Set(r.g.match(/[^\s()]+/g).filter((x) => !/^(and|or)$/i.test(x)))) { if (!geneRxns.has(g)) geneRxns.set(g, []); geneRxns.get(g).push(j); } });

async function screen(ctx) {
  let off = new Set(), block = null;
  if (ctx) {
    off = offGenes(ctx, model.genes);
    block = blockedReactions(model, off, { complexes: variant === "complex" });
    if (variant === "complex") off = new Set([...off].filter((g) => !isMitoGene(g)));
  }
  const cfg = { ...cfg0, context: block ? { id: ctx.id, mode: "block", block } : undefined };
  const f = (await sim.runSteady(cfg, sup, st)).frames[0];
  if (f.status !== "optimal" || !(f.objective > 1e-9)) return { growth: 0, status: f.status, essential: [], lowered: [] };
  const blocked = new Set(block || []);
  const groups = new Map();
  for (const [sym, rx] of geneRxns) {
    if (off.has(sym) || !rx.some((j) => Math.abs(f.v[j]) > 1e-9)) continue;
    const levels = new Map([...off].map((g) => [g, 0])); levels.set(sym, 0);
    const lost = rx.filter((j) => !blocked.has(j) && gprActivity(model.rxns[j].g, levels) === 0).sort((a, b) => a - b);
    if (!lost.length) continue;
    const key = lost.join(",");
    if (!groups.has(key)) groups.set(key, { rxns: Int32Array.from(lost), genes: [] });
    groups.get(key).genes.push(sym);
  }
  const list = [...groups.values()];
  // reference for the ratio: plain FBA optimum (the screen is plain FBA too)
  const ref = (await sim.screenKnockouts(cfg, sup, [new Int32Array(0)]))[0].objective;
  const res = await sim.screenKnockouts(cfg, sup, list.map((g) => g.rxns));
  const essential = [], lowered = [];
  res.forEach((r, i) => { const x = r.objective / ref; for (const g of list[i].genes) { if (x < 0.01) essential.push(g); else if (x < 0.95) lowered.push(g); } });
  return { growth: +ref.toFixed(5), relaxed: !!f.relaxed, lps: list.length, essential: essential.sort(), lowered: lowered.sort() };
}

const out = {};
const t0 = Date.now();
if (variant === "generic") {
  out.generic = await screen(null);
} else {
  const byId = new Map(ctxs.contexts.map((c) => [c.id, c]));
  // systematic sample: every STEP-th line in ACH order (STEP=1: all lines), then split into n shards
  const STEP = +(process.env.STEP || 1);
  const lines = Object.keys(D.lines).sort().filter((_, i) => i % STEP === 0).slice(0, +maxLines).filter((_, i) => i % +n === +k);
  for (const [i, ach] of lines.entries()) {
    out[ach] = await screen(byId.get(D.lines[ach].context));
    if (i % 10 === 0) console.error(`${variant} ${k}/${n}: ${i + 1}/${lines.length} lines, ${((Date.now() - t0) / 1000).toFixed(0)} s`);
  }
}
fs.writeFileSync(R + `tools/validation/depmap/pred_${variant}${tag}_${k}of${n}.json`, JSON.stringify(out));
console.error(`done ${variant} ${k}/${n} in ${((Date.now() - t0) / 1000).toFixed(0)} s`);
