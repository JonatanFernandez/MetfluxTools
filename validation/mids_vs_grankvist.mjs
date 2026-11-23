// Validation 2: predicted mass-isotopologue distributions (MIDs) against those measured in human liver
// slices by Grankvist et al., Nat Metab 6:1963 (2024), Supplementary Data 3 (sheet "MIDs": 24-h cultured
// tissue, medium in which glucose and all 20 amino acids are uniformly 13C-labelled; MIDs not corrected for
// natural 13C). Donors D5, D8, D10, D16 and D20 (fasted medium), 3 slices each; "D20 fed" is left out.
//
// The prediction is what the application shows with: Human-GEM, HPA "liver" tissue context, Williams E
// (Grankvist) medium, the "liver" measured exchange set, biomass objective + pFBA, the 21 nutrients U-13C.
// Tissue extracts mix compartments, so each metabolite's MID is the average over its compartments weighted
// by production flux. Natural 13C (1.07 %) is added to the unlabelled carbons of each predicted isotopologue.
//
//   node tools/validation/mids_vs_grankvist.mjs <grankvist_mids.json>
import fs from "node:fs";
const R = new URL("../../", import.meta.url).pathname;
const { default: loadHighs } = await import(R + "vendor/highs/highs.mjs");
const { FbaProblem } = await import(R + "js/core/fba.js");
const { presetMedium } = await import(R + "js/core/medium.js");
const { measuredPreset, measuredNgam } = await import(R + "js/core/measured.js");
const { Simulator } = await import(R + "js/sim.js");
const { ModelIndex } = await import(R + "js/modelindex.js");
const { offGenes, blockedReactions } = await import(R + "js/core/contexts.js");
const { IsotopologueSampler } = await import(R + "js/core/isotopologues.js");
globalThis.atob ??= (s) => Buffer.from(s, "base64").toString("binary");

const measured = JSON.parse(fs.readFileSync(process.argv[2]));
const rd = (f) => JSON.parse(fs.readFileSync(R + "data/human-gem/" + f));
const model = rd("model.json"), structures = rd("structures.json"), atommaps = rd("atommaps.json"), layout = rd("layout.json");
const ctxData = rd("contexts.json");
const highs = await loadHighs({ print: () => {}, printErr: () => {} });
let problem = null;
const backend = async (msg) => { if (msg.type === "init") { problem = new FbaProblem(highs, msg.net); return true; } return problem.solve(msg.params); };
const sim = new Simulator(model, structures, atommaps, backend), index = new ModelIndex(model, layout);

const LABELLED = ["glucose", "alanine", "arginine", "asparagine", "aspartate", "cysteine", "glutamate", "glutamine", "glycine",
  "histidine", "isoleucine", "leucine", "lysine", "methionine", "phenylalanine", "proline", "serine", "threonine",
  "tryptophan", "tyrosine", "valine"];
const supplies = presetMedium(model, "williamsE").map((it) => {
  const b = model.bases[model.mets[it.met].b], st = structures[b.id];
  const lab = LABELLED.includes(b.n) && st ? st.e.map((el) => (el === "C" ? 1 : 0)) : null;
  return { ...it, label: lab, fraction: 1 };
});
const missing = LABELLED.filter((n) => !supplies.some((it) => it.label && model.bases[model.mets[it.met].b].n === n));
if (missing.length) console.log("labelled nutrients not in the medium:", missing.join(", "));

const ctx = ctxData.contexts.find((c) => c.id === "hpa_tissue:liver");
const meas = { items: measuredPreset(model, sim.exchanges, "liver").map((x) => ({ met: x.met, rate: x.rate })), tol: 0.2, ngam: measuredNgam("liver") };
const context = { id: ctx.id, mode: "block", block: blockedReactions(model, offGenes(ctx, model.genes)) };
const res = await sim.runSteady({ objective: index.biomassRxn, maximize: true, parsimonious: true, measured: meas, context }, supplies, structures);
const f = res.frames[0];
console.log(`FBA: ${f.status}${f.relaxed ? " (measured set relaxed)" : ""}, objective ${(+f.objective).toPrecision(3)}`);

if (process.argv.includes("--fluxes")) { // exchanges and a few central reactions of the solution
  for (const [met, e] of sim.exchanges) { const up = f.v[e.rxn] * e.dir; if (Math.abs(up) > 1e-5) console.log("  exchange", model.bases[model.mets[met].b].n.padEnd(28), up.toFixed(4)); }
  model.rxns.forEach((r, j) => { if (/citrate synthase|pyruvate dehydrogenase|ATP-citrate|carnitine O-palmitoyl|glutaminase|glutamine synthetase/i.test(r.n) && Math.abs(f.v[j]) > 1e-6) console.log("  rxn", r.id, r.n.slice(0, 50), f.v[j].toFixed(4)); });
}
const sampler = new IsotopologueSampler(sim.tracer, index, f.v, f.inputs, f.e, new Set(["C"]));
const NAMES = { "ala-L": "alanine", "arg-L": "arginine", "asp-L": "aspartate", cit: "citrate", "citr-L": "citrulline",
  "glc-D": "glucose", "gln-L": "glutamine", "glu-L": "glutamate", gly: "glycine", glyc3p: "sn-glycerol-3-phosphate",
  "mal-L": "malate", orn: "ornithine", "ser-L": "serine" };

const binom = (n, k, p) => { let c = 1; for (let i = 0; i < k; i++) c = (c * (n - i)) / (i + 1); return c * p ** k * (1 - p) ** (n - k); };
const addNatural = (dist, nC) => { // labelled atoms stay labelled; each of the other carbons is 13C with p = 1.07 %
  const out = new Array(nC + 1).fill(0);
  dist.forEach((p, k) => { for (let j = 0; j + k <= nC; j++) out[k + j] += p * binom(nC - k, j, 0.0107); });
  return out;
};
const meanFrac = (d) => d.reduce((s, p, k) => s + k * p, 0) / (d.length - 1);
const pattern = (d) => { const s = d.slice(1).reduce((a, b) => a + b, 0); return d.slice(1).map((p) => p / s); };

const rows = [];
for (const [id, name] of Object.entries(NAMES)) {
  const ms = measured[id];
  const meas = ms.map(([, v]) => v.slice(0, 15).reduce((a, b) => a + b, 0) / 15); // D5..D20 fasted
  const sd = ms.map(([, v]) => { const x = v.slice(0, 15), m = x.reduce((a, b) => a + b, 0) / 15; return Math.sqrt(x.reduce((s, y) => s + (y - m) ** 2, 0) / 14); });
  const tot = meas.reduce((a, b) => a + b, 0);
  const mN = meas.map((p) => p / tot), nC = mN.length - 1;
  let acc = null, wsum = 0;
  for (let m = 0; m < model.mets.length; m++) {
    if (model.bases[model.mets[m].b].n !== name || model.mets[m].c === "e") continue;
    const w = sampler._producers(sampler.ctxs[0], m).total;
    if (!(w > 0)) continue;
    const r = sampler.mid(m, 4000, Infinity);
    if (!r) continue;
    acc ??= new Array(r.dist.length).fill(0);
    r.dist.forEach((p, k) => { acc[k] += w * p; });
    wsum += w;
  }
  if (!acc) { console.log(name, "not produced in the solution"); continue; }
  const pred = addNatural(acc.map((p) => p / wsum).slice(0, nC + 1), nC);
  rows.push({ id, name, nC, meas: mN, sd, pred });
}

const pearson = (x, y) => { const n = x.length, mx = x.reduce((a, b) => a + b) / n, my = y.reduce((a, b) => a + b) / n;
  let sxy = 0, sxx = 0, syy = 0; for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; } return sxy / Math.sqrt(sxx * syy); };
const f1 = (x) => (100 * x).toFixed(1).padStart(5);
console.log("\nmetabolite                 13C fraction meas / pred   M+0 meas / pred   M+n meas / pred   pattern distance");
const allM = [], allP = [], patM = [], patP = [];
for (const r of rows) {
  const pm = pattern(r.meas), pp = pattern(r.pred), pd = 0.5 * pm.reduce((s, p, k) => s + Math.abs(p - pp[k]), 0);
  r.patternDist = pd;
  console.log(`${(r.name + " C" + r.nC).padEnd(26)} ${f1(meanFrac(r.meas))} / ${f1(meanFrac(r.pred))}          ${f1(r.meas[0])} / ${f1(r.pred[0])}     ${f1(r.meas[r.nC])} / ${f1(r.pred[r.nC])}     ${pd.toFixed(2)}`);
  console.log(`   measured  ${r.meas.map((p, k) => `M+${k} ${f1(p)}±${(100 * r.sd[k]).toFixed(1)}`).join("  ")}`);
  console.log(`   predicted ${r.pred.map((p, k) => `M+${k} ${f1(p)}`).join("  ")}`);
  r.meas.forEach((p, k) => { allM.push(p); allP.push(r.pred[k]); });
  pm.forEach((p, k) => { patM.push(p); patP.push(pp[k]); });
}
const mae = (a, b) => a.reduce((s, x, i) => s + Math.abs(x - b[i]), 0) / a.length;
console.log(`\nall isotopologue fractions (n=${allM.length}): Pearson r = ${pearson(allM, allP).toFixed(2)}, mean |error| = ${(100 * mae(allM, allP)).toFixed(1)} points`);
console.log(`labelled-isotopologue pattern, M+0 excluded (n=${patM.length}): Pearson r = ${pearson(patM, patP).toFixed(2)}, mean |error| = ${(100 * mae(patM, patP)).toFixed(1)} points`);
const fm = rows.map((r) => meanFrac(r.meas)), fp = rows.map((r) => meanFrac(r.pred));
console.log(`13C fraction per metabolite (n=${rows.length}): Pearson r = ${pearson(fm, fp).toFixed(2)}, mean |error| = ${(100 * mae(fm, fp)).toFixed(1)} points`);
fs.writeFileSync(R + "tools/validation/mids_vs_grankvist.json", JSON.stringify(rows, null, 1));
