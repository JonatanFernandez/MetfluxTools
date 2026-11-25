// Does an unlabelled endogenous amino acid pool explain the labelling of human liver slices?
//
// Grankvist et al. (Nat Metab 6:1963, 2024) labelled glucose AND all 20 amino acids of the medium, so every
// unlabelled carbon the tissue shows has to come from inside it. Without such a source the simulation predicts
// nine of eleven metabolites fully labelled (M+0 = 0) where 25-65 % is measured. This script asks whether one
// number -- the rate of protein turnover -- accounts for that, and whether the number it needs is physiological.
//
// PROTEIN TURNOVER IS A CYCLE, NOT A SOURCE. A liver slice that is not growing degrades and resynthesises
// protein at the same rate, so there is no net amino acid gain: modelling the pool as a one-way influx would
// invent mass that the model must then dispose of through urea, and the measured urea output (0.012) already
// balances the measured amino acid uptake. Both halves are therefore added:
//
//     MFL_P1_deg :      -> sum_a w_a AA_a      old protein, carbon unlabelled
//     MFL_P1_syn : sum_a w_a AA_a ->           new protein, drawn from the mixed pool
//
// both fixed at lb = ub = F. Net mass flux is zero and the measured exchanges are untouched; what changes is
// that the free amino acid pool is now fed partly by unlabelled carbon, which is exactly the dilution.
//
// ONE COUPLED COLUMN, NOT TWENTY INJECTIONS. The engine's `sources` mechanism cannot express this: fba.js
// gives a source only an upper bound (`up[n2+i] = s.max`, lower bound left at 0) and zero pFBA cost, so twenty
// independent injections would be twenty free parameters that pFBA would prefer over paying for medium uptake.
// A single coupled reaction with a fixed flux is one degree of freedom and obeys the protein composition.
// The reactions are added to the in-memory model only; data/human-gem/model.json is never touched.
//
// Usage: node tools/validation/liver_pools.mjs tools/validation/grankvist_mids.json [--sweep a,b,c] [--fluxes]
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

const donorArg = process.argv.find((a) => a.startsWith("--donors="));
const donors = donorArg ? donorArg.split("=")[1].split(",").map(Number) : [0, 1, 2, 3, 4];
const measured = JSON.parse(fs.readFileSync(process.argv[2] || R + "tools/validation/grankvist_mids.json"));
const rd = (f) => JSON.parse(fs.readFileSync(R + "data/human-gem/" + f));
const model = rd("model.json"), structures = rd("structures.json"), atommaps = rd("atommaps.json"), layout = rd("layout.json");
const ctxData = rd("contexts.json");

// ---------------------------------------------------------------- the pool, from the model's own protein
// The amino acids and proportions are those of the biomass reaction, so the pool releases what this model
// says liver protein is made of. Nothing here is fitted.
const AA_NAMES = new Set(["alanine", "arginine", "asparagine", "aspartate", "cysteine", "glutamine", "glutamate",
  "glycine", "histidine", "isoleucine", "leucine", "lysine", "methionine", "phenylalanine", "proline", "serine",
  "threonine", "tryptophan", "tyrosine", "valine"]);
const baseName = (i) => model.bases[model.mets[i].b].n;
const bio = model.rxns.find((r) => r.id === "MAR04413");
const poolMets = [], poolW = [];
bio.m.forEach((i, k) => {
  const n = baseName(i).replace(/^L-/, "").toLowerCase();
  if (bio.c[k] < 0 && AA_NAMES.has(n)) { poolMets.push(i); poolW.push(-bio.c[k]); }
});
const AA_TOTAL = poolW.reduce((s, x) => s + x, 0);        // mmol amino acid per gDW of biomass protein
const w = poolW.map((x) => x / AA_TOTAL);                  // composition, sums to 1

// The degrading half must reach the TRACER, and the tracer is given only exchanges and injections
// (sim.js, traceFlux). A plain added reaction carries flux in the LP but propagates no label -- it has no
// substrate, so there is no atom to draw an enrichment from and the amino acid it makes stays as labelled as
// the medium's. The degradation is therefore supplied as one pinned injection per amino acid (min = max =
// F*w_a, which needs the `min` field added to fba.js), each unlabelled by virtue of carrying no label.
// Twenty pinned columns are still ONE degree of freedom: all twenty are determined by F and the fixed
// composition, and being pinned they cannot be chosen by the solver or preferred by pFBA's zero cost.
//
// The resynthesising half is a single coupled reaction, since it only has to remove what was released.
const nRxn0 = model.rxns.length;
model.rxns.push({ id: "MFL_P1_syn", n: "protein synthesis (drain of the mixed amino acid pool)", m: poolMets.slice(), c: w.map((x) => -x), lb: 0, ub: 0, ss: bio.ss, g: "", ec: [] });
const SYN = nRxn0;

const highs = await loadHighs({ print: () => {}, printErr: () => {} });
let problem = null;
const backend = async (msg) => { if (msg.type === "init") { problem = new FbaProblem(highs, msg.net); return true; } return problem.solve(msg.params); };
const sim = new Simulator(model, structures, atommaps, backend), index = new ModelIndex(model, layout);

const LABELLED = ["glucose", "alanine", "arginine", "asparagine", "aspartate", "cysteine", "glutamate", "glutamine",
  "glycine", "histidine", "isoleucine", "leucine", "lysine", "methionine", "phenylalanine", "proline", "serine",
  "threonine", "tryptophan", "tyrosine", "valine"];
const supplies = presetMedium(model, "williamsE").map((it) => {
  const b = model.bases[model.mets[it.met].b], st = structures[b.id];
  const lab = LABELLED.includes(b.n) && st ? st.e.map((el) => (el === "C" ? 1 : 0)) : null;
  return { ...it, label: lab, fraction: 1 };
});
const ctx = ctxData.contexts.find((c) => c.id === "hpa_tissue:liver");
const meas = { items: measuredPreset(model, sim.exchanges, "liver").map((x) => ({ met: x.met, rate: x.rate })), tol: 0.2, ngam: measuredNgam("liver") };
const context = { id: ctx.id, mode: "block", block: blockedReactions(model, offGenes(ctx, model.genes)) };

// what the medium supplies, for the dilution a given F predicts a priori
const aaUptake = meas.items.reduce((s, it) => {
  const n = baseName(it.met).replace(/^L-/, "").toLowerCase();
  return AA_NAMES.has(n) && it.rate > 0 ? s + it.rate : s;
}, 0);

// ---------------------------------------------------------------- measured data, as in mids_vs_grankvist.mjs
const NAMES = { "ala-L": "alanine", "arg-L": "arginine", "asp-L": "aspartate", cit: "citrate", "citr-L": "citrulline",
  "glc-D": "glucose", "gln-L": "glutamine", "glu-L": "glutamate", gly: "glycine", "mal-L": "malate", orn: "ornithine" };
// Glucose and citrate fail in the opposite direction (predicted LESS labelled than measured): their error is in
// the flux solution, not the dilution, so they are reported apart and act as the control -- an unlabelled pool
// must not improve them. If it does, the fit is absorbing flux error and the result cannot be trusted.
const FLUX_CONTROL = new Set(["glucose", "citrate"]);

/** Mean MID over a chosen set of donors. The first 15 columns are D5, D8, D10, D16, D20 fasted, 3 slices each. */
function donorMean(ms, donorIdx) {
  const cols = donorIdx.flatMap((d) => [3 * d, 3 * d + 1, 3 * d + 2]);
  const raw = ms.map(([, v]) => cols.reduce((a, c) => a + v[c], 0) / cols.length);
  const tot = raw.reduce((a, b) => a + b, 0);
  return raw.map((p) => p / tot);
}
const addNatural = (dist, nC) => { // labelled atoms stay labelled; every other carbon is 13C with p = 1.07 %
  const out = new Array(nC + 1).fill(0);
  dist.forEach((p, k) => { for (let j = 0; j + k <= nC; j++) out[k + j] += p * binom(nC - k, j, 0.0107); });
  return out;
};
const binom = (n, k, p) => { let c = 1; for (let i = 0; i < k; i++) c = (c * (n - i)) / (i + 1); return c * p ** k * (1 - p) ** (n - k); };
const pearson = (a, b) => {
  const ma = a.reduce((s, x) => s + x, 0) / a.length, mb = b.reduce((s, x) => s + x, 0) / b.length;
  let sa = 0, sb = 0, sab = 0;
  for (let i = 0; i < a.length; i++) { const da = a[i] - ma, db = b[i] - mb; sa += da * da; sb += db * db; sab += da * db; }
  return sab / Math.sqrt(sa * sb || 1e-12);
};
const mae = (a, b) => (a.length ? a.reduce((s, x, i) => s + Math.abs(x - b[i]), 0) / a.length : NaN);
const pattern = (mid) => { const s = mid.slice(1).reduce((t, x) => t + x, 0); return s > 1e-9 ? mid.slice(1).map((x) => x / s) : mid.slice(1).map(() => 0); };

/** Run the whole prediction with the pool at flux F and score it against the measurements. */
async function evaluate(F) {
  // The Simulator snapshots the bounds in its constructor (sim.js:16), so set them there.
  sim.baseLb[SYN] = sim.baseUb[SYN] = F;
  // one pinned, unlabelled injection per amino acid, in the proportions of protein
  const pool = poolMets.map((met, k) => ({ met, conc: 1, vmax: F * w[k], km: 0, min: F * w[k], label: null, fraction: 1 }));
  const supply = F > 0 ? supplies.concat(pool) : supplies;
  const res = await sim.runSteady({ objective: index.biomassRxn, maximize: true, parsimonious: true, measured: meas, context }, supply, structures);
  const f = res.frames[0];
  if (!/optimal/i.test(f.status)) return { F, status: f.status };
  const sampler = new IsotopologueSampler(sim.tracer, index, f.v, f.inputs, f.e, new Set(["C"]));
  const rows = [];
  for (const [id, name] of Object.entries(NAMES)) {
    const ms = measured[id];
    if (!ms) continue;
    const mN = donorMean(ms, donors), nC = mN.length - 1;
    // tissue extracts mix compartments: average each compartment's MID weighted by its production flux
    let acc = null, wsum = 0;
    for (let m = 0; m < model.mets.length; m++) {
      if (model.bases[model.mets[m].b].n !== name || model.mets[m].c === "e") continue;
      const wgt = sampler._producers(sampler.ctxs[0], m).total;
      if (!(wgt > 0)) continue;
      const r = sampler.mid(m, 4000, Infinity);
      if (!r) continue;
      acc ??= new Array(r.dist.length).fill(0);
      r.dist.forEach((p, k) => { acc[k] += wgt * p; });
      wsum += wgt;
    }
    if (!acc) continue;
    const pred = addNatural(acc.map((p) => p / wsum).slice(0, nC + 1), nC);
    // every donor separately, for leave-one-donor-out scoring later
    const perDonor = [0, 1, 2, 3, 4].map((d) => donorMean(ms, [d]));
    rows.push({ name, meas: mN, pred, perDonor });
  }
  return { F, status: f.status, objective: +f.objective, rows, v: f.v };
}

// ---------------------------------------------------------------- the sweep
const pct = (F) => (100 * F * 24 / AA_TOTAL);                     // F -> % of protein turned over per day
const sweepArg = process.argv.find((a) => a.startsWith("--sweep"));
const GRID = sweepArg ? sweepArg.split("=")[1].split(",").map(Number)
  : [0, 0.004, 0.008, 0.0133, 0.020, 0.0267, 0.040];               // 0, and 1.5 / 3 / 5 / 7.5 / 10 / 15 %/day

console.log(`liver protein pool, from the model's own biomass: ${poolMets.length} amino acids, ${AA_TOTAL.toFixed(3)} mmol AA/gDW`);
console.log(`measured amino acid uptake from the medium      : ${aaUptake.toFixed(5)} mmol/gDW/h`);
console.log(`physiological liver protein turnover 3-10 %/day : F = ${(AA_TOTAL * 0.03 / 24).toFixed(5)} - ${(AA_TOTAL * 0.10 / 24).toFixed(5)} mmol AA/gDW/h`);
console.log(`\n  F (mmol/gDW/h)   %/day   predicted dilution   status`);
for (const F of GRID) {
  const d = F / (F + aaUptake);
  console.log(`  ${F.toFixed(5).padStart(12)}   ${pct(F).toFixed(1).padStart(5)}   ${(100 * d).toFixed(0).padStart(17)} %`);
}

const results = [];
for (const F of GRID) {
  const r = await evaluate(F);
  if (!r.rows) { console.log(`\nF = ${F}: ${r.status}`); results.push(r); continue; }
  const dil = r.rows.filter((x) => !FLUX_CONTROL.has(x.name));
  const ctl = r.rows.filter((x) => FLUX_CONTROL.has(x.name));
  const flat = (rs, pick) => rs.flatMap(pick);
  const allM = flat(r.rows, (x) => x.meas), allP = flat(r.rows, (x) => x.pred);
  const dilM = flat(dil, (x) => x.meas), dilP = flat(dil, (x) => x.pred);
  const ctlM = flat(ctl, (x) => x.meas), ctlP = flat(ctl, (x) => x.pred);
  const patM = flat(r.rows, (x) => pattern(x.meas)), patP = flat(r.rows, (x) => pattern(x.pred));
  r.score = {
    all: 100 * mae(allM, allP), dilution: 100 * mae(dilM, dilP), control: 100 * mae(ctlM, ctlP),
    pattern: 100 * mae(patM, patP), r: pearson(allM, allP),
  };
  results.push(r);
  console.log(`\nF = ${F.toFixed(5)} (${pct(F).toFixed(1)} %/day), growth ${r.objective.toPrecision(3)}`);
  console.log(`   all ${r.score.all.toFixed(1)} pts (r = ${r.score.r.toFixed(2)})   dilution-limited (${dil.length} mets) ${r.score.dilution.toFixed(1)}   flux control (glucose, citrate) ${r.score.control.toFixed(1)}   pattern without M+0 ${r.score.pattern.toFixed(1)}`);
  console.log("   " + dil.map((x) => `${x.name.slice(0, 4)} M+0 ${(100 * x.pred[0]).toFixed(0)}/${(100 * x.meas[0]).toFixed(0)}`).join("  "));
}

// ---------------------------------------------------------------- verdict
const ok = results.filter((r) => r.score);
if (ok.length > 1) {
  const base = ok.find((r) => r.F === 0);
  const inRange = ok.filter((r) => pct(r.F) >= 3 && pct(r.F) <= 10);
  const bestPhys = inRange.length ? inRange.reduce((a, b) => (b.score.dilution < a.score.dilution ? b : a)) : null;
  const bestAny = ok.reduce((a, b) => (b.score.dilution < a.score.dilution ? b : a));
  const f2 = (x) => x.toFixed(1).padStart(5);
  console.log(`\n${"=".repeat(92)}`);
  console.log(`  F (%/day)        all   dilution-limited   flux control   pattern without M+0`);
  for (const r of ok) console.log(`  ${pct(r.F).toFixed(1).padStart(6)}      ${f2(r.score.all)}            ${f2(r.score.dilution)}          ${f2(r.score.control)}                 ${f2(r.score.pattern)}`);

  console.log(`\nWhat the three columns are for:`);
  console.log(`  dilution-limited  the nine metabolites the pool is meant to fix. It should fall.`);
  console.log(`  flux control      glucose and citrate, whose error is of the opposite sign. It must NOT fall:`);
  console.log(`                    an unlabelled pool cannot help them, so if it does, the fit is absorbing flux error.`);
  console.log(`  pattern (no M+0)  the shape of the labelled isotopologues. A pool adds M+0 only, so it should barely`);
  console.log(`                    move; what is left needs reversible exchange fluxes that one pFBA solution lacks.`);

  if (bestPhys) {
    const gainPhys = base.score.dilution - bestPhys.score.dilution, gainAll = base.score.dilution - bestAny.score.dilution;
    console.log(`\nInside the physiological range (3-10 %/day): best at ${pct(bestPhys.F).toFixed(1)} %/day,`);
    console.log(`  dilution-limited ${base.score.dilution.toFixed(1)} -> ${bestPhys.score.dilution.toFixed(1)} points (${(100 * gainPhys / gainAll).toFixed(0)} % of everything the whole grid achieves).`);
  }
  const ctlOk = bestAny.score.control >= base.score.control - 0.5;
  console.log(`\nControl: ${base.score.control.toFixed(1)} -> ${bestAny.score.control.toFixed(1)} across the grid -- ${ctlOk ? "never improves, as required." : "IMPROVES, which invalidates the result: the fit is absorbing flux error."}`);
  console.log(`Monotone beyond the physiological range: the residual M+0 that steady-state turnover cannot make is`);
  console.log(`  what isotopic non-stationarity at 24 h would add, so the tail is expected and should not be fitted away.`);
  fs.writeFileSync(R + "tools/validation/liver_pools.json", JSON.stringify(
    ok.map((r) => ({ F: r.F, pctPerDay: pct(r.F), score: r.score, rows: r.rows.map((x) => ({ name: x.name, pred: x.pred, perDonor: x.perDonor })) })), null, 1));
  console.log(`\nPredictions written to tools/validation/liver_pools.json for cross-validation.`);
}
