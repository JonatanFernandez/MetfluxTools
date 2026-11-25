// Tests of the exometabolomics modules: js/core/exodata.js (rates from a time course) and
// js/core/fitexchange.js (fit of measured rates to the network).
//   node tools/test_exofit.mjs
import fs from "node:fs";
import loadHighs from "../vendor/highs/highs.mjs";
import { packNetwork } from "../js/core/fba.js";
import { presetMedium, findExchanges, applyMedium, uptakeCapacity } from "../js/core/medium.js";
import { measuredPreset } from "../js/core/measured.js";
import { biomassIntegral, rateFromCourse, ratesFromReplicates, growthRate, toGdw, decayConstant,
  parseExoTable, groupForRates } from "../js/core/exodata.js";
import { ExchangeFit, fitBounds, exchangeBalance } from "../js/core/fitexchange.js";

let bad = 0;
const ok = (c, msg) => { console.log((c ? "ok  " : "FAIL") + " " + msg); if (!c) bad++; };
const near = (a, b, tol, msg) => ok(Math.abs(a - b) <= tol, `${msg} (${a.toPrecision(6)} vs ${b.toPrecision(6)}, tol ${tol})`);

// ---------------------------------------------------------------- 1. biomass integral ----
console.log("\n== biomass integral (logarithmic mean per sub-interval)");
for (const doublings of [1, 2, 3]) {
  const dt = 24, mu = (Math.log(2) * doublings) / dt, X1 = 1e-4, X2 = X1 * 2 ** doublings;
  const exact = (X2 - X1) / mu;
  near(biomassIntegral([0, dt], [X1, X2]), exact, exact * 1e-12, `${doublings} doubling(s), 2 points: exact`);
  const arith = ((X1 + X2) / 2) * dt;
  console.log(`     the arithmetic mean would overestimate ∫X dt by ${(100 * (arith / exact - 1)).toFixed(1)} %, so underestimate q by as much`);
}
{ // adding an intermediate sample must not make it worse (it does with plain trapezoids)
  const mu = Math.log(2) / 24, X = (t) => 1e-4 * Math.exp(mu * t);
  const exact = (X(48) - X(0)) / mu;
  const t3 = [0, 24, 48], x3 = t3.map(X);
  let trap = 0; for (let i = 1; i < 3; i++) trap += ((x3[i] + x3[i - 1]) / 2) * (t3[i] - t3[i - 1]);
  near(biomassIntegral(t3, x3), exact, exact * 1e-12, "3 points on an exponential curve: still exact");
  console.log(`     plain trapezoids over the same 3 points would be ${(100 * (trap / exact - 1)).toFixed(2)} % off: the extra sample would have hurt`);
  near(biomassIntegral([0, 48], [X(0), X(48)]), biomassIntegral(t3, x3), exact * 1e-12, "no jump between 2 and 3 time points");
}
{ // a culture that stops growing: follow the measured curve, not an exponential
  const t = [0, 12, 24, 36], X = [1e-4, 2e-4, 3e-4, 3.1e-4];
  const got = biomassIntegral(t, X);
  const asExp = (X[3] - X[0]) / (Math.log(X[3] / X[0]) / 36);
  ok(got > 0 && Math.abs(asExp / got - 1) > 0.05, `confluent culture: ${got.toExponential(3)}; as an exponential it would be ${(100 * (asExp / got - 1)).toFixed(1)} % off`);
}

// ---------------------------------------------------------------- 2. rate recovery ----
console.log("\n== rate from a time course");
{
  const V = 2, dt = 24, mu = Math.log(4) / dt, X1 = 1e-4, X2 = X1 * 4;
  const qTrue = 0.35, integral = (X2 - X1) / mu;
  const C1 = 5.5, C2 = C1 - (qTrue * integral * 1000) / V;
  const r = rateFromCourse({ times: [0, dt], conc: [C1, C2], biomass: [X1, X2], volume: V });
  near(r.rate, qTrue, 1e-9, "uptake recovered from concentrations");
  ok(r.flags.includes("twoPoints"), "flagged as a two-point estimate");
  const s = rateFromCourse({ times: [0, dt], conc: [0, (0.8 * integral * 1000) / V], biomass: [X1, X2], volume: V });
  near(s.rate, -0.8, 1e-9, "secretion is negative");
  near(growthRate([0, dt], [X1, X2]).mu, mu, 1e-12, "growth rate from the biomass series");
}
{ // evaporation: what counts is the amount in the well, not the concentration
  const dt = 24, X1 = 1e-4, X2 = 2e-4, mu = Math.log(2) / dt, integral = (X2 - X1) / mu;
  const qTrue = 0.5, V = [2, 1.8]; // 10 % lost to evaporation
  const n1 = 5.5 * V[0], n2 = n1 - qTrue * integral * 1000;
  const r = rateFromCourse({ times: [0, dt], conc: [5.5, n2 / V[1]], biomass: [X1, X2], volume: V });
  near(r.rate, qTrue, 1e-9, "a shrinking volume is handled through the amount");
  const ignoring = rateFromCourse({ times: [0, dt], conc: [5.5, n2 / V[1]], biomass: [X1, X2], volume: 2 });
  console.log(`     taking the volume as constant would give ${ignoring.rate.toPrecision(3)} instead of ${qTrue} (${(100 * (ignoring.rate / qTrue - 1)).toFixed(0)} %)`);
}

// ---------------------------------------------------------------- 3. cell-free control ----
console.log("\n== cell-free control (first-order decay)");
{
  const V = 2, dt = 48, k = 0.004, mu = Math.log(4) / dt, X1 = 1e-4;
  const X = (t) => X1 * Math.exp(mu * t);
  const integral = (X(dt) - X1) / mu, qTrue = 0.08;
  const N = 9600, h = dt / N, keep = [0, N / 2, N];
  let C = 2.0; const times = [], conc = [];
  for (let i = 0; i <= N; i++) {
    const t = i * h;
    if (keep.includes(i)) { times.push(t); conc.push(C); }
    C += h * (-k * C - (qTrue * X(t) * 1000) / V);
  }
  const blank = { times: [0, dt], conc: [2.0, 2.0 * Math.exp(-k * dt)], volume: V };
  near(decayConstant([blank]).k, k, 1e-9, "k recovered from the control");
  const rep = { times, conc, biomass: times.map(X), volume: V };
  const got = ratesFromReplicates({ replicates: [rep], blanks: [blank] });
  near(got.rate, qTrue, qTrue * 0.01, "uptake recovered with the first-order correction (within 1 %)");
  const naive = ((conc[0] * V - conc[2] * V) / 1000 - (blank.conc[0] * V - blank.conc[1] * V) / 1000) / integral;
  console.log(`     subtracting the control's own loss gives ${naive.toPrecision(3)} instead of ${qTrue} (${(100 * (naive / qTrue - 1)).toFixed(1)} %): it over-corrects`);
  const noBlank = ratesFromReplicates({ replicates: [rep], blanks: [] });
  ok(noBlank.rate > qTrue && noBlank.flags.includes("noBlank"), `with no control the uptake is overestimated: ${noBlank.rate.toPrecision(3)} vs ${qTrue}`);
  ok(Math.abs(got.rate - qTrue) < Math.abs(naive - qTrue), "the first-order correction beats subtracting the control");
  // evaporation in the control: fitting ln C instead of ln(C·V) would hide part of the decay
  const vol2 = [V, V * 0.9];
  const evap = { times: [0, dt], conc: blank.conc.map((c, i) => (c * V) / vol2[i]), volume: vol2 };
  near(decayConstant([evap]).k, k, 1e-9, "k is right even when the control evaporates");
  const lnCOnly = Math.max(0, -(Math.log(evap.conc[1]) - Math.log(evap.conc[0])) / dt);
  console.log(`     fitting ln C instead of ln(C·V) would give k = ${lnCOnly.toPrecision(3)} instead of ${k}: evaporation masks the decay`);
}

// ---------------------------------------------------------------- 4. replicates ----
console.log("\n== replicates (a seeding difference must not become variance)");
{
  const dt = 24, mu = Math.log(2) / dt, qTrue = 0.4, V = 2;
  const wells = (X1, q) => { const X2 = X1 * 2, I = (X2 - X1) / mu;
    return { times: [0, dt], conc: [5.5, 5.5 - (q * I * 1000) / V], biomass: [X1, X2], volume: V }; };
  const reps = [wells(1e-4, qTrue), wells(1.4e-4, qTrue)]; // seeded 40 % apart, same specific rate
  const got = ratesFromReplicates({ replicates: reps, blanks: [] });
  near(got.rate, qTrue, 1e-9, "rate per replicate, then averaged");
  ok(got.sd < 1e-9, `a 40 % difference in seeding gives sd = ${got.sd.toExponential(2)}: it cancels, as it must`);
  const real = ratesFromReplicates({ replicates: [wells(1e-4, 0.4), wells(1.4e-4, 0.5)], blanks: [] });
  ok(real.sd > 0.01, `a genuine difference between wells does show: sd = ${real.sd.toPrecision(2)}`);
  const mixed = [reps[0], { times: [0, 12, dt], conc: [5.5, 5.2, reps[1].conc[1]], biomass: [1.4e-4, 1.4e-4 * Math.SQRT2, 2.8e-4], volume: V }];
  ok(ratesFromReplicates({ replicates: mixed, blanks: [] }).n === 2, "replicates sampled at different times are both kept");
  near(toGdw(1e6, "cells"), 1e6 * 300e-12, 1e-18, "1e6 cells -> gDW at 300 pg/cell");
  near(toGdw(100, "ugProtein"), 100e-6 / 0.59, 1e-12, "100 µg protein -> gDW at 0.59 g protein/gDW");
  near(toGdw(1e8, "cells", { gdwPerCell: 2.5e-13 }), 2.5e-5, 1e-18, "the conversion is configurable (yeast, E. coli)");
}

// ---------------------------------------------------------------- 5. the input table ----
console.log("\n== input table");
{
  const V = 2, dt = 24, c1 = 1e6, c2 = 2e6;
  const g1 = toGdw(c1, "cells"), g2 = toGdw(c2, "cells");
  const I = biomassIntegral([0, dt], [g1, g2]);
  const qGlc = 0.4, qLac = -0.7;
  const conc = (q, c0) => c0 - (q * I * 1000) / V;
  const tsv = [
    "condition\treplicate\ttime_h\tmetabolite\tconc_mM\tvolume_mL\tcells",
    `control\t1\t0\tglucose\t5.5\t${V}\t${c1}`,
    `control\t1\t${dt}\tglucose\t${conc(qGlc, 5.5)}\t${V}\t${c2}`,
    `control\t2\t0\tglucose\t5.5\t${V}\t${c1}`,
    `control\t2\t${dt}\tglucose\t${conc(qGlc, 5.5)}\t${V}\t${c2}`,
    `control\t1\t0\tL-lactate\t0\t${V}\t${c1}`,
    `control\t1\t${dt}\tL-lactate\t${conc(qLac, 0)}\t${V}\t${c2}`,
    `blank\t1\t0\tglucose\t5.5\t${V}\t0`,
    `blank\t1\t${dt}\tglucose\t5.45\t${V}\t0`,
    `blank\t2\t0\tglucose\t5.5\t${V}\t0`,
    `blank\t2\t${dt}\tglucose\t5.46\t${V}\t0`,
    `control\t1\t0\tstray\t1.0\t${V}\t${c1}`,
    `control\t1\t${dt}\t\t1.0\t${V}\t${c2}`,
    `control\t1\t${dt}\tnegative\t-1\t${V}\t${c2}`,
    `control\t3\t0\tglucose\t5.5\t${V}\t${c1}`,
    `control\t3\t0\tglucose\t5.5\t${V}\t${c1}`,
    `control\t3\t${dt}\tglucose\t${conc(qGlc, 5.5)}\t${V}\t${c2}`,
  ].join("\n");
  const parsed = parseExoTable(tsv);
  ok(parsed.errors.length === 0, `parses with no errors (biomass unit: ${parsed.biomassUnit})`);
  ok(parsed.warnings.some((w) => /stray/.test(w)), "one time point only -> warning");
  ok(parsed.warnings.some((w) => /no metabolite name/.test(w)), "a row with no metabolite name -> warning");
  ok(parsed.warnings.some((w) => /negative/.test(w)), "a negative concentration -> warning");
  ok(parsed.warnings.some((w) => /two rows at t = 0/.test(w)), "a repeated time point is averaged, with a warning");
  const groups = groupForRates(parsed);
  const glc = groups.find((g) => g.metabolite === "glucose");
  ok(glc.replicates.length === 3 && glc.blanks.length === 2, `glucose: ${glc.replicates.length} replicates, ${glc.blanks.length} controls pooled`);
  const r = ratesFromReplicates(glc);
  ok(r.k > 0, `the controls give k = ${r.k.toPrecision(3)} 1/h`);
  near(r.rate, qGlc, qGlc * 0.05, "glucose uptake recovered from the table");
  near(ratesFromReplicates(groups.find((g) => g.metabolite === "L-lactate")).rate, qLac, 1e-6, "lactate comes out as secretion");
  ok(!groups.some((g) => g.metabolite === ""), "the nameless row did not create a series");
  const noBio = parseExoTable("time_h\tmetabolite\tconc_mM\n0\tglucose\t5.5\n24\tglucose\t3.0");
  ok(noBio.errors.length === 1 && /biomass/.test(noBio.errors[0]), `a missing biomass column is an error, not a crash: "${noBio.errors[0]}"`);
  const mg = parseExoTable("time_h\tmetabolite\tconc_mM\tprotein_mg\n0\tglucose\t5.5\t0.1\n24\tglucose\t3\t0.2");
  ok(mg.errors.length === 0 && mg.biomassUnit === "mgProtein", "protein in mg is accepted too");
}

// ---------------------------------------------------------------- 6. the consistency fit ----
console.log("\n== fit of measured rates to Human-GEM");
const model = JSON.parse(fs.readFileSync(new URL("../data/human-gem/model.json", import.meta.url)));
const highs = await loadHighs({ print: () => {}, printErr: () => {} });
const fitter = new ExchangeFit(highs, packNetwork(model));
const ex = findExchanges(model);
const name = (m) => model.bases[model.mets[m].b].n;
const growthRxn = model.rxns.findIndex((r) => r.id === "MFL_biomass_nocof");
const medium = presetMedium(model, "defaultSerum");

function setup(set, mu, { sigmaRel = 0.1, muSigma = 0.1, trace = 0.001, scale = {} } = {}) {
  const items = measuredPreset(model, ex, set);
  const measurements = items.map((x) => ({ rxn: ex.get(x.met).rxn, dir: ex.get(x.met).dir, met: x.met, name: name(x.met),
    rate: x.rate * (scale[name(x.met)] ?? 1), sigma: Math.max(1e-4, Math.abs(x.rate) * sigmaRel) }));
  if (mu !== null) measurements.push({ rxn: growthRxn, dir: 1, met: -1, name: "growth", rate: mu, sigma: Math.abs(mu) * muSigma });
  const base = applyMedium(model, ex, medium, model.rxns.map((r) => r.lb), model.rxns.map((r) => r.ub));
  const { lb, ub } = fitBounds(model, ex, base.lb, base.ub,
    { measuredMets: new Set(items.map((x) => x.met)), presentMets: new Set(medium.filter((q) => q.unlimited || q.conc > 0).map((q) => q.met)),
      trace, openRxns: [growthRxn] });
  return { lb, ub, measurements };
}
const worstOf = (res, m) => m.find((x) => x.rxn === res.residuals[0].rxn)?.name;

{ // the trap this used to fall into: the medium's generic ceilings
  const it = medium.find((q) => name(q.met) === "glucose");
  const cap = uptakeCapacity(it.conc, it.vmax, it.km);
  const items = measuredPreset(model, ex, "hek293");
  const measuredGlc = items.find((x) => name(x.met) === "glucose").rate;
  ok(measuredGlc > cap, `the standard medium caps glucose at ${cap.toFixed(3)} while the HEK293 cells take up ${measuredGlc.toFixed(3)}: the ceiling has to be freed`);
  const base = applyMedium(model, ex, medium, model.rxns.map((r) => r.lb), model.rxns.map((r) => r.ub));
  const freed = fitBounds(model, ex, base.lb, base.ub, { measuredMets: new Set(items.map((x) => x.met)), openRxns: [growthRxn] });
  const e = ex.get(it.met);
  ok(freed.ub[e.rxn] >= 1000 && freed.lb[e.rxn] <= -1000, "fitBounds frees every measured exchange in both directions");
  // and unmeasured organic secretion is capped, as sim.js does
  const sorb = [...ex].find(([m]) => name(m) === "sorbitol");
  if (sorb) ok(Math.abs(freed.lb[sorb[1].rxn]) <= 0.001 + 1e-12 || Math.abs(freed.ub[sorb[1].rxn]) <= 0.001 + 1e-12,
    "an unmeasured organic metabolite may only leave in traces");
}

{ // HEK293: 22 exchanges including O2, cells grew at 0.028 1/h (docs, section 5)
  const s = setup("hek293", 0.028);
  const res = fitter.fit(s);
  ok(res.status === "optimal", `HEK293 set: status ${res.status}, weighted deviation ${res.deviation.toPrecision(4)}`);
  console.log(`     growth fitted ${res.v[growthRxn].toPrecision(3)} 1/h (measured 0.028); worst residuals: ${res.residuals.slice(0, 3).map((r) => `${s.measurements.find((x) => x.rxn === r.rxn)?.name} ${r.residual >= 0 ? "+" : ""}${r.residual.toPrecision(3)}`).join(", ")}`);
  ok(res.v[growthRxn] > 0.01, "the fit keeps the cell growing (growth is one of the measurements)");

  const bal = exchangeBalance(model, ex, res.v);
  ok(bal.rq > 0.5 && bal.rq < 2, `the respiratory quotient is physiological: CO2/O2 = ${bal.rq.toFixed(2)}`);
  console.log(`     carbon closure out/in = ${bal.closure.C.toFixed(2)} (the rest goes into biomass); formulas not counted: ${bal.unknown.length}`);

  // A rate that is five times too high is matched EXACTLY and the deviation does not move: the surplus
  // carbon leaves as CO2 without the oxygen to make it. The deviation alone cannot see this; the
  // respiratory quotient can, which is why the fit must always be read together with the balance.
  const broken = setup("hek293", 0.028, { scale: { glucose: 5 } });
  const res2 = fitter.fit(broken);
  const glc = broken.measurements.find((x) => x.name === "glucose");
  near(glc.dir * res2.v[glc.rxn], glc.rate, 1e-3, "a 5x glucose uptake is reproduced exactly by the fit");
  ok(Math.abs(res2.deviation - res.deviation) < 1e-3, `so the weighted deviation does not change (${res2.deviation.toPrecision(4)}): it is blind to this error`);
  const bal2 = exchangeBalance(model, ex, res2.v);
  ok(bal2.rq > 5, `but the respiratory quotient gives it away: CO2/O2 = ${bal2.rq.toFixed(1)} against ${bal.rq.toFixed(2)} (O2 pinned at ${bal2.o2.toFixed(3)} while CO2 reaches ${bal2.co2.toFixed(2)})`);
  ok(bal2.warnings.length >= 2 && bal.warnings.length === 0, `the balance raises ${bal2.warnings.length} warnings here and none on the good fit:\n       - ${bal2.warnings.join("\n       - ")}`);
  // the electron budget: the surplus is exactly the proton efflux, which is the charge defect of the model
  const extraGap = bal2.electrons.gap - bal.electrons.gap;
  near(extraGap, -bal2.hplus, Math.abs(bal2.hplus) * 0.05, "the extra electron gap equals the proton efflux (the charge-unbalanced cycle)");
  console.log(`     electron gap ${bal.electrons.gap.toFixed(1)} -> ${bal2.electrons.gap.toFixed(1)} (the baseline gap is the carbon kept as biomass); H+ efflux ${(-bal2.hplus).toFixed(1)}`);

  const noG = fitter.fit(setup("hek293", null));
  console.log(`     without the growth measurement: growth = ${noG.v[growthRxn].toPrecision(3)} 1/h, deviation ${noG.deviation.toPrecision(4)}`);
  ok(noG.deviation <= res.deviation + 1e-6, "dropping a measurement cannot worsen the fit (it is a relaxation)");

  const open = fitter.fit(setup("hek293", 0.028, { trace: Infinity }));
  ok(open.deviation <= res.deviation + 1e-6, `leaving unmeasured secretion free drops the deviation ${res.deviation.toPrecision(4)} -> ${open.deviation.toPrecision(4)}: that gap is what an unmeasured product would explain`);
}

{ // an impossible demand: more growth than the measured nutrients can support
  const s = setup("hek293", 0.40, { muSigma: 0.05 });
  const res = fitter.fit(s);
  ok(res.status === "optimal", `an impossible growth returns a diagnosis, not an infeasibility: ${res.status}`);
  console.log(`     growth demanded 0.40 1/h -> fitted ${res.v[growthRxn].toPrecision(3)}; worst residual: ${worstOf(res, s.measurements)} (${res.residuals[0].residual.toPrecision(3)})`);
  ok(res.residuals[0].weighted > 0, "the conflict is reported on specific measurements");
}

{ // a metabolite measured only as a product must not be freed for uptake
  const items = measuredPreset(model, ex, "hek293");
  const base = applyMedium(model, ex, medium, model.rxns.map((r) => r.lb), model.rxns.map((r) => r.ub));
  const lac = [...ex].find(([m]) => name(m) === "L-lactate");
  const inMedium = new Set(medium.filter((q) => q.unlimited || q.conc > 0).map((q) => q.met));
  ok(!inMedium.has(lac[0]), "lactate is not in the standard medium (it is only secreted)");
  const strict = fitBounds(model, ex, base.lb, base.ub,
    { measuredMets: new Set(items.map((x) => x.met)), presentMets: inMedium, openRxns: [growthRxn] });
  const upBound = lac[1].dir < 0 ? -strict.lb[lac[1].rxn] : strict.ub[lac[1].rxn];
  ok(upBound <= 0, `its uptake stays shut (max uptake ${upBound}) while secretion is free`);
  const loose = fitBounds(model, ex, base.lb, base.ub, { measuredMets: new Set(items.map((x) => x.met)), openRxns: [growthRxn] });
  ok((lac[1].dir < 0 ? -loose.lb[lac[1].rxn] : loose.ub[lac[1].rxn]) >= 1000, "without presentMets it would be freed in both directions");
}

{ // reproducibility: stage 2 must pick the same optimum every time
  const s = setup("hek293", 0.028);
  const a = fitter.fit(s), b = fitter.fit(s);
  ok(a.residuals.every((r, i) => Math.abs(r.residual - b.residuals[i].residual) < 1e-9), "the fit is reproducible run to run");
}

// ---------------------------------------------------------------- 7. the dialog's own sequence ----
// The exact steps the importer performs: paste -> parse -> rates -> resolve names -> fit.
console.log("\n== the importer, end to end");
{
  const V = 2, dt = 24, c1 = 1e6, c2 = 2.2e6;
  const DATA = { glucose: [5.5, 3.9], glutamine: [2.0, 1.45], "L-lactate": [0, 2.6], alanine: [0, 0.22], HMDB0000148: [0.1, 0.33] };
  const rows = ["condition\treplicate\ttime_h\tmetabolite\tconc_mM\tvolume_mL\tcells"];
  for (const rep of [1, 2]) for (const [met, [a, b]] of Object.entries(DATA)) {
    const jitter = rep === 2 ? 1.04 : 1;
    rows.push(`control\t${rep}\t0\t${met}\t${a}\t${V}\t${c1}`);
    rows.push(`control\t${rep}\t${dt}\t${met}\t${(a + (b - a) * jitter).toFixed(4)}\t${V}\t${c2}`);
  }
  rows.push(`blank\t1\t0\tglutamine\t2.0\t${V}\t0`, `blank\t1\t${dt}\tglutamine\t1.90\t${V}\t0`);
  const parsed = parseExoTable(rows.join("\n"), { gdwPerCell: 300e-12 });
  ok(parsed.errors.length === 0 && parsed.warnings.length === 0, `a realistic table parses clean (biomass as ${parsed.biomassUnit})`);
  const groups = groupForRates(parsed);
  // the importer's identifier lookup
  const index = new Map();
  for (const [met] of ex) { const b = model.bases[model.mets[met].b];
    for (const q of [b.n, b.id, b.kegg, b.hmdb, b.chebi, b.bigg, b.mnx]) if (q) index.set(String(q).toLowerCase(), met); }
  const rates = groups.map((g) => { const r = ratesFromReplicates(g); return { ...r, label: g.metabolite, met: index.get(g.metabolite.toLowerCase()) }; })
    .filter((r) => Number.isFinite(r.rate));
  ok(rates.every((r) => r.met !== undefined), "every metabolite is found, including the one given as an HMDB id");
  ok(name(rates.find((r) => r.label === "HMDB0000148").met) === "glutamate", "HMDB0000148 resolves to glutamate");
  const glc = rates.find((r) => r.label === "glucose"), lac = rates.find((r) => r.label === "L-lactate");
  ok(glc.rate > 0 && lac.rate < 0, `signs are right: glucose ${glc.rate.toFixed(3)} (uptake), lactate ${lac.rate.toFixed(3)} (secretion)`);
  ok(!rates.find((r) => r.label === "glutamine").flags.includes("noBlank") && glc.flags.includes("noBlank"),
    "the control is used for glutamine only, and glucose is flagged as having none");
  const mu = groups[0].replicates.map((r) => growthRate(r.times, r.biomass).mu).reduce((a, b) => a + b, 0) / groups[0].replicates.length;
  const measurements = rates.map((r) => ({ rxn: ex.get(r.met).rxn, dir: ex.get(r.met).dir, met: r.met, rate: r.rate,
    sigma: r.sd > 0 ? r.sd : Math.max(1e-4, Math.abs(r.rate) * 0.1) }));
  measurements.push({ rxn: growthRxn, dir: 1, met: -1, rate: mu, sigma: Math.abs(mu) * 0.1 });
  const b0 = applyMedium(model, ex, medium, model.rxns.map((r) => r.lb), model.rxns.map((r) => r.ub));
  const bb = fitBounds(model, ex, b0.lb, b0.ub, { measuredMets: new Set(rates.map((r) => r.met)),
    presentMets: new Set(medium.filter((q) => q.unlimited || q.conc > 0).map((q) => q.met)), openRxns: [growthRxn] });
  const res = fitter.fit({ lb: bb.lb, ub: bb.ub, measurements });
  const bal = exchangeBalance(model, ex, res.v);
  ok(res.status === "optimal" && res.deviation < 1e-6, `the rates are consistent: deviation ${res.deviation.toPrecision(3)}`);
  near(res.v[growthRxn], mu, mu * 1e-3, "and the growth measured from the cell counts is reproduced");
  ok(bal.warnings.length === 0, `no balance warning; RQ ${bal.rq.toFixed(2)}, carbon out/in ${Math.round(100 * bal.closure.C)} %`);
}

// ------------------------------------------------- 8. units against published numbers ----
// The likeliest silent error in this pipeline is a factor of 1000 or of a cell mass. This anchors it to
// three independently published numbers for HEK293 (Dietmair et al., PLoS ONE 7:e43394, 2012): the glucose
// uptake of 528 µmol/gDW/h and the lactate secretion of 695 (Table 2), the dry mass of 514 pg/cell
// (Table S1) and the growth rate of 0.028 1/h. A 48 h course that those rates would produce is built and
// run back through the importer, which must return the published rates.
console.log("\n== units, against the published HEK293 numbers");
{
  const Q_GLC = 0.528, Q_LAC = -0.695, GDW_CELL = 514e-12, MU = 0.028, dt = 48, V = 20; // mL, a shake flask
  const n1 = 5e5, n2 = n1 * Math.exp(MU * dt);           // cells/mL -> the flask holds n·V
  const X = (n) => n * V * GDW_CELL;                      // gDW in the flask
  const integral = (X(n2) - X(n1)) / MU;                  // exact for exponential growth
  const conc = (q, c0) => c0 - (q * integral * 1000) / V; // mM
  const tsv = ["condition\treplicate\ttime_h\tmetabolite\tconc_mM\tvolume_mL\tcells",
    `hek293\t1\t0\tglucose\t25\t${V}\t${n1 * V}`,
    `hek293\t1\t${dt}\tglucose\t${conc(Q_GLC, 25)}\t${V}\t${n2 * V}`,
    `hek293\t1\t0\tL-lactate\t0\t${V}\t${n1 * V}`,
    `hek293\t1\t${dt}\tL-lactate\t${conc(Q_LAC, 0)}\t${V}\t${n2 * V}`].join("\n");
  const parsed = parseExoTable(tsv, { gdwPerCell: GDW_CELL });
  ok(parsed.errors.length === 0, "the table parses with the published cell mass of 514 pg/cell");
  const got = groupForRates(parsed).map((g) => ({ name: g.metabolite, rate: ratesFromReplicates(g).rate }));
  const glc = got.find((x) => x.name === "glucose").rate, lac = got.find((x) => x.name === "L-lactate").rate;
  near(glc, Q_GLC, Q_GLC * 1e-9, "the published glucose uptake comes back in mmol/gDW/h");
  near(lac, Q_LAC, Math.abs(Q_LAC) * 1e-9, "and the published lactate secretion, with its sign");
  near(growthRate([0, dt], [X(n1), X(n2)]).mu, MU, 1e-12, "and the published growth rate of 0.028 1/h");
  console.log(`     glucose fell from 25 to ${conc(Q_GLC, 25).toFixed(2)} mM over ${dt} h: ${glc.toPrecision(3)} mmol/gDW/h, the published 0.528`);
  // the error a wrong cell mass would cause, which is what this test exists to catch
  const wrong = groupForRates(parseExoTable(tsv, { gdwPerCell: 300e-12 })).find((g) => g.metabolite === "glucose");
  console.log(`     read with 300 pg/cell instead it would give ${ratesFromReplicates(wrong).rate.toPrecision(3)} (${(100 * (514 / 300 - 1)).toFixed(0)} % high): the conversion has to come from the paper`);
}

console.log(bad ? `\n${bad} check(s) FAILED` : "\nall passed");
process.exit(bad ? 1 : 0);
