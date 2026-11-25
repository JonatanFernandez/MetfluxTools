// Cross-validation of the unlabelled protein pool (liver_pools.mjs).
//
// The prediction for a given turnover rate F depends on F alone -- never on which donors or which metabolites
// were used to choose it. So every held-out test can be scored from the predictions liver_pools.mjs already
// wrote, without solving anything again. Three tests, in increasing order of how much they prove:
//
//   1. Leave one donor out. Choose F on four donors, score the fifth. Shows F is not fitted to one donor's
//      noise. Weak on its own: the five donors are biological replicates of one experiment.
//   2. Leave metabolites out. Choose F on four amino acids and predict five others -- including the urea-cycle
//      intermediates -- without refitting. This is the real test of the mechanism: one number, a protein
//      composition that was never fitted, and metabolites the number never saw.
//   3. No fitting at all. Take F from the literature range of liver protein turnover and score everything.
//      Zero free parameters. If this works the result does not depend on an optimiser at all.
//
// Usage: node tools/validation/liver_pools_cv.mjs
import fs from "node:fs";
const R = new URL("../../", import.meta.url).pathname;
const runs = JSON.parse(fs.readFileSync(R + "tools/validation/liver_pools.json"));

// Glucose and citrate are excluded throughout: their error has the opposite sign and is a flux problem, not a
// dilution one. liver_pools.mjs scores them separately as the control.
const CONTROL = new Set(["glucose", "citrate"]);
const metNames = runs[0].rows.map((r) => r.name).filter((n) => !CONTROL.has(n));

const mae = (a, b) => a.reduce((s, x, i) => s + Math.abs(x - b[i]), 0) / a.length;
/** Mean |error| over the isotopologue fractions of `mets`, measured on `donors`. */
function score(run, mets, donors) {
  const m = [], p = [];
  for (const row of run.rows) {
    if (!mets.includes(row.name)) continue;
    const obs = row.perDonor[0].map((_, k) =>
      donors.reduce((s, d) => s + row.perDonor[d][k], 0) / donors.length);
    obs.forEach((x, k) => { m.push(x); p.push(row.pred[k]); });
  }
  return 100 * mae(m, p);
}
const pick = (mets, donors) => runs.reduce((a, b) => (score(b, mets, donors) < score(a, mets, donors) ? b : a));
const f1 = (x) => x.toFixed(1).padStart(5);
const ALL = [0, 1, 2, 3, 4];
const noPool = runs.find((r) => r.F === 0);

console.log(`${runs.length} turnover rates, ${metNames.length} dilution-limited metabolites, 5 donors`);
console.log(`no pool at all: ${f1(score(noPool, metNames, ALL))} points\n`);

console.log("1. LEAVE ONE DONOR OUT");
console.log("   held-out donor   F chosen on the other four   error on the held-out donor   no pool");
let heldSum = 0, baseSum = 0;
for (const d of ALL) {
  const train = ALL.filter((x) => x !== d);
  const best = pick(metNames, train);
  const held = score(best, metNames, [d]), base = score(noPool, metNames, [d]);
  heldSum += held; baseSum += base;
  console.log(`   D${d + 1}               ${best.pctPerDay.toFixed(1).padStart(8)} %/day            ${f1(held)}                  ${f1(base)}`);
}
console.log(`   mean over held-out donors: ${f1(heldSum / 5)} against ${f1(baseSum / 5)} with no pool`);

console.log("\n2. LEAVE METABOLITES OUT  (choose F on four amino acids, predict five others, no refitting)");
const TRAIN_MET = ["alanine", "glutamate", "glutamine", "glycine"].filter((n) => metNames.includes(n));
const TEST_MET = metNames.filter((n) => !TRAIN_MET.includes(n));
const bestMet = pick(TRAIN_MET, ALL);
console.log(`   chosen on   ${TRAIN_MET.join(", ")}`);
console.log(`   predicting  ${TEST_MET.join(", ")}`);
console.log(`   F = ${bestMet.pctPerDay.toFixed(1)} %/day`);
console.log(`   on the metabolites it never saw: ${f1(score(bestMet, TEST_MET, ALL))} points, against ${f1(score(noPool, TEST_MET, ALL))} with no pool`);
for (const n of TEST_MET) {
  const row = bestMet.rows.find((r) => r.name === n), base = noPool.rows.find((r) => r.name === n);
  const obs = row.perDonor[0].map((_, k) => ALL.reduce((s, d) => s + row.perDonor[d][k], 0) / 5);
  console.log(`     ${n.padEnd(12)} M+0 measured ${(100 * obs[0]).toFixed(0).padStart(3)} %   predicted ${(100 * row.pred[0]).toFixed(0).padStart(3)} %   (was ${(100 * base.pred[0]).toFixed(0)} % with no pool)`);
}

console.log("\n3. NO FITTING  (F straight from the literature range of liver protein turnover, 3-10 %/day)");
for (const r of runs.filter((x) => x.pctPerDay >= 3 && x.pctPerDay <= 10)) {
  console.log(`   ${r.pctPerDay.toFixed(1).padStart(5)} %/day : ${f1(score(r, metNames, ALL))} points   (no pool ${f1(score(noPool, metNames, ALL))})`);
}
console.log("\n   Zero free parameters: the rate is read off the literature, the composition off the model's own");
console.log("   biomass, and nothing is chosen to fit these measurements.");
