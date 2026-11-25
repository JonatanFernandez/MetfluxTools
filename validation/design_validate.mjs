// Validation Suite 5: Tracer Experiment Designer validation against classic 13C literature designs
// Uses the exact same scoring algorithm as js/ui/design.js:
//   - IsotopologueSampler + whole-cell production-weighted MIDs + Euclidean separation / sd (1 mol%)
//   - Ranks all available commercial tracers from data/tracers.json across candidate measurable metabolites
// Literature benchmarks:
//   1. Glycolysis vs. oxidative PPP (low vs high G6PDH flux; Lee et al. 1998, Am J Physiol):
//      Does [1,2-13C2]glucose rank top for distinguishing lactate/pyruvate/R5P?
//   2. Reductive carboxylation vs. oxidative TCA (ETC/Complex I block or hypoxia; Metallo et al. 2012, Nature;
//      Mullen et al. 2012, Nature): Do [U-13C5]glutamine / [1-13C]glutamine rank top for citrate/malate/aspartate?
//   3. SDHB knockout vs. WT (Cardaci et al. 2015, Nature Cell Biol): Top distinguishing tracers & metabolites.
//
// Run: node tools/validation/design_validate.mjs
import fs from "node:fs";
const R = new URL("../../", import.meta.url).pathname;
const { default: loadHighs } = await import(R + "vendor/highs/highs.mjs");
const { FbaProblem } = await import(R + "js/core/fba.js");
const { presetMedium } = await import(R + "js/core/medium.js");
const { Simulator } = await import(R + "js/sim.js");
const { ModelIndex } = await import(R + "js/modelindex.js");
const { IsotopologueSampler } = await import(R + "js/core/isotopologues.js");
const { tracerLabel } = await import(R + "js/core/tracerlib.js");
const { gprActivity } = await import(R + "js/core/gpr.js");

const rd = (f) => JSON.parse(fs.readFileSync(R + f, "utf8"));
const model = rd("data/human-gem/model.json");
const st = rd("data/human-gem/structures.json");
const am = rd("data/human-gem/atommaps.json");
const layout = rd("data/human-gem/layout.json");
const lib = rd("data/tracers.json");
const index = new ModelIndex(model, layout);

const highs = await loadHighs({ print: () => {}, printErr: () => {} });
let pb;
const backend = async (m) => {
  if (m.type === "init") { pb = new FbaProblem(highs, m.net); return true; }
  const r = pb.solve(m.params); r.ms = 0; return r;
};
const sim = new Simulator(model, st, am, backend);
const sup = presetMedium(model, "defaultSerum");
const objIdx = model.rxns.findIndex((r) => r.id === "MFL_biomass_nocof");

const MEASURED_NAMES = ["L-lactate", "pyruvate", "alanine", "serine", "glycine", "citrate", "AKG", "succinate", "fumarate", "malate", "aspartate", "glutamate", "glutamine", "ribose-5-phosphate"];
const nameToBase = new Map(model.bases.map((b, i) => [b.n.toLowerCase(), i]));
const bases = MEASURED_NAMES.map((n) => [n, nameToBase.get(n.toLowerCase())]).filter(([, b]) => b !== undefined);

const inMedium = new Set(sup.filter((it) => sim.exchanges.has(it.met)).map((it) => model.bases[model.mets[it.met].b].n));
// Score all 13C tracers whose substrate is in the medium
const candidateTracers = lib.tracers.filter((tr) => inMedium.has(tr.met) && tr.label?.C);

function wholeCellMid(sampler, v, base) {
  const mets = index.metsOfBase[base].filter((m) => !model.mets[m].removed);
  const acc = []; let wsum = 0;
  for (const m of mets) {
    let prod = 0;
    for (let p = index.rxnStart[m]; p < index.rxnStart[m + 1]; p++) {
      const x = v[index.rxnList[p]] * index.rxnCoef[p];
      if (x > 0) prod += x;
    }
    if (!(prod > 1e-9)) continue;
    const r = sampler.mid(m, 1200, 120);
    if (!r) continue;
    r.dist.forEach((p, k) => { acc[k] = (acc[k] || 0) + prod * p; });
    wsum += prod;
  }
  return wsum > 0 ? acc.map((x) => (x || 0) / wsum) : null;
}
const euclid = (p, q) => { let s = 0; for (let k = 0; k < Math.max(p.length, q.length); k++) { const d = (p[k] || 0) - (q[k] || 0); s += d * d; } return Math.sqrt(s); };

function rankTracers(sits, tracers, sd = 0.01) {
  const rows = [];
  for (const tr of tracers) {
    const base = model.bases.findIndex((b) => b.n === tr.met);
    if (base < 0) continue;
    const lab = tracerLabel(tr, st[model.bases[base].id]);
    if (!lab) continue;
    const per = sits.map((sit) => {
      const labels = new Map();
      for (const it of sit.sup) if (model.mets[it.met].b === base) labels.set(it.met, lab.label);
      const tr2 = sim.traceFlux(sit.v, sit.src, labels);
      const sampler = new IsotopologueSampler(sim.tracer, index, sit.v, tr2.inputs, tr2.e, lab.elements);
      const m = new Map();
      for (const [, b] of bases) m.set(b, wholeCellMid(sampler, sit.v, b));
      return m;
    });
    const cells = bases.map(([n, b]) => {
      const p = per[0].get(b), q = per[1].get(b);
      return { n, s: p && q ? euclid(p, q) / sd : null };
    });
    const top3 = cells.filter((c) => c.s !== null).sort((a, b) => b.s - a.s).slice(0, 3);
    const total = top3.reduce((s, x) => s + x.s, 0);
    rows.push({ id: tr.id, name: tr.name, met: tr.met, total: +total.toFixed(1), top3: top3.map((c) => ({ met: c.n, sep: +c.s.toFixed(1) })) });
  }
  rows.sort((a, b) => b.total - a.total);
  return rows;
}

const results = [];

const { measuredPreset, measuredNgam } = await import(R + "js/core/measured.js");
const meas = { items: measuredPreset(model, sim.exchanges, "nci60").map((x) => ({ met: x.met, rate: x.rate })), tol: 0.2, ngam: measuredNgam("nci60") };

// Case 1: Oxidative PPP vs Glycolysis (Mode B FVA on cytosolic 6-phosphogluconate dehydrogenase PGD / MAR04474; Lee et al. 1998)
console.log("== Case 1: Glycolysis vs Oxidative PPP (Mode B FVA on cytosolic PGD / MAR04474; Lee et al. 1998) ==");
{
  const pgdIdx = model.rxns.findIndex((r) => r.id === "MAR04474");
  const glcTracers = candidateTracers.filter((tr) => tr.met === "glucose");
  const cfgLow = { objective: objIdx, maximize: true, parsimonious: true, measured: meas };
  const fLow = (await sim.runSteady(cfgLow, sup, st)).frames[0];
  const fvaPgd = await sim.fluxRanges(cfgLow, sup, [pgdIdx], fLow.v, { frac: 0.9, solutions: true, cycles: 0.05 });
  const [vmin, vmax] = fvaPgd.solutions[0];
  console.log(`  PGD rxn=${model.rxns[pgdIdx].id} (${model.rxns[pgdIdx].n}): FVA=[${fvaPgd.ranges[0][0].toFixed(4)}, ${fvaPgd.ranges[0][1].toFixed(4)}]`);
  const ranked = rankTracers([{ name: "Low oxPPP", v: vmin, src: fLow, sup }, { name: "High oxPPP", v: vmax, src: fLow, sup }], glcTracers);
  const rank12 = ranked.findIndex((x) => x.id === "glc_12") + 1;
  const rank1 = ranked.findIndex((x) => x.id === "glc_1") + 1;
  console.log(`  Top 5 glucose tracers (${glcTracers.length} tested):`, ranked.slice(0, 5).map((x, i) => `#${i+1} ${x.name}(${x.total}) [${x.top3.map((c)=>`${c.met}:${c.sep}x`).join(",")}]`).join(" | "));
  console.log(`  [1,2-13C2]glucose rank=#${rank12}/${ranked.length}, [1-13C]glucose rank=#${rank1}/${ranked.length}`);
  results.push({ case: "Glycolysis vs oxPPP (Mode B PGD/MAR04474; Lee 1998)", fvaRange: fvaPgd.ranges[0].map((v) => +v.toFixed(4)), targetTracer: "glc_12", rank12, rank1, nTracers: ranked.length, top5: ranked.slice(0, 5) });
}

// Case 2: Reductive carboxylation vs Oxidative TCA (Mode B FVA on cytosolic IDH1 / MAR00710; Metallo 2012, Mullen 2012)
console.log("\n== Case 2: Reductive Carboxylation vs Oxidative TCA (Mode B FVA on IDH1; Metallo 2012, Mullen 2012) ==");
{
  const idh1Idx = model.rxns.findIndex((r) => r.id === "MAR00710");
  const cfg = { objective: objIdx, maximize: true, parsimonious: true, measured: meas };
  const f0 = (await sim.runSteady(cfg, sup, st)).frames[0];
  const fvaIdh = await sim.fluxRanges(cfg, sup, [idh1Idx], f0.v, { frac: 0.9, solutions: true, cycles: 0.05 });
  const [vmin, vmax] = fvaIdh.solutions[0];
  console.log(`  IDH1 (${model.rxns[idh1Idx].id}): FVA=[${fvaIdh.ranges[0][0].toFixed(4)}, ${fvaIdh.ranges[0][1].toFixed(4)}]`);
  const keyTracers = candidateTracers.filter((tr) => ["glc_u", "glc_12", "glc_1", "gln_u", "gln_1", "gln_5", "lac_u", "pyr_u"].includes(tr.id));
  const ranked = rankTracers([{ name: "Reductive IDH1", v: vmin, src: f0, sup }, { name: "Oxidative IDH1", v: vmax, src: f0, sup }], keyTracers);
  const rankGlnU = ranked.findIndex((x) => x.id === "gln_u") + 1;
  const rankGln1 = ranked.findIndex((x) => x.id === "gln_1") + 1;
  const rankGln5 = ranked.findIndex((x) => x.id === "gln_5") + 1;
  console.log(`  Top 5 tracers (${keyTracers.length} tested):`, ranked.slice(0, 5).map((x, i) => `#${i+1} ${x.name}(${x.total}) [${x.top3.map((c)=>`${c.met}:${c.sep}x`).join(",")}]`).join(" | "));
  console.log(`  [U-13C5]gln=#${rankGlnU}, [1-13C]gln=#${rankGln1}, [5-13C]gln=#${rankGln5}`);
  results.push({ case: "Reductive carboxylation (Mode B IDH1; Metallo 2012)", fvaRange: fvaIdh.ranges[0].map((v) => +v.toFixed(4)), glnURank: rankGlnU, gln1Rank: rankGln1, gln5Rank: rankGln5, nTracers: ranked.length, top5: ranked.slice(0, 5) });
}

// Case 3: SDHB KO vs WT (Cardaci et al. 2015) — exact [rxnIdx, 0] perturb tuple format!
console.log("\n== Case 3: SDHB knockout vs WT (Cardaci et al. 2015) ==");
{
  const cfgWT = { objective: objIdx, maximize: true, parsimonious: true, measured: meas };
  const fWT = (await sim.runSteady(cfgWT, sup, st)).frames[0];
  const sdhRxns = [];
  model.rxns.forEach((r, j) => { if (r.g && r.g.includes("SDHB") && gprActivity(r.g, new Map([["SDHB", 0]])) === 0) sdhRxns.push(j); });
  const cfgKO = { ...cfgWT, perturb: sdhRxns.map((j) => [j, 0]) };
  const fKO = (await sim.runSteady(cfgKO, sup, st)).frames[0];
  console.log(`  Growth WT=${fWT.objective.toFixed(4)}, SDHB KO=${fKO.objective.toFixed(4)}; SDH (${model.rxns[sdhRxns[0]].id}) WT=${fWT.v[sdhRxns[0]].toFixed(4)}, KO=${fKO.v[sdhRxns[0]].toFixed(4)}`);
  const keyTracers = candidateTracers.filter((tr) => ["glc_u", "glc_12", "glc_1", "gln_u", "gln_1", "gln_5", "lac_u", "pyr_u"].includes(tr.id));
  const ranked = rankTracers([{ name: "WT", v: fWT.v, src: fWT, sup }, { name: "SDHB KO", v: fKO.v, src: fKO, sup }], keyTracers);
  console.log(`  Top 5 tracers (${keyTracers.length} tested):`, ranked.slice(0, 5).map((x, i) => `#${i+1} ${x.name}(${x.total}) [${x.top3.map((c)=>`${c.met}:${c.sep}x`).join(",")}]`).join(" | "));
  results.push({ case: "SDHB KO vs WT (Cardaci 2015)", sdhFluxWt: +fWT.v[sdhRxns[0]].toFixed(4), sdhFluxKo: +fKO.v[sdhRxns[0]].toFixed(4), nTracers: ranked.length, top5: ranked.slice(0, 5) });
}

fs.writeFileSync(R + "tools/validation/design_results.json", JSON.stringify(results, null, 2));
