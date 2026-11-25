// Validation Suite 7: Cross-model Performance & Stoichiometric Mass-Balance Verification Benchmark
// Measures across ALL 8 bundled GEMs (Human-GEM, Mouse-GEM, Rat-GEM, Zebrafish-GEM, Fruitfly-GEM,
// Worm-GEM, yeast-GEM, E. coli iML1515):
//   - JSON model size & load time
//   - FBA solve time & stoichiometric steady-state residual ||S v||_inf (must be < 1e-7)
//   - pFBA solve time
//   - Pathway FVA time (central glycolysis/TCA subsystem)
//   - Full single-gene essentiality screen time
//   - Positional 13C tracing time + 10-metabolite isotopologue sampling time
//   - ULM regulatory inference time (for models with regulons.json)
//
// Run: node tools/validation/benchmark_perf.mjs
import fs from "node:fs";
import { performance } from "node:perf_hooks";
const R = new URL("../../", import.meta.url).pathname;
const { default: loadHighs } = await import(R + "vendor/highs/highs.mjs");
const { FbaProblem } = await import(R + "js/core/fba.js");
const { presetMedium, defaultMediumFor } = await import(R + "js/core/medium.js");
const { Simulator } = await import(R + "js/sim.js");
const { ModelIndex } = await import(R + "js/modelindex.js");
const { IsotopologueSampler } = await import(R + "js/core/isotopologues.js");
const { gprActivity } = await import(R + "js/core/gpr.js");
const { ulm } = await import(R + "js/core/reginfer.js");

const MODELS = [
  { slug: "human-gem",     medium: "defaultSerum" },
  { slug: "mouse-gem",     medium: "defaultSerum" },
  { slug: "rat-gem",       medium: "defaultSerum" },
  { slug: "zebrafish-gem", medium: "default" },
  { slug: "fruitfly-gem",  medium: "schneider" },
  { slug: "worm-gem",      medium: "default" },
  { slug: "yeast-gem",     medium: "verduyn" },
  { slug: "ecoli-gem",     medium: "m9glc" },
];

const highs = await loadHighs({ print: () => {}, printErr: () => {} });
const rows = [];

function maxStoichResidual(model, v) {
  const res = new Float64Array(model.mets.length);
  for (let j = 0; j < model.rxns.length; j++) {
    const vj = v[j];
    if (!vj) continue;
    const r = model.rxns[j];
    for (let k = 0; k < r.m.length; k++) res[r.m[k]] += r.c[k] * vj;
  }
  let mx = 0;
  for (let i = 0; i < res.length; i++) {
    const a = Math.abs(res[i]);
    if (a > mx) mx = a;
  }
  return mx;
}

console.log("Model          |  Mets |  Rxns | Genes | Load ms | FBA ms | pFBA ms | FVA(pw) ms | Screen ms (LPs) | Trace ms | MID10 ms | ||Sv||_inf | Obj");
console.log("-".repeat(138));

for (const m of MODELS) {
  const dir = R + `data/${m.slug}/`;
  const tLoad0 = performance.now();
  const model = JSON.parse(fs.readFileSync(dir + "model.json", "utf8"));
  const st = JSON.parse(fs.readFileSync(dir + "structures.json", "utf8"));
  const am = JSON.parse(fs.readFileSync(dir + "atommaps.json", "utf8"));
  const layout = JSON.parse(fs.readFileSync(dir + "layout.json", "utf8"));
  const index = new ModelIndex(model, layout);
  const loadMs = performance.now() - tLoad0;

  let pb;
  const backend = async (msg) => {
    if (msg.type === "init") { pb = new FbaProblem(highs, msg.net); return true; }
    const r = pb.solve(msg.params); r.ms = 0; return r;
  };
  const sim = new Simulator(model, st, am, backend);
  const sup = presetMedium(model, m.medium);
  // Use MAR13082 for mouse/rat as the app does for generic biomass
  let objIdx = model.rxns.findIndex((r) => r.obj);
  if ((m.slug === "mouse-gem" || m.slug === "rat-gem")) {
    const alt = model.rxns.findIndex((r) => r.id === "MAR13082");
    if (alt >= 0) objIdx = alt;
  }

  // 1. Plain FBA
  const tFba0 = performance.now();
  const fFba = (await sim.runSteady({ objective: objIdx, maximize: true, parsimonious: false }, sup, st)).frames[0];
  const fbaMs = performance.now() - tFba0;
  const svRes = fFba.status === "optimal" ? maxStoichResidual(model, fFba.v) : NaN;

  // 2. pFBA
  const tPfba0 = performance.now();
  const fPfba = (await sim.runSteady({ objective: objIdx, maximize: true, parsimonious: true }, sup, st)).frames[0];
  const pfbaMs = performance.now() - tPfba0;

  // 3. Pathway FVA (pick Glycolysis / Gluconeogenesis or first active pathway with 10-40 rxns)
  const pwIdx = model.subsystems.findIndex((s) => /glycolysis/i.test(s));
  const pwRxns = [];
  model.rxns.forEach((r, j) => { if (r.sub === pwIdx && Math.abs(fPfba.v[j]) > 1e-9) pwRxns.push(j); });
  const fvaList = pwRxns.slice(0, 20);
  const tFva0 = performance.now();
  if (fvaList.length) await sim.fluxRanges({ objective: objIdx, maximize: true, parsimonious: true }, sup, fvaList, fPfba.v, { cycles: 0.05 });
  const fvaMs = performance.now() - tFva0;

  // 4. Full gene essentiality screen (same grouped LP strategy as app.js)
  const geneRxns = new Map();
  model.rxns.forEach((r, j) => {
    if (!r.g) return;
    for (const tok of new Set(r.g.match(/[^\s()]+/g).filter((x) => !/^(and|or)$/i.test(x)))) {
      if (!geneRxns.has(tok)) geneRxns.set(tok, []);
      geneRxns.get(tok).push(j);
    }
  });
  const tScr0 = performance.now();
  const groups = new Map();
  for (const [sym, rx] of geneRxns) {
    if (!rx.some((j) => Math.abs(fPfba.v[j]) > 1e-9)) continue;
    const levels = new Map([[sym, 0]]);
    const lost = rx.filter((j) => gprActivity(model.rxns[j].g, levels) === 0).sort((a, b) => a - b);
    if (!lost.length) continue;
    const key = lost.join(",");
    if (!groups.has(key)) groups.set(key, Int32Array.from(lost));
  }
  const grpList = [...groups.values()];
  if (grpList.length) await sim.screenKnockouts({ objective: objIdx, maximize: true }, sup, grpList);
  const scrMs = performance.now() - tScr0;

  // 5. Positional 13C tracing ([U-13C6]glucose) + 10-metabolite MID sampling
  const glcSup = sup.find((s) => /^glucose$/i.test(model.bases[model.mets[s.met].b].n));
  let traceMs = NaN, midMs = NaN;
  if (glcSup) {
    const base = model.mets[glcSup.met].b;
    const sObj = st[model.bases[base].id];
    if (sObj?.atoms) {
      const lab = sObj.atoms.map((el, i) => (el === "C" ? { atom: i, iso: "13C", frac: 1 } : null)).filter(Boolean);
      const tTr0 = performance.now();
      const tr2 = sim.traceFlux(fPfba.v, fPfba, new Map([[glcSup.met, lab]]));
      traceMs = performance.now() - tTr0;
      const sampler = new IsotopologueSampler(sim.tracer, index, fPfba.v, tr2.inputs, tr2.e, new Set(["C"]));
      // Pick 10 active intracellular labelled metabolites
      const activeMets = [];
      for (let i = 0; i < model.mets.length && activeMets.length < 10; i++) {
        if (tr2.e[i]?.some((x) => x > 0.05) && model.mets[i].c !== model.compartments.findIndex((c) => /extracellular|external|medium/i.test(c))) {
          activeMets.push(i);
        }
      }
      const tMid0 = performance.now();
      for (const mi of activeMets) sampler.mid(mi, 1000, 100);
      midMs = performance.now() - tMid0;
    }
  }

  // 6. ULM regulatory inference timing (if regulons.json exists)
  let ulmMs = null, nTfs = 0;
  if (fs.existsSync(dir + "regulons.json")) {
    const reg = JSON.parse(fs.readFileSync(dir + "regulons.json", "utf8"));
    const vals = new Map();
    for (const g of reg.genes) vals.set(g, (Math.sin(vals.size * 12.9898) * 43758.5453) % 2 - 1);
    const tU0 = performance.now();
    const u = ulm(reg, vals, { minTargets: 5 });
    ulmMs = +(performance.now() - tU0).toFixed(1);
    nTfs = u.length;
  }

  console.log(
    `${m.slug.padEnd(14)} | ${String(model.mets.length).padStart(5)} | ${String(model.rxns.length).padStart(5)} | ${String(model.genes.length).padStart(5)} | ${loadMs.toFixed(0).padStart(7)} | ${fbaMs.toFixed(0).padStart(6)} | ${pfbaMs.toFixed(0).padStart(7)} | ${fvaMs.toFixed(0).padStart(10)} | ${(scrMs.toFixed(0) + " (" + grpList.length + ")").padStart(15)} | ${(Number.isFinite(traceMs) ? traceMs.toFixed(0) : "NA").padStart(8)} | ${(Number.isFinite(midMs) ? midMs.toFixed(0) : "NA").padStart(8)} | ${svRes.toExponential(1).padStart(10)} | ${fFba.objective.toFixed(4)}`
  );
  rows.push({
    slug: m.slug, medium: m.medium,
    mets: model.mets.length, rxns: model.rxns.length, genes: model.genes.length,
    loadMs: +loadMs.toFixed(1), fbaMs: +fbaMs.toFixed(1), pfbaMs: +pfbaMs.toFixed(1),
    fvaMs: +fvaMs.toFixed(1), fvaRxns: fvaList.length,
    screenMs: +scrMs.toFixed(1), screenLps: grpList.length,
    traceMs: Number.isFinite(traceMs) ? +traceMs.toFixed(1) : null,
    mid10Ms: Number.isFinite(midMs) ? +midMs.toFixed(1) : null,
    ulmMs, nTfs,
    maxStoichResidual: svRes,
    objective: +fFba.objective.toFixed(5),
    status: fFba.status,
  });
}

const maxResAll = Math.max(...rows.map((r) => r.maxStoichResidual));
console.log(`\nMax stoichiometric residual ||S v||_inf across all 8 GEMs: ${maxResAll.toExponential(2)} (${maxResAll < 1e-7 ? "PASS < 1e-7" : "CHECK"})`);
fs.writeFileSync(R + "tools/validation/benchmark_perf_results.json", JSON.stringify(rows, null, 2));
