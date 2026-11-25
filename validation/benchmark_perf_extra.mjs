// Fast addendum to benchmark_perf.mjs: measures pathway FVA, 13C tracing, 10-metabolite MID sampling,
// and ULM timings across all 8 GEMs (using exact r.s subsystem index & sObj.a atom list) and merges
// into tools/validation/benchmark_perf_results.json without re-running the 8 full KO screens.
import fs from "node:fs";
import { performance } from "node:perf_hooks";
const R = new URL("../../", import.meta.url).pathname;
const { default: loadHighs } = await import(R + "vendor/highs/highs.mjs");
const { FbaProblem } = await import(R + "js/core/fba.js");
const { presetMedium } = await import(R + "js/core/medium.js");
const { Simulator } = await import(R + "js/sim.js");
const { ModelIndex } = await import(R + "js/modelindex.js");
const { IsotopologueSampler } = await import(R + "js/core/isotopologues.js");
const { ulm } = await import(R + "js/core/reginfer.js");

const rows = JSON.parse(fs.readFileSync(R + "tools/validation/benchmark_perf_results.json", "utf8"));
const highs = await loadHighs({ print: () => {}, printErr: () => {} });

for (const row of rows) {
  const dir = R + `data/${row.slug}/`;
  const model = JSON.parse(fs.readFileSync(dir + "model.json", "utf8"));
  const st = JSON.parse(fs.readFileSync(dir + "structures.json", "utf8"));
  const am = JSON.parse(fs.readFileSync(dir + "atommaps.json", "utf8"));
  const layout = JSON.parse(fs.readFileSync(dir + "layout.json", "utf8"));
  const index = new ModelIndex(model, layout);

  let pb;
  const backend = async (msg) => {
    if (msg.type === "init") { pb = new FbaProblem(highs, msg.net); return true; }
    const r = pb.solve(msg.params); r.ms = 0; return r;
  };
  const sim = new Simulator(model, st, am, backend);
  const sup = presetMedium(model, row.medium);
  let objIdx = model.rxns.findIndex((r) => r.obj);
  if ((row.slug === "mouse-gem" || row.slug === "rat-gem")) {
    const alt = model.rxns.findIndex((r) => r.id === "MAR13082");
    if (alt >= 0) objIdx = alt;
  }
  const cfg = { objective: objIdx, maximize: true, parsimonious: true };
  const fPfba = (await sim.runSteady(cfg, sup, st)).frames[0];

  // Pathway FVA: Glycolysis (r.ss === pwIdx) up to 15 active reactions
  const pwIdx = model.subsystems.findIndex((s) => /glycolysis/i.test(s));
  const pwRxns = [];
  model.rxns.forEach((r, j) => { if (r.ss === pwIdx && Math.abs(fPfba.v[j]) > 1e-9) pwRxns.push(j); });
  const fvaList = pwRxns.slice(0, 15);
  const tFva0 = performance.now();
  if (fvaList.length) await sim.fluxRanges(cfg, sup, fvaList, fPfba.v, { cycles: 0.05 });
  row.fvaMs = +(performance.now() - tFva0).toFixed(1);
  row.fvaRxns = fvaList.length;

  // Positional 13C tracing ([U-13C6]glucose) + 10-metabolite MID sampling
  const glcSup = sup.find((s) => /^glucose$/i.test(model.bases[model.mets[s.met].b].n));
  if (glcSup) {
    const base = model.mets[glcSup.met].b;
    const sObj = st[model.bases[base].id];
    const atoms = sObj?.e;
    if (atoms) {
      const lab = new Float32Array(atoms.length);
      for (let i = 0; i < atoms.length; i++) if (atoms[i] === "C") lab[i] = 1;
      const tTr0 = performance.now();
      const tr2 = sim.traceFlux(fPfba.v, fPfba, new Map([[glcSup.met, lab]]));
      row.traceMs = +(performance.now() - tTr0).toFixed(1);
      const sampler = new IsotopologueSampler(sim.tracer, index, fPfba.v, tr2.inputs, tr2.e, new Set(["C"]));
      const activeMets = [];
      for (let i = 0; i < model.mets.length && activeMets.length < 10; i++) {
        if (model.mets[i].c === "e" || model.mets[i].c === "s") continue;
        let prod = 0;
        for (let p = index.rxnStart[i]; p < index.rxnStart[i + 1]; p++) {
          const x = fPfba.v[index.rxnList[p]] * index.rxnCoef[p];
          if (x > 0) prod += x;
        }
        if (prod > 1e-4 && st[model.bases[model.mets[i].b].id]?.e?.includes("C")) activeMets.push(i);
      }
      const tMid0 = performance.now();
      for (const mi of activeMets) sampler.mid(mi, 1000, 100);
      row.mid10Ms = +(performance.now() - tMid0).toFixed(1);
    }
  }
  console.log(`${row.slug.padEnd(14)} | FVA(${String(row.fvaRxns).padStart(2)} rxns)=${String(row.fvaMs).padStart(7)} ms | Trace=${String(row.traceMs).padStart(6)} ms | MID10=${String(row.mid10Ms).padStart(6)} ms | ULM=${row.ulmMs != null ? row.ulmMs + " ms (" + row.nTfs + " TFs)" : "—"}`);
}
fs.writeFileSync(R + "tools/validation/benchmark_perf_results.json", JSON.stringify(rows, null, 2));
