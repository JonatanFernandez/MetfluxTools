// Cell-type contexts: decode data/<slug>/contexts.json, block reactions, and check with the real Simulator that
//  - every context with a measured set reproduces it (not relaxed) with growth > 0 where the set grows,
//  - "penalty" mode keeps the generic objective value,
//  - a sample of contexts still grows on the standard medium without measured fluxes.
//   node tools/test_contexts.mjs [human-gem|mouse-gem] [sample size]
const R = new URL("../", import.meta.url).pathname;
const fs = await import("node:fs");
const { default: loadHighs } = await import(R + "vendor/highs/highs.mjs");
const { FbaProblem } = await import(R + "js/core/fba.js");
const { presetMedium } = await import(R + "js/core/medium.js");
const { measuredPreset, measuredNgam } = await import(R + "js/core/measured.js");
const { Simulator } = await import(R + "js/sim.js");
const { ModelIndex } = await import(R + "js/modelindex.js");
const { offGenes, blockedReactions } = await import(R + "js/core/contexts.js");
globalThis.atob ??= (s) => Buffer.from(s, "base64").toString("binary");

const slug = process.argv[2] || "human-gem", nSample = +(process.argv[3] || 6);
const rd = (f) => JSON.parse(fs.readFileSync(R + `data/${slug}/` + f));
const model = rd("model.json"), structures = rd("structures.json"), atommaps = rd("atommaps.json"), layout = rd("layout.json");
const data = rd("contexts.json");
const highs = await loadHighs({ print: () => {}, printErr: () => {} });
let problem = null;
const backend = async (msg) => { if (msg.type === "init") { problem = new FbaProblem(highs, msg.net); return true; } return problem.solve(msg.params); };
const sim = new Simulator(model, structures, atommaps, backend), index = new ModelIndex(model, layout);
const bio = index.biomassRxn;
let fails = 0;
const check = (ok, msg) => { console.log((ok ? "ok   " : "FAIL ") + msg); if (!ok) fails++; };
const run = async (ctx, mode, set, medium = "default", objective = bio) => {
  const meas = { items: measuredPreset(model, sim.exchanges, set).map((x) => ({ met: x.met, rate: x.rate })), tol: 0.2, ngam: measuredNgam(set) };
  const context = ctx ? { id: ctx.id, mode, block: blockedReactions(model, offGenes(ctx, model.genes)) } : undefined;
  const res = await sim.runSteady({ objective, maximize: true, parsimonious: true, measured: meas, context }, presetMedium(model, medium), structures);
  return res.frames[0];
};

check(data.contexts.length > 0 && data.meta?.license, `${data.contexts.length} contexts, licence: ${data.meta.license.slice(0, 40)}…`);
check(data.contexts.every((c) => !c.failed.length), `every context passes its tasks (${data.contexts.filter((c) => c.failed.length).length} do not)`);
const atp = model.rxns.findIndex((r) => r.id === "MAR03964");
for (const c of data.contexts.filter((c) => c.pairs.length)) {
  for (const set of c.pairs) {
    const medium = set === "liver" ? "williamsE" : "default", obj = set === "fibroQ" ? atp : bio;
    const g = await run(null, "block", set, medium, obj), f = await run(c, "block", set, medium, obj);
    check(f.status === "optimal" && !f.relaxed, `${c.name} + ${set}: ${f.status}${f.relaxed ? " (relaxed)" : ""}, objective ${(+f.objective).toPrecision(3)} (generic ${(+g.objective).toPrecision(3)})`);
  }
}
const g0 = await run(null, "block", "none");
const sample = data.contexts.filter((c, i) => i % Math.max(1, Math.floor(data.contexts.length / nSample)) === 0).slice(0, nSample);
for (const c of sample) {
  const f = await run(c, "block", "none"), p = await run(c, "penalty", "none");
  check(f.status === "optimal" && f.objective > 0, `${c.name}: growth ${(+f.objective).toPrecision(3)} h-1 blocked (generic ${(+g0.objective).toPrecision(3)})`);
  check(Math.abs(p.objective - g0.objective) < 1e-6 * Math.max(1, g0.objective), `${c.name}: penalty mode keeps the objective (${(+p.objective).toPrecision(4)})`);
}
console.log(fails ? `${fails} failed` : "all passed");
process.exit(fails ? 1 : 0);
