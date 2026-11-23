// Dynamic runs with scheduled medium events, and isotopologues over time.
// Checks: label only appears after a pulse, a chase removes it again, and the time-aware isotopologue
// sampler agrees with the positional tracer (mean labelled atoms per molecule).
import fs from "node:fs";
import loadHighs from "../vendor/highs/highs.mjs";
import { FbaProblem } from "../js/core/fba.js";
import { presetMedium } from "../js/core/medium.js";
import { Simulator } from "../js/sim.js";
import { IsotopologueSampler } from "../js/core/isotopologues.js";
import { ModelIndex } from "../js/modelindex.js";
const slug = process.argv[2] || "human-gem";
const rd = (f) => JSON.parse(fs.readFileSync(new URL(`../data/${slug}/${f}`, import.meta.url)));
const model = rd("model.json"), structures = rd("structures.json"), maps = rd("atommaps.json"), layout = rd("layout.json");
const highs = await loadHighs({ print: () => {}, printErr: () => {} });
let problem = null;
const backend = async (msg) => {
  if (msg.type === "init") { problem = new FbaProblem(highs, msg.net); return true; }
  const t0 = Date.now(), res = problem.solve(msg.params); res.ms = Date.now() - t0; return res;
};
const sim = new Simulator(model, structures, maps, backend), index = new ModelIndex(model, layout);
const metByName = (name, c) => model.mets.findIndex((m) => m.c === c && model.bases[m.b].n === name);
const cfg = { objective: model.rxns.findIndex((r) => r.obj), maximize: true, parsimonious: true };
const glc = metByName("glucose", "e"), stG = structures[model.bases[model.mets[glc].b].id];
const allC = Float64Array.from(stG.e, (el) => (el === "C" ? 1 : 0));
const medium = (labelled) => presetMedium(model).map((it) => ({ ...it, label: labelled && it.met === glc ? allC : null, fraction: 1 }));
const watch = [["glucose-6-phosphate", "c"], ["pyruvate", "c"], ["citrate", "m"], ["glutamate", "c"]].map(([n, c]) => [n + "[" + c + "]", metByName(n, c)]);
const fmt = (r) => (r ? r.dist.map((p, i) => (p > 0.015 ? `M+${i}:${(100 * p).toFixed(0)}` : null)).filter(Boolean).join(" ") + ` (mean ${r.mean.toFixed(2)} / pos ${r.positional.toFixed(2)})` : "—");
let bad = 0;
async function scenario(title, supplies, events, checks) {
  const t0 = Date.now();
  const dyn = { hours: 24, biomass: 0.5, pool: 0.0005, frames: 30, events };
  const res = await sim.runDynamic(cfg, supplies, structures, dyn);
  const F = res.frames, last = F.length - 1;
  console.log(`\n== ${title}: ${F.length} frames, ${res.solves} FBA, ${((Date.now() - t0) / 1000).toFixed(1)} s, last status ${F[last].status}, biomass ${F[last].biomass.toFixed(3)}`);
  console.log("   glucose mM:", [0, Math.floor(last / 3), Math.floor((2 * last) / 3), last].map((i) => `${F[i].t.toFixed(2)}h=${F[i].conc.get(glc).toFixed(2)}`).join("  "), "| frames with label:", F.filter((f) => f.e).length);
  const S = IsotopologueSampler.dynamic(sim.tracer, index, F, F[last].elements, dyn.pool, res.instant);
  const at = (h) => { let b = 0; F.forEach((f, i) => { if (Math.abs(f.t - h) < Math.abs(F[b].t - h)) b = i; }); return b; };
  let worst = 0;
  for (const h of checks.times) {
    const i = at(h);
    console.log(`   t = ${F[i].t.toFixed(3)} h`);
    for (const [name, m] of watch) {
      const r = F[i].e ? S.mid(m, 3000, 1500, i) : null;
      console.log(`      ${name.padEnd(24)} ${fmt(r)}`);
      // Right after the start or an event the positional tracer takes one coarse implicit step (10 s, longer than
      // the turnover of the fastest pools) while the sampler uses exact waiting times: only later frames are compared.
      const settled = [0, ...events.map((ev) => ev.t)].every((te) => F[i].t < te || F[i].t - te > 0.3);
      if (r && settled) worst = Math.max(worst, Math.abs(r.mean - r.positional) / Math.max(0.15, 0.06 * r.atoms));
    }
  }
  console.log(`   sampler vs positional tracer: worst deviation ${worst.toFixed(2)} (should be < 1)`);
  if (worst >= 1) bad++;
  const g6p = watch[0][1], frac = (i) => (F[i].metE ? F[i].metE[g6p] : 0);
  for (const [h, lo, hi, what] of checks.g6p) {
    const x = frac(at(h)), ok = x >= lo && x <= hi;
    console.log(`   G6P enrichment at ${h} h = ${x.toFixed(3)} (${what}) ${ok ? "ok" : "FAILED"}`);
    if (!ok) bad++;
  }
  return res;
}
await scenario("label from the start", medium(true), [], { times: [0.05, 1, 24], g6p: [[24, 0.95, 1.001, "fully labelled"]] });
await scenario("pulse: 5.5 mM [U-13C]glucose added at 4 h", medium(false), [{ t: 4, type: "add", met: glc, conc: 5.5, vmax: 0.6, km: 1.5, label: allC, fraction: 1 }],
  { times: [4.05, 6, 24], g6p: [[3.9, 0, 0.0001, "no label before the pulse"], [24, 0.4, 0.6, "about half: 5.5 labelled on ~5.4 unlabelled"]] });
await scenario("chase: labelled medium replaced by unlabelled at 4 h", medium(true), [{ t: 4, type: "medium", items: presetMedium(model) }],
  { times: [3.9, 4.05, 24], g6p: [[3.9, 0.95, 1.001, "labelled before the chase"], [24, 0, 0.02, "label washed out"]] });
console.log(bad ? `\nFAILED (${bad})` : "\nOK");
process.exit(bad ? 1 : 0);
