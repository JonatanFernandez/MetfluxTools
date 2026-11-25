// Reproduces the charge-balance defect of Human-GEM 2.1.1 described in
// tools/upstream/human-gem-charge-balance/README.md: a cycle of half-reactions whose net equation is
// NAD(P)H -> NAD(P)+ + H+, which lets the model oxidise NAD(P)H with no electron acceptor.
//   node tools/validation/humangem_charge_cycle.mjs
import fs from "node:fs";
import loadHighs from "../../vendor/highs/highs.mjs";
import { packNetwork } from "../../js/core/fba.js";
import { presetMedium, findExchanges, applyMedium } from "../../js/core/medium.js";
import { measuredPreset } from "../../js/core/measured.js";
import { ExchangeFit, fitBounds, exchangeBalance } from "../../js/core/fitexchange.js";

const R = new URL("../../", import.meta.url).pathname;
const model = JSON.parse(fs.readFileSync(R + "data/human-gem/model.json"));
const thermo = JSON.parse(fs.readFileSync(R + "data/human-gem/thermo.json"));
const highs = await loadHighs({ print: () => {}, printErr: () => {} });
const fitter = new ExchangeFit(highs, packNetwork(model));
const ex = findExchanges(model);
const nm = (m) => model.bases[model.mets[m].b].n, ff = (m) => model.bases[model.mets[m].b].f;
const growthRxn = model.rxns.findIndex((r) => r.id === "MFL_biomass_nocof");
const medium = presetMedium(model, "defaultSerum");
const items = measuredPreset(model, ex, "hek293");
const base = applyMedium(model, ex, medium, model.rxns.map((r) => r.lb), model.rxns.map((r) => r.ub));
const { lb, ub } = fitBounds(model, ex, base.lb, base.ub,
  { measuredMets: new Set(items.map((x) => x.met)), presentMets: new Set(medium.filter((q) => q.unlimited || q.conc > 0).map((q) => q.met)), openRxns: [growthRxn] });
// the glucose uptake set 5x too high: the surplus carbon has to leave somehow
const measurements = items.map((x) => ({ rxn: ex.get(x.met).rxn, dir: ex.get(x.met).dir, met: x.met, name: nm(x.met),
  rate: x.rate * (nm(x.met) === "glucose" ? 5 : 1), sigma: Math.max(1e-4, Math.abs(x.rate) * 0.1) }));
measurements.push({ rxn: growthRxn, dir: 1, met: -1, name: "growth", rate: 0.028, sigma: 0.0028 });
const res = fitter.fit({ lb, ub, measurements });
const bal = exchangeBalance(model, ex, res.v);
console.log(`fit: deviation ${res.deviation.toFixed(4)} (every measured rate reproduced), growth ${res.v[growthRxn].toFixed(4)} 1/h`);
console.log(`balance: RQ = ${bal.rq.toFixed(1)}, O2 ${bal.o2.toFixed(3)}, CO2 ${bal.co2.toFixed(2)}, H+ ${bal.hplus.toFixed(2)}, electron gap ${bal.electrons.gap.toFixed(1)}`);
for (const w of bal.warnings) console.log("  warning: " + w);

const IDS = ["MAR03167", "MAR06606", "MAR03265", "MAR06396", "MAR03868", "MAR08620"];
const atoms = (f) => { const o = {}; if (!f) return null; for (const m of String(f).matchAll(/([A-Z][a-z]?)(\d*)/g)) if (m[1]) o[m[1]] = (o[m[1]] || 0) + (m[2] ? +m[2] : 1); return o; };
const net = {};
console.log("\nthe reactions, with their flux in this solution:");
for (const id of IDS) {
  const j = model.rxns.findIndex((q) => q.id === id), rx = model.rxns[j], v = res.v[j];
  const eq = rx.m.map((mi, q) => `${rx.c[q] > 0 ? "+" : ""}${rx.c[q]} ${nm(mi)}`).join(" ");
  const b = {}; let unk = false;
  rx.m.forEach((mi, q) => { const a = atoms(ff(mi)); if (!a) { unk = true; return; } for (const [e, n] of Object.entries(a)) b[e] = (b[e] || 0) + rx.c[q] * n; });
  const off = Object.entries(b).filter(([, x]) => Math.abs(x) > 1e-9);
  console.log(`  ${id}  v = ${v.toFixed(3).padStart(8)}  elements ${unk ? "R group" : off.length ? off.map(([e, x]) => e + ":" + x).join(",") : "balanced"}  dG ${thermo.rxns[id] ? thermo.rxns[id][0] + " ± " + thermo.rxns[id][1] : "ABSENT"}`);
  console.log(`      ${eq}`);
  if (Math.abs(v) > 1e-9) rx.m.forEach((mi, q) => { const k = nm(mi); net[k] = (net[k] || 0) + v * rx.c[q]; });
}
console.log("\nnet of the six reactions at their fluxes:");
console.log("  " + Object.entries(net).filter(([, x]) => Math.abs(x) > 1e-6).map(([k, x]) => `${x > 0 ? "+" : ""}${x.toFixed(2)} ${k}`).join("   "));
console.log("  -> NAD(P)H is oxidised and a proton appears; the two electrons have no acceptor.");
