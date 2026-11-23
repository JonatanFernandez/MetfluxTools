// Numbers for the docs examples on comparison and model editing (not a test).
import fs from "node:fs";
import { createRequire } from "node:module";
import loadHighs from "../vendor/highs/highs.mjs";
import { FbaProblem } from "../js/core/fba.js";
import { presetMedium } from "../js/core/medium.js";
import { measuredPreset, measuredNgam } from "../js/core/measured.js";
import { Simulator } from "../js/sim.js";
import { ModelIndex } from "../js/modelindex.js";
import { IsotopologueSampler } from "../js/core/isotopologues.js";
import { applyEdits, modelLookups } from "../js/core/edits.js";
import { setRDKit, smilesToStructure, mapFromMappedSmiles } from "../js/core/chem.js";
const rd = (f) => JSON.parse(fs.readFileSync(new URL(`../data/human-gem/${f}`, import.meta.url)));
const model = rd("model.json"), structures = rd("structures.json"), atommaps = rd("atommaps.json"), layout = rd("layout.json");
setRDKit(await createRequire(import.meta.url)("../vendor/rdkit/RDKit_minimal.js")());
const highs = await loadHighs({ print: () => {}, printErr: () => {} });
let problem = null;
const backend = async (msg) => { if (msg.type === "init") { problem = new FbaProblem(highs, msg.net); return true; } const r = problem.solve(msg.params); r.ms = 0; return r; };
const glcB = model.bases.findIndex((b) => b.n === "glucose");
async function run(mdl, st, am, { o2 = null, objective = null, measured = true } = {}) {
  const sim = new Simulator(mdl, st, am, backend), index = new ModelIndex(mdl, layout);
  const glc = mdl.mets.findIndex((m) => m.c === "e" && m.b === glcB);
  const allC = Float64Array.from(st[mdl.bases[glcB].id].e, (el) => (el === "C" ? 1 : 0));
  const sup = presetMedium(mdl).map((it) => ({ ...it, label: it.met === glc ? allC : null, fraction: 1, conc: o2 !== null && mdl.bases[mdl.mets[it.met].b].n === "O2" ? o2 : it.conc }));
  const meas = measured ? { items: measuredPreset(mdl, sim.exchanges, "nci60").map((x) => ({ met: x.met, rate: x.rate })), tol: 0.2, ngam: measuredNgam("nci60") } : undefined;
  const cfg = { objective: objective ?? mdl.rxns.findIndex((r) => r.obj), maximize: true, parsimonious: true, measured: meas };
  const res = await sim.runSteady(cfg, sup, st);
  const f = res.frames[0];
  const smp = f.e ? new IsotopologueSampler(sim.tracer, index, f.v, f.inputs, f.e, f.elements) : null;
  const ex = (name) => { for (const [m, e] of sim.exchanges) if (mdl.bases[mdl.mets[m].b].n === name) return f.v[e.rxn] * e.dir; return 0; };
  const mid = (name, c) => { const m = mdl.mets.findIndex((x) => x.c === c && mdl.bases[x.b].n === name); const r = smp?.mid(m, 4000, 400); return r ? r.dist.map((p, i) => (p > 0.02 ? `M+${i} ${(100 * p).toFixed(0)}%` : null)).filter(Boolean).join(", ") : "-"; };
  return { f, sim, ex, mid, relaxed: f.relaxed };
}
const show = (tag, R) => console.log(`${tag}: growth ${R.f.objective.toFixed(4)}${R.relaxed ? " (relaxed)" : ""}  glc ${R.ex("glucose").toFixed(3)}  lac ${R.ex("L-lactate").toFixed(3)}  O2 ${R.ex("O2").toFixed(3)} | lactate[c] ${R.mid("L-lactate", "c")} | citrate[m] ${R.mid("citrate", "m")} | fumarate[m] ${R.mid("fumarate", "m")} | malate[c] ${R.mid("malate", "c")}`);
const A = await run(model, structures, atommaps); show("A control", A);
const B = await run(model, structures, atommaps, { o2: 0.0003 }); show("B hypoxia", B);
// biggest exchange changes
const rows = [];
for (const [m, e] of A.sim.exchanges) { const a = A.f.v[e.rxn] * e.dir, b = B.f.v[e.rxn] * e.dir; if (Math.abs(a - b) > 1e-4) rows.push([model.bases[model.mets[m].b].n, a.toFixed(3), b.toFixed(3)]); }
rows.sort((x, y) => Math.abs(y[2] - y[1]) - Math.abs(x[2] - x[1])); console.log("exchanges A -> B:", rows.slice(0, 8).map((r) => r.join(" ")).join(" | "));
// editing: acetoin pathway
const L = modelLookups(model), id = (n, c) => model.mets.find((m) => m.c === c && model.bases[m.b].n === n).id;
const ac = smilesToStructure("CC(O)C(C)=O");
let ops = [{ op: "addMet", comp: "c", name: "acetoin", formula: ac.formula, charge: 0, smiles: ac.smiles, structure: ac.st }];
let E = applyEdits({ model, structures, atommaps }, ops);
const acC = E.model.mets[E.model.mets.length - 1].id, acE = acC.slice(0, -1) + "e";
ops.push({ op: "addMet", base: acC.slice(0, -1), comp: "e" });
const mets = [[id("pyruvate", "c"), -2], [id("H+", "c"), -1], [acC, 1], [id("CO2", "c"), 2]];
E = applyEdits({ model, structures, atommaps }, ops);
const parts = mets.map(([mid, c], loc) => { const m = E.model.mets.find((x) => x.id === mid); return { loc, side: Math.sign(c), coef: c, st: E.structures[E.model.bases[m.b].id] }; });
const map = mapFromMappedSmiles("[CH3:1][C:2](=[O:3])[C:8](=[O:7])[O-:9].[CH3:5][C:4](=[O:6])[C:11](=[O:10])[O-:12].[H+]>>[CH3:1][C:2](=[O:3])[CH:4]([CH3:5])[OH:6].[O:7]=[C:8]=[O:9].[O:10]=[C:11]=[O:12]", parts);
ops.push({ op: "addRxn", id: "USR_R001", name: "acetolactate synthase + decarboxylase", subsystem: "User pathway", mets, lb: 0, ub: 1000, gpr: "alsS and alsD", atommap: { k: "user", s: map.s } });
ops.push({ op: "addRxn", id: "USR_R002", name: "acetoin transport", subsystem: "User pathway", mets: [[acC, -1], [acE, 1]], lb: -1000, ub: 1000, atommap: { k: "id", s: [[0, 1, -1]] } });
ops.push({ op: "addRxn", id: "USR_R003", name: "acetoin exchange", subsystem: "Exchange/demand reactions", mets: [[acE, -1]], lb: 0, ub: 1000 });
E = applyEdits({ model, structures, atommaps }, ops);
console.log("edit results:", E.results.map((r) => r.ok).join(","));
const C = await run(E.model, E.structures, E.atommaps); show("C acetoin pathway, growth objective", C); console.log("   acetoin secretion", C.ex("acetoin").toFixed(4));
const D = await run(E.model, E.structures, E.atommaps, { objective: E.model.rxns.findIndex((r) => r.id === "USR_R003") }); console.log(`D max acetoin: ${D.f.objective.toFixed(3)}${D.relaxed ? " (relaxed)" : ""}  glc ${D.ex("glucose").toFixed(3)} lac ${D.ex("L-lactate").toFixed(3)} | acetoin[c] ${D.mid("acetoin", "c")} | pyruvate[c] ${D.mid("pyruvate", "c")}`);
// remove the four L-lactate dehydrogenases
const ldh = ["MAR04280", "MAR04281", "MAR04388", "MAR06018"].map((x) => ({ op: "removeRxn", id: x }));
const F = applyEdits({ model, structures, atommaps }, ldh);
const G = await run(F.model, F.structures, F.atommaps); show("G no LDH", G);
for (const [m, e] of G.sim.exchanges) { const a = A.f.v[e.rxn] * e.dir, b = G.f.v[e.rxn] * e.dir; if (Math.abs(a - b) > 0.02) console.log("   ", model.bases[model.mets[m].b].n, a.toFixed(3), "->", b.toFixed(3)); }
