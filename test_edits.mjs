// Model edits: validation rules, applying an edit list, FBA and isotope tracing on the edited model.
//   node tools/test_edits.mjs [slug]
import fs from "node:fs";
import { createRequire } from "node:module";
import loadHighs from "../vendor/highs/highs.mjs";
import { FbaProblem } from "../js/core/fba.js";
import { presetMedium } from "../js/core/medium.js";
import { Simulator } from "../js/sim.js";
import { ModelIndex } from "../js/modelindex.js";
import { applyEdits, checkReaction, gprError, modelIssues, modelLookups, editedLayout, hasErrors } from "../js/core/edits.js";
import { setRDKit, smilesToStructure, mcsMap, mapFromMappedSmiles } from "../js/core/chem.js";

const slug = process.argv[2] || "human-gem";
const rd = (f) => JSON.parse(fs.readFileSync(new URL(`../data/${slug}/${f}`, import.meta.url)));
const model = rd("model.json"), structures = rd("structures.json"), atommaps = rd("atommaps.json"), layout = rd("layout.json");
const require = createRequire(import.meta.url);
setRDKit(await require("../vendor/rdkit/RDKit_minimal.js")());
const highs = await loadHighs({ print: () => {}, printErr: () => {} });
let bad = 0;
const check = (ok, what) => { console.log(`${ok ? "ok  " : "FAIL"} ${what}`); if (!ok) bad++; };
const keys = (iss) => iss.map((x) => x.level[0] + ":" + x.key);

// ---- GPR syntax
check(gprError("A and (B or C)") === null, "GPR: valid rule");
check(gprError("A and (B or C") === "gpr.parens", "GPR: unbalanced parenthesis");
check(gprError("A B") === "gpr.missingOp", "GPR: missing operator");
check(gprError("A and") === "gpr.operator", "GPR: dangling operator");

// ---- reaction checks
const L = modelLookups(model);
const metId = (name, c) => model.mets.find((m) => m.c === c && model.bases[m.b].n === name).id;
const pyr = metId("pyruvate", "c"), lac = metId("L-lactate", "c"), nadh = metId("NADH", "c"), nad = metId("NAD+", "c"), h = metId("H+", "c"), co2 = metId("CO2", "c");
let r = checkReaction(model, L, { id: "T1", name: "x", mets: [[pyr, -1], [nadh, -1], [h, -1], [lac, 1], [nad, 1]], lb: 0, ub: 1000 }, { structures });
check(r.issues.some((x) => x.key === "chk.duplicate"), "duplicate of lactate dehydrogenase detected: " + keys(r.issues).join(" "));
r = checkReaction(model, L, { id: "T2", name: "x", mets: [[pyr, -1], [lac, 1]], lb: 0, ub: 1000 }, { structures });
check(r.issues.some((x) => x.key === "chk.hUnbalanced") && !hasErrors(r.issues), "pyruvate -> lactate: only H unbalanced (warning)");
r = checkReaction(model, L, { id: "T3", name: "x", mets: [[pyr, -1], [co2, 1]], lb: 0, ub: 1000 }, { structures });
check(r.issues.some((x) => x.key === "chk.unbalanced" && x.level === "error"), "pyruvate -> CO2: carbon imbalance is an error");
r = checkReaction(model, L, { id: model.rxns[0].id, name: "x", mets: [[pyr, -1], [lac, 1]], lb: 5, ub: 1 }, { structures });
check(r.issues.some((x) => x.key === "chk.rxnDupId") && r.issues.some((x) => x.key === "chk.lbAboveUb"), "duplicate id and lb > ub");

// ---- a new compound from SMILES and a heterologous pathway: 2 pyruvate -> acetoin + 2 CO2 (bacterial ALS + ALDC)
const acetoin = smilesToStructure("CC(O)C(C)=O");
check(acetoin.formula === "C4H8O2" && acetoin.st.e.length === 6, "acetoin read from SMILES: " + acetoin.formula);
const ops = [
  { op: "addMet", comp: "c", name: "acetoin", formula: acetoin.formula, charge: 0, smiles: acetoin.smiles, structure: acetoin.st },
];
let E = applyEdits({ model, structures, atommaps }, ops);
const acC = E.model.mets[E.model.mets.length - 1].id;
ops.push({ op: "addMet", comp: "e", base: acC.slice(0, -1) });
const acE = acC.slice(0, -1) + "e";
const h2 = metId("H+", "c");
const draft = { op: "addRxn", id: "USR_ALS", name: "acetolactate synthase + decarboxylase (lumped)", subsystem: "User pathway", mets: [[pyr, -2], [h2, -1], [acC, 1], [co2, 2]], lb: 0, ub: 1000, gpr: "alsS and alsD" };
E = applyEdits({ model, structures, atommaps }, ops);
const parts = draft.mets.map(([id, c], loc) => { const m = E.model.mets.find((x) => x.id === id); return { loc, side: Math.sign(c), coef: c, st: E.structures[E.model.bases[m.b].id] }; });
const segs = mcsMap(parts);
const mappedC = segs.filter((s) => s[1] === 2).reduce((a, s) => a + s[2].length / 2, 0);
check(mappedC >= 5, `MCS maps ${mappedC} of acetoin's 6 heavy atoms`);
draft.atommap = { k: "mcs", s: segs };
ops.push(draft);
ops.push({ op: "addRxn", id: "USR_ACT", name: "acetoin transport", subsystem: "User pathway", mets: [[acC, -1], [acE, 1]], lb: -1000, ub: 1000, atommap: { k: "id", s: [[0, 1, -1]] } });
ops.push({ op: "addRxn", id: "USR_EXA", name: "acetoin exchange", subsystem: "Exchange/demand reactions", mets: [[acE, -1]], lb: 0, ub: 1000 });
// remove both lactate dehydrogenases (cytosol) and a removal that must be refused
const bOf = (id) => model.mets[L.met.get(id)].b;
const ldh = model.rxns.filter((x) => x.m.some((m) => model.mets[m].b === bOf(pyr)) && x.m.some((m) => model.mets[m].b === bOf(lac))).map((x) => x.id);
for (const id of ldh) ops.push({ op: "removeRxn", id });
ops.push({ op: "removeMet", id: pyr });
E = applyEdits({ model, structures, atommaps }, ops);
console.log("   results:", E.results.map((x, i) => `${ops[i].op}:${x.ok ? "ok" : "ERR"}[${keys(x.issues).join(",")}]`).join("  "));
check(E.results.slice(0, -1).every((x) => x.ok), "every valid operation applied");
check(!E.results[E.results.length - 1].ok && E.results[E.results.length - 1].issues[0].key === "chk.metInUse", "removing a metabolite still in use is refused");
check(E.model.mets.length === model.mets.length + 2 && E.model.rxns.length === model.rxns.length + 3, "indices of the published model unchanged, additions appended");
console.log("   model issues:", keys(modelIssues(E)).join(" ") || "none");
const lay = editedLayout(layout, E.model);
const userSs = E.model.subsystems.indexOf("User pathway");
check(!!lay.subs[userSs] && lay.subs[userSs].r.length === 2, "new pathway gets a map with its 2 reactions");
check(ldh.every((id) => !Object.values(lay.subs).some((Lx) => Lx.r.includes(L.rxn.get(id)))), "removed reactions left out of every map");

// ---- FBA + tracing on the edited model
let problem = null;
const backend = async (msg) => {
  if (msg.type === "init") { problem = new FbaProblem(highs, msg.net); return true; }
  const t0 = Date.now(), res = problem.solve(msg.params); res.ms = Date.now() - t0; return res;
};
const cfg = { objective: model.rxns.findIndex((x) => x.obj), maximize: true, parsimonious: true };
const glc = model.mets.findIndex((m) => m.c === "e" && model.bases[m.b].n === "glucose");
const allC = Float64Array.from(structures[model.bases[model.mets[glc].b].id].e, (el) => (el === "C" ? 1 : 0));
const supplies = (mdl) => presetMedium(mdl).map((it) => ({ ...it, label: it.met === glc ? allC : null, fraction: 1 }));
async function solve(mdl, st, am, extraCfg = {}) {
  const sim = new Simulator(mdl, st, am, backend);
  const res = await sim.runSteady({ ...cfg, ...extraCfg }, supplies(mdl), st);
  return { sim, f: res.frames[0], index: new ModelIndex(mdl, layout) };
}
const base = await solve(model, structures, atommaps);
const ed = await solve(E.model, E.structures, E.atommaps);
const exFlux = (R, name) => { for (const [m, ex] of R.sim.exchanges) if (R.sim.model.bases[R.sim.model.mets[m].b].n === name) return -R.f.v[ex.rxn] * ex.dir; return 0; };
console.log(`   growth: published ${base.f.objective.toFixed(4)}, edited ${ed.f.objective.toFixed(4)}`);
console.log(`   lactate secretion: published ${exFlux(base, "L-lactate").toFixed(3)}, edited ${exFlux(ed, "L-lactate").toFixed(3)}; acetoin secretion edited ${exFlux(ed, "acetoin").toFixed(3)}`);
check(base.f.status === "optimal" && ed.f.status === "optimal", "both models solve");
console.log(`   removed lactate dehydrogenases: ${ldh.join(", ")}`);
// the new pathway is connected: acetoin can be made and secreted, with carbons from glucose
const prod = await solve(E.model, E.structures, E.atommaps, { objective: E.model.rxns.findIndex((x) => x.id === "USR_EXA") });
const acIdx = E.model.mets.findIndex((m) => m.id === acC);
console.log(`   maximal acetoin secretion ${prod.f.objective.toFixed(3)} mmol/gDW/h`);
check(prod.f.status === "optimal" && prod.f.objective > 0.1, "the new pathway can carry flux to the new exchange");
const pyrIdx = L.met.get(pyr);
check(prod.f.metE[acIdx] > 0.3 && Math.abs(prod.f.metE[acIdx] - prod.f.metE[pyrIdx]) < 0.03, `acetoin carries the label of its precursor pyruvate (${(100 * prod.f.metE[acIdx]).toFixed(1)} % vs ${(100 * prod.f.metE[pyrIdx]).toFixed(1)} %)`);
// lactate: the cell can no longer make it
const lacEx = E.model.rxns.findIndex((x) => x.m.length === 1 && E.model.bases[E.model.mets[x.m[0]].b].n === "L-lactate");
const noLac = await solve(E.model, E.structures, E.atommaps, { objective: lacEx, maximize: model.rxns[lacEx].c[0] < 0 });
const withLac = await solve(model, structures, atommaps, { objective: lacEx, maximize: model.rxns[lacEx].c[0] < 0 });
console.log(`   maximal lactate secretion: published ${withLac.f.objective.toFixed(3)}, edited ${noLac.f.objective.toFixed(3)}`);
check(Math.abs(withLac.f.objective) > 0.1 && Math.abs(noLac.f.objective) < Math.abs(withLac.f.objective), "removing the dehydrogenases lowers the lactate the cell can secrete");

// ---- an atom-mapped reaction SMILES gives the same carbons as the curated map of an existing reaction
const pr = mapFromMappedSmiles("[CH3:1][C:2](=[O:3])[C:4](=[O:5])[OH:6]>>[CH3:1][CH:2]([OH:3])[C:4](=[O:5])[OH:6]",
  [{ loc: 0, side: -1, coef: 1, st: structures[model.bases[model.mets[L.met.get(pyr)].b].id] }, { loc: 1, side: 1, coef: 1, st: structures[model.bases[model.mets[L.met.get(lac)].b].id] }]);
check(pr.mapped === 6 && pr.s.length === 1, `mapped SMILES pyruvate -> lactate: ${pr.mapped}/${pr.total} atoms`);

// stoichiometric copies: RXNMapper's map of 2 pyruvate + H+ -> acetoin + 2 CO2 (from tools/map_reaction.py)
const rx = mapFromMappedSmiles("[CH3:1][C:2](=[O:3])[C:8](=[O:7])[O-:9].[CH3:5][C:4](=[O:6])[C:11](=[O:10])[O-:12].[H+]>>[CH3:1][C:2](=[O:3])[CH:4]([CH3:5])[OH:6].[O:7]=[C:8]=[O:9].[O:10]=[C:11]=[O:12]", parts);
check(rx.mapped === 12 && rx.notes.length === 0 && rx.s.length === 4, `mapped SMILES with two pyruvate copies: ${rx.mapped}/${rx.total} atoms, ${rx.s.length} segments`);
const E2 = applyEdits({ model, structures, atommaps }, [...ops.slice(0, 2), { ...draft, atommap: { k: "user", s: rx.s } }, ...ops.slice(3, 5)]);
check(E2.results.every((x) => x.ok) && !E2.results[2].issues.some((x) => x.key === "chk.mapPartial"), "that map passes the atom-map checks with every product carbon covered");
const prod2 = await solve(E2.model, E2.structures, E2.atommaps, { objective: E2.model.rxns.findIndex((x) => x.id === "USR_EXA") });
const ac2 = E2.model.mets.findIndex((m) => m.id === acC), co2c = L.met.get(co2);
console.log(`   with the RXNMapper map: acetoin ${(100 * prod2.f.metE[ac2]).toFixed(1)} %, pyruvate ${(100 * prod2.f.metE[L.met.get(pyr)]).toFixed(1)} %, CO2[c] ${(100 * prod2.f.metE[co2c]).toFixed(1)} %`);

console.log(bad ? `\n${bad} check(s) FAILED` : "\nall checks passed");
process.exit(bad ? 1 : 0);
