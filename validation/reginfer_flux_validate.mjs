// Does the flux axis of the cross-omic read reproduce biology we already know, and does it refuse to give a
// p-value when it should?
//
// Hypoxia is the clean case. HIF1A rises, it activates the glycolytic enzymes, and the cell ferments: the
// transcriptome and the flux have a textbook expectation to agree on. So:
//
//   1. GSE101526 (Galbraith et al. 2017, HCT116 hypoxia vs normoxia) gives the transcriptome.
//   2. Human-GEM solved with oxygen free, then with oxygen restricted, gives the flux change.
//   3. The enzymes HIF1A is inferred to drive should sit on reactions that carry MORE flux under hypoxia.
//
// And the control that matters as much as the result: because that flux change is imposed on the model
// rather than measured, the joint test must REFUSE to report a p-value. A number there would be the model
// agreeing with itself.
//
//   node tools/validation/reginfer_flux_validate.mjs
import fs from "node:fs";
import zlib from "node:zlib";
const R = new URL("../../", import.meta.url).pathname;
const { default: loadHighs } = await import(R + "vendor/highs/highs.mjs");
const { FbaProblem } = await import(R + "js/core/fba.js");
const { presetMedium } = await import(R + "js/core/medium.js");
const { measuredPreset, measuredNgam } = await import(R + "js/core/measured.js");
const { Simulator } = await import(R + "js/sim.js");
const { ModelIndex } = await import(R + "js/modelindex.js");
const { parseRnaDiffTable, inferRegulatorsMultiOmic, jointFluxTest, regulonRoutes } = await import(R + "js/core/reginfer.js");
globalThis.atob ??= (s) => Buffer.from(s, "base64").toString("binary");

const rd = (f) => JSON.parse(fs.readFileSync(R + "data/human-gem/" + f));
const model = rd("model.json"), structures = rd("structures.json"), atommaps = rd("atommaps.json"), layout = rd("layout.json");
const regulons = rd("regulons.json");
let bad = 0;
const ok = (c, msg) => { console.log((c ? "ok  " : "FAIL") + " " + msg); if (!c) bad++; };

// ---------------------------------------------------------------- the transcriptome
const txt = zlib.gunzipSync(fs.readFileSync(R + "tools/cache/geo_human/GSE101526_counts.txt.gz")).toString("utf8");
const lines = txt.split(/\r?\n/).filter(Boolean);
const hdr = lines[0].split("\t");
const idxNorm = hdr.map((h, i) => (/HCT116wt_norm/i.test(h) ? i : -1)).filter((i) => i > 0);
const idxHyp = hdr.map((h, i) => (/HCT116wt_hypox/i.test(h) ? i : -1)).filter((i) => i > 0);
const libSize = hdr.map(() => 0);
const rawRows = [];
for (let i = 1; i < lines.length; i++) {
  const c = lines[i].split("\t");
  if (!c[0] || c[0].startsWith("__")) continue;
  const vals = c.slice(1).map(Number);
  vals.forEach((v, k) => { if (Number.isFinite(v)) libSize[k + 1] += v; });
  rawRows.push([c[0], vals]);
}
const rows = [], vars = [];
for (const [sym, vals] of rawRows) {
  const cpm = (idx) => idx.map((i) => Math.log2((vals[i - 1] / libSize[i]) * 1e6 + 0.5));
  const gN = cpm(idxNorm), gH = cpm(idxHyp);
  const mean = (a) => a.reduce((s, x) => s + x, 0) / a.length;
  const mN = mean(gN), mH = mean(gH);
  if (mN < 0 && mH < 0) continue;
  const s2 = (gN.reduce((s, x) => s + (x - mN) ** 2, 0) + gH.reduce((s, x) => s + (x - mH) ** 2, 0)) / (gN.length + gH.length - 2);
  rows.push({ sym, lfc: mH - mN, s2, nN: gN.length, nH: gH.length });
  vars.push(s2);
}
vars.sort((a, b) => a - b);
const s02 = vars[Math.floor(vars.length / 2)] || 0.1, d0 = 4;
const tsv = ["gene\tlog2FC\tstat"];
for (const r of rows) {
  const df = r.nN + r.nH - 2;
  const modVar = (d0 * s02 + df * r.s2) / (d0 + df);
  tsv.push(`${r.sym}\t${r.lfc.toFixed(4)}\t${(r.lfc / Math.sqrt(modVar * (1 / r.nN + 1 / r.nH))).toFixed(4)}`);
}
const parsed = parseRnaDiffTable(tsv.join("\n"));
console.log(`GSE101526: ${parsed.map.size} expressed genes, ${idxHyp.length} hypoxic vs ${idxNorm.length} normoxic samples\n`);

// ---------------------------------------------------------------- the flux change
const highs = await loadHighs({ print: () => {}, printErr: () => {} });
let problem = null;
const backend = async (m) => { if (m.type === "init") { problem = new FbaProblem(highs, m.net); return true; } return problem.solve(m.params); };
const sim = new Simulator(model, structures, atommaps, backend), index = new ModelIndex(model, layout);
const base = presetMedium(model, "defaultSerum").map((it) => ({ ...it, label: null, fraction: 1 }));
const o2 = model.rxns.findIndex((r) => r.id === "MAR09048");       // O2 exchange
const cfg = { objective: index.biomassRxn, maximize: true, parsimonious: true };

async function solve(o2conc) {
  // o2conc === null leaves the medium as it is (normoxia); a number restricts oxygen
  const supplies = o2conc === null ? base
    : base.map((it) => (model.bases[model.mets[it.met].b].n === "O2" ? { ...it, conc: o2conc } : it));
  const res = await sim.runSteady(cfg, supplies, structures);
  const f = res.frames[0];
  return /optimal/i.test(f.status) ? f : null;
}
// The oxygen cap must actually bite: at the medium's Vmax 2.5 and Km 0.005, a concentration of 0.004 still
// allows 1.11 mmol/gDW/h while the cell only wants 0.448, so nothing would change. 0.0005 caps it at 0.23.
const fNorm = await solve(null), fHyp = await solve(0.0005);
ok(!!fNorm && !!fHyp, `both states solve (normoxia growth ${fNorm?.objective.toPrecision(3)}, hypoxia ${fHyp?.objective.toPrecision(3)})`);
if (!fNorm || !fHyp) process.exit(1);
const upO2 = (f) => (o2 >= 0 ? -f.v[o2] : NaN);
console.log(`  oxygen uptake: ${upO2(fNorm).toPrecision(3)} -> ${upO2(fHyp).toPrecision(3)} mmol/gDW/h`);
ok(upO2(fHyp) < upO2(fNorm), "the hypoxic state really does take up less oxygen");

// ---------------------------------------------------------------- the read, with and without the flux
const noFlux = inferRegulatorsMultiOmic(model, null, regulons, null, null, parsed.map, { minTargets: 5 });
const rank = (res, tf) => res.tfs.findIndex((x) => x.tf === tf) + 1;
ok(rank(noFlux, "HIF1A") === 1, `HIF1A is still rank ${rank(noFlux, "HIF1A")} of ${noFlux.tfs.length} from the transcriptome alone`);

const withFlux = inferRegulatorsMultiOmic(model, null, regulons, null, null, parsed.map, {
  minTargets: 5, flux: { v: fHyp.v, vRef: fNorm.v, source: "contrast" },
});
const hif = withFlux.tfs.find((x) => x.tf === "HIF1A");
ok(rank(withFlux, "HIF1A") === rank(noFlux, "HIF1A"), "adding the flux does not disturb the ranking, which is scored on the transcriptome");
console.log(`\nHIF1A: ${hif.metGenes.length} enzymes confirmed by mRNA, ${hif.metRxns.size} reactions`);
console.log(`  against the flux: ${hif.fluxAgree} agree, ${hif.fluxDisagree} disagree, ${hif.fluxIdle} did not move`);
ok(hif.fluxAgree + hif.fluxDisagree > 0, "some of HIF1A's enzymes sit on reactions whose flux moved");
ok(hif.fluxAgree > hif.fluxDisagree,
  `more of them move the way HIF1A predicts than against it (${hif.fluxAgree} vs ${hif.fluxDisagree}) -- hypoxia should push HIF1A's targets up`);
for (const g of hif.metGenes.filter((x) => x.fluxAgree).slice(0, 8)) {
  console.log(`    ${g.gene.padEnd(10)} mRNA ${g.log2fc > 0 ? "+" : ""}${g.log2fc.toFixed(2)}  expected ${g.expected > 0 ? "up" : "down"}  flux ${g.fluxDir > 0 ? "up" : g.fluxDir < 0 ? "down" : "--"}  ${g.fluxAgree === 1 ? "agree" : "DISAGREE"}`);
}

// ---------------------------------------------------------------- the control that matters
const jt = jointFluxTest(hif, withFlux, 2000);
ok(jt.p === null && jt.reason === "fluxNotIndependent",
  `the joint test refuses a p-value on an imposed flux change (reason: ${jt.reason})`);
const asIfMeasured = inferRegulatorsMultiOmic(model, null, regulons, null, null, parsed.map, {
  minTargets: 5, flux: { v: fHyp.v, vRef: fNorm.v, source: "fitted" },
});
const jt2 = jointFluxTest(asIfMeasured.tfs.find((x) => x.tf === "HIF1A"), asIfMeasured, 2000);
console.log(`\n  (labelled as measured instead, the same data would give p = ${jt2.p === null ? "(" + jt2.reason + ")" : jt2.p.toFixed(4)};`);
console.log(`   that is exactly why the label is not the user's to invent.)`);

// ---------------------------------------------------------------- routes
const { routes, scattered } = regulonRoutes(hif.metRxns, model, index.currency, 2);
console.log(`\nHIF1A's reactions form ${routes.length} contiguous route(s) and ${scattered.length} scattered steps`);
if (routes.length) console.log(`  largest: ${routes[0].length} reactions -- ${routes[0].slice(0, 6).map((j) => model.rxns[j].id).join(", ")}`);
ok(routes.length > 0, "the regulon is not entirely scattered");

console.log(bad ? `\n${bad} check(s) FAILED` : "\nall passed");
process.exit(bad ? 1 : 0);
