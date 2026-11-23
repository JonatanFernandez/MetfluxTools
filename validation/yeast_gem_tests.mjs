// Validation 3: Saccharomyces cerevisiae. The benchmarks that yeast-GEM 9.1.1 runs on every release
// (code/python/yeastgem/model_tests, results in data/testResults/README.md of the release), repeated on
// MetfluxLive's converted model (data/yeast-gem, metabolites renamed to Human-GEM names by convert_foreign.py)
// with MetfluxLive's own LP code (js/core/fba.js on HiGHS) and GPR evaluation (js/core/gpr.js).
// Reaction (r_xxxx) and metabolite (s_xxxx) identifiers are kept by the conversion, so each protocol is applied
// to the same reactions as upstream. Data files: tools/validation/yeast/ (from the yeast-GEM 9.1.1 release, CC BY 4.0).
//
//   1  Growth in chemostats (Österlund et al. 2013, BMC Syst Biol 7:36; data/physiology/chemostatData_Tobias2013.tsv):
//      32 chemostats, N-/C-limited, aerobic/anaerobic. Uptake rates fixed, growth maximised. R² about the identity.
//   2  Gene essentiality (Stanford deletion collection, SGD verified ORFs): Kennedy synthetic complete medium,
//      single-gene deletions, viable if growth >= 1e-6 x wild type.
//   3  Intracellular fluxes under anaerobiosis (13C-MFA: Jouhten 2008, Celton 2012, Nissen 1997, Wasylenko 2014;
//      data/physiology/flux_data_anaerobic.tsv): glucose uptake fixed, growth maximised, fluxes per 100 glucose,
//      fold error. Plus the anaerobic exchanges of Sjöberg et al. 2024 at 23 mmol/gDW/h glucose.
//
// Upstream protocols use plain FBA (cobrapy optimize); so does this script. Test 3 is repeated with pFBA, the
// default of the app. The biomass rescaling of upstream (raven_toolbox.biomass.scale_biomass) is reimplemented
// below with the same mass rules (ids.yml); the H+ rebalancing of the pseudoreactions after rescaling is left out.
//
//   node tools/validation/yeast_gem_tests.mjs            -> prints a summary, writes yeast_gem_tests.json
import fs from "node:fs";
const R = new URL("../../", import.meta.url).pathname, V = R + "tools/validation/yeast/";
const { default: loadHighs } = await import(R + "vendor/highs/highs.mjs");
const { FbaProblem, packNetwork } = await import(R + "js/core/fba.js");
const { gprActivity } = await import(R + "js/core/gpr.js");

const base = JSON.parse(fs.readFileSync(R + "data/yeast-gem/model.json"));
const highs = await loadHighs({ print: () => {}, printErr: () => {} });
const tsv = (f) => fs.readFileSync(V + f, "utf8").trim().split("\n").map((l) => l.split("\t"));
const rIdx = new Map(base.rxns.map((r, j) => [r.id, j])), mIdx = new Map(base.mets.map((m, i) => [m.id, i]));
const rx = (id) => { const j = rIdx.get(id); if (j === undefined) throw new Error("no reaction " + id); return j; };
const mt = (id) => { const i = mIdx.get(id); if (i === undefined) throw new Error("no metabolite " + id); return i; };
const GLC = rx("r_1714"), O2 = rx("r_1992"), NH3 = rx("r_1654"), GROWTH = rx("r_2111");
const exchanges = base.rxns.map((r, j) => (r.m.length === 1 && base.mets[r.m[0]].c === "e" ? j : -1)).filter((j) => j >= 0);

// ------------------------------------------------------------- model variants
const clone = (m) => ({ ...m, rxns: m.rxns.map((r) => ({ ...r, m: [...r.m], c: [...r.c] })) });
function setCoef(r, met, value) {
  const k = r.m.indexOf(met);
  if (k >= 0) { if (value === 0) { r.m.splice(k, 1); r.c.splice(k, 1); } else r.c[k] = value; }
  else if (value !== 0) { r.m.push(met); r.c.push(value); }
}
const coef = (r, met) => { const k = r.m.indexOf(met); return k >= 0 ? r.c[k] : 0; };
// masses (g/mol) from formulas; "R" (the tRNA or polymer rest) is ignored, as in RAVEN
const AW = { C: 12.011, H: 1.008, N: 14.007, O: 15.999, P: 30.974, S: 32.06, Fe: 55.845 };
function mw(f) {
  let s = 0;
  for (const [, el, n] of f.matchAll(/([A-Z][a-z]?)(\d*)/g)) if (AW[el]) s += AW[el] * (n ? +n : 1);
  return s;
}
const H = mt("s_0794"); // H+ [c], not counted as mass
const STRAT = { "protein pseudoreaction": 2.016, "RNA pseudoreaction": 18.015, "carbohydrate pseudoreaction": 0 };
const pseudo = (m, name) => m.rxns.find((r) => r.n === name);
/** g/gDW of a biomass component: substrates x (MW - loss) / 1000. */
function mass(m, name) {
  const r = pseudo(m, name), loss = STRAT[name];
  let s = 0;
  r.m.forEach((i, k) => { if (r.c[k] < 0 && i !== H) s += -r.c[k] * (mw(m.bases[m.mets[i].b].f) - loss) / 1000; });
  return s;
}
/** Rescale a pseudoreaction so the component weighs `target` g/gDW: every coefficient except the component itself
 *  (the product, +1) and H+ is multiplied, so the tRNAs released by the protein pseudoreaction follow the charged
 *  tRNAs it uses. */
function scaleTo(m, name, target) {
  const r = pseudo(m, name), f = target / mass(m, name);
  r.c = r.c.map((c, k) => (r.m[k] === H || (c > 0 && m.bases[m.mets[r.m[k]].b].n === name.replace(" pseudoreaction", "")) ? c : c * f));
}
// upstream scale_biomass(model, "RNA", x, balance_out="carbohydrate"): the carbohydrate takes up the mass difference
function scaleBalanced(m, name, target, other) {
  const before = mass(m, name);
  scaleTo(m, name, target);
  scaleTo(m, other, mass(m, other) + (before - target));
}
// anaerobic.yml of yeast-GEM 9.1.1
const AA = tsv("aminoAcid_Bjorkeroth2020.tsv").slice(1);
function anaerobic(m) {
  setCoef(m.rxns[rx("r_4598")], mt("s_3714"), 0); // heme a out of the cofactor pseudoreaction
  // amino acid ratio of anaerobic growth, then protein back to its previous mass
  const p = pseudo(m, "protein pseudoreaction"), before = mass(m, "protein pseudoreaction");
  for (const [, sub, prod, , , ana] of AA) { setCoef(p, mt(sub), -+ana); setCoef(p, mt(prod), +ana); }
  scaleTo(m, "protein pseudoreaction", before);
  const b = m.rxns[rx("r_4041")];
  setCoef(b, mt("s_0689"), coef(b, mt("s_0689")) + 0.08); // FADH2
  setCoef(b, mt("s_0687"), coef(b, mt("s_0687")) - 0.08); // FAD
  setCoef(b, H, coef(b, H) - 0.16);
  const set = (id, lb, ub) => { const r = m.rxns[rx(id)]; if (lb !== undefined) r.lb = lb; if (ub !== undefined) r.ub = ub; };
  set("r_1992", 0);
  for (const id of ["r_1757", "r_1915", "r_2106", "r_2134", "r_1994", "r_2189", "r_1967", "r_1548"]) set(id, -1000);
  set("r_2137", 0); set("r_0714", 0, 0); set("r_0659", 0, 0);
  return m;
}
// growth.m / growth.py, N-limited chemostats
function nLimited(m) {
  scaleTo(m, "protein pseudoreaction", 0.28);
  scaleBalanced(m, "RNA pseudoreaction", 0.0329, "carbohydrate pseudoreaction");
  for (const id of ["r_0472", "r_0501", "r_0507", "r_0509"]) m.rxns[rx(id)].ub = 1000;
  return m;
}

// ------------------------------------------------------------- solving
const problems = new WeakMap();
function lp(m) { if (!problems.has(m)) problems.set(m, new FbaProblem(highs, packNetwork(m))); return problems.get(m); }
const bounds = (m) => [Float64Array.from(m.rxns, (r) => r.lb), Float64Array.from(m.rxns, (r) => r.ub)];
const solve = (m, lb, ub, parsimonious = false) => lp(m).solve({ lb, ub, objective: GROWTH, maximize: true, parsimonious });

const out = { model: "yeast-GEM 9.1.1 as converted for MetfluxLive (data/yeast-gem)", solver: "HiGHS via js/core/fba.js" };

// ------------------------------------------------------------- 1. chemostats
{
  const rows = tsv("chemostatData_Tobias2013.tsv").slice(1).map((r) => r.map(Number));
  const groups = [["N-limited aerobic", 0, 9, false, true], ["C-limited aerobic", 9, 20, false, false], ["C-limited anaerobic", 20, 26, true, false], ["N-limited anaerobic", 26, 32, true, true]];
  const obs = [], pred = [], detail = [];
  for (const [label, a, b, ana, nlim] of groups) {
    let m = clone(base);
    if (ana) m = anaerobic(m);
    if (nlim) m = nLimited(m);
    for (const row of rows.slice(a, b)) {
      const [lb, ub] = bounds(m);
      [GLC, O2, NH3].forEach((j, k) => { const x = row[k]; if (Math.abs(x) === 1000) lb[j] = -x; else lb[j] = ub[j] = -x; });
      const s = solve(m, lb, ub);
      const g = s.status === "optimal" ? Math.abs(s.v[GROWTH]) : 0;
      obs.push(row[3]); pred.push(g); detail.push({ condition: label, glucose: row[0], o2: row[1], nh3: row[2], measured: row[3], predicted: +g.toFixed(4), status: s.status });
    }
  }
  const mean = obs.reduce((s, x) => s + x, 0) / obs.length;
  const r2 = 1 - obs.reduce((s, x, i) => s + (x - pred[i]) ** 2, 0) / obs.reduce((s, x) => s + (x - mean) ** 2, 0);
  const byCond = {};
  for (const d of detail) (byCond[d.condition] ||= []).push(d);
  out.growth = { r2: +r2.toFixed(4), upstream_r2: 0.9019, n: obs.length, perCondition: Object.fromEntries(Object.entries(byCond).map(([k, v]) => [k, +(v.reduce((s, d) => s + Math.abs(d.predicted - d.measured) / d.measured, 0) / v.length).toFixed(3)])), detail };
  console.log(`1 chemostat growth: R² = ${r2.toFixed(4)} (yeast-GEM 9.1.1: 0.9019), n = ${obs.length}`);
  for (const [k, v] of Object.entries(out.growth.perCondition)) console.log(`    ${k}: mean relative error ${(100 * v).toFixed(1)} %`);
}

// ------------------------------------------------------------- 2. gene essentiality
{
  const m = clone(base);
  for (const j of exchanges) { m.rxns[j].lb = 0; m.rxns[j].ub = 1000; }
  const lbOf = (id, v) => { const j = rIdx.get(id); if (j !== undefined) m.rxns[j].lb = v; };
  ["r_1604", "r_1639", "r_1873", "r_1879", "r_1880", "r_1881", "r_1671", "r_1883", "r_1757", "r_1891", "r_1889", "r_1810", "r_1993", "r_1893", "r_1897", "r_1947", "r_1899", "r_1900", "r_1902", "r_1967", "r_1903", "r_1548", "r_1904", "r_2028", "r_2038", "r_1906", "r_2067", "r_1911", "r_1912", "r_1913", "r_2090", "r_1914", "r_2106"].forEach((id) => lbOf(id, -0.5));
  lbOf("r_1714", -20);
  ["r_1672", "r_1654", "r_1992", "r_2005", "r_2060", "r_1861", "r_1832", "r_2100", "r_4593", "r_4595", "r_4596", "r_4597", "r_2049", "r_4594", "r_4600", "r_2020"].forEach((id) => lbOf(id, -1000));
  const [lb0, ub0] = bounds(m);
  const wt = solve(m, lb0, ub0);
  const lines = (f) => new Set(fs.readFileSync(V + f, "utf8").split("\n").map((s) => s.trim()).filter(Boolean));
  const inviable = lines("inviable_orfs.txt"), verified = lines("verified_orfs.txt");
  const viable = new Set();
  const gprRxns = base.rxns.map((r, j) => (r.g ? j : -1)).filter((j) => j >= 0);
  for (const [orf, sym] of base.genes) {
    const key = sym || orf, levels = new Map([[key, 0]]);
    const lb = Float64Array.from(lb0), ub = Float64Array.from(ub0);
    let hit = false;
    for (const j of gprRxns) if (base.rxns[j].g.split(/[\s()]+/).includes(key) && gprActivity(base.rxns[j].g, levels) === 0) { lb[j] = ub[j] = 0; hit = true; }
    if (!hit) { viable.add(orf); continue; }
    const s = solve(m, lb, ub);
    if (s.status === "optimal" && s.objective / wt.objective >= 1e-6) viable.add(orf);
  }
  const genes = new Set(base.genes.map((g) => g[0]));
  const inter = (a, b) => new Set([...a].filter((x) => b.has(x)));
  const expInv = inter(inter(genes, inviable), verified), expVia = inter(new Set([...genes].filter((g) => !inviable.has(g))), verified);
  const modVia = inter(viable, verified), modInv = inter(new Set([...genes].filter((g) => !viable.has(g))), verified);
  const tp = inter(expVia, modVia).size, tn = inter(expInv, modInv).size, fp = inter(expInv, modVia).size, fn = inter(expVia, modInv).size;
  const acc = (tp + tn) / (tp + tn + fp + fn), mcc = (tp * tn - fp * fn) / Math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
  // gene by gene against the classification that yeast-GEM publishes
  const cls = (g) => (expVia.has(g) ? (modVia.has(g) ? "TP" : "FN") : expInv.has(g) ? (modInv.has(g) ? "TN" : "FP") : null);
  const up = tsv("essentialGenes.tsv").slice(1);
  let same = 0; const differ = [];
  for (const [g, c] of up) { if (cls(g) === c) same++; else differ.push(`${g} ${c}->${cls(g)}`); }
  out.essentiality = { accuracy: +acc.toFixed(4), mcc: +mcc.toFixed(3), tp, tn, fp, fn, upstream: { accuracy: 0.9015, tp: 933, tn: 65, fp: 94, fn: 15 }, sameClassAsUpstream: `${same}/${up.length}`, differ };
  console.log(`2 gene essentiality: accuracy ${acc.toFixed(4)} (yeast-GEM: 0.9015), MCC ${mcc.toFixed(3)}; non-essential right ${tp}, essential right ${tn}, essential called viable ${fp}, viable called essential ${fn} (yeast-GEM: 933/65/94/15)`);
  console.log(`    same class as yeast-GEM's own run for ${same} of ${up.length} genes${differ.length ? ": differ " + differ.slice(0, 10).join(", ") : ""}`);
}

// ------------------------------------------------------------- 3. anaerobic fluxes (13C-MFA) and exchanges
{
  const m = anaerobic(clone(base));
  const data = tsv("flux_data_anaerobic.tsv"); // name_in, name_out, exp_low, target_flux, experimental_flux, dataset, rxn_id, ...
  const sets = new Map();
  for (const r of data) (sets.get(r[5]) || sets.set(r[5], []).get(r[5])).push(r);
  const run = (parsimonious) => {
    const pairs = [];
    for (const [ds, rows] of sets) {
      const g = rows.find((r) => r[6] === "r_1714");
      if (!g) continue;
      const [lb, ub] = bounds(m);
      lb[GLC] = ub[GLC] = -+g[3];
      const s = solve(m, lb, ub, parsimonious);
      if (s.status !== "optimal" || !s.v[GLC]) continue;
      for (const r of rows) { const j = rIdx.get(r[6]); if (j === undefined) continue; pairs.push({ ds, rxn: r[6], name: r[7], measured: Math.abs(+r[4]), predicted: Math.abs(-100 * s.v[j] / s.v[GLC]) }); }
    }
    const ok = pairs.filter((p) => p.measured > 1e-9 && p.predicted > 1e-9);
    const lr = ok.map((p) => Math.abs(Math.log10(p.predicted / p.measured))).sort((a, b) => a - b);
    const med = lr.length % 2 ? lr[(lr.length - 1) / 2] : (lr[lr.length / 2 - 1] + lr[lr.length / 2]) / 2;
    return { medianFoldError: +(10 ** med).toFixed(4), meanFoldError: +(10 ** (lr.reduce((s, x) => s + x, 0) / lr.length)).toFixed(4), within2fold: +(lr.filter((x) => x < Math.log10(2)).length / lr.length).toFixed(3), unpredicted: pairs.length - ok.length, n: pairs.length, worst: ok.map((p) => ({ ...p, fold: +(10 ** Math.abs(Math.log10(p.predicted / p.measured))).toFixed(2) })).sort((a, b) => b.fold - a.fold).slice(0, 8).map((p) => `${p.ds} ${p.rxn} ${p.name}: measured ${p.measured} predicted ${p.predicted.toFixed(2)} (${p.fold}x)`) };
  };
  const fba = run(false), pfba = run(true);
  // exchanges at 23 mmol/gDW/h glucose (Sjöberg et al. 2024)
  const ex = tsv("exchange_data_anaerobic.tsv").slice(1);
  const [lb, ub] = bounds(m);
  lb[GLC] = ub[GLC] = -+ex[0][2];
  const s = solve(m, lb, ub);
  const exRows = ex.map((r) => { const v = Math.abs(s.v[rx(r[3])]), mval = +r[5], sd = +r[6]; return { exchange: r[4], measured: mval, sd, predicted: +v.toFixed(3), within: Math.abs(v - mval) <= sd }; });
  const mre = exRows.reduce((t, r) => t + Math.abs(r.predicted - r.measured) / r.measured, 0) / exRows.length;
  out.anaerobic = { fluxes_FBA: fba, fluxes_pFBA: pfba, upstream_medianFoldError: 1.0615, exchanges: exRows, exchangeMeanRelativeError: +mre.toFixed(4), upstream_exchangeMRE: 0.0645 };
  console.log(`3 anaerobic 13C-MFA fluxes (${fba.n} values, ${sets.size} data sets): median fold error ${fba.medianFoldError} with FBA (yeast-GEM: 1.0615), ${pfba.medianFoldError} with pFBA; within 2-fold ${(100 * fba.within2fold).toFixed(0)} % / ${(100 * pfba.within2fold).toFixed(0)} %; not predicted ${fba.unpredicted} / ${pfba.unpredicted}`);
  console.log(`  anaerobic exchanges at 23 mmol/gDW/h glucose: mean relative error ${mre.toFixed(4)} (yeast-GEM: 0.0645): ` + exRows.map((r) => `${r.exchange} ${r.measured}±${r.sd} -> ${r.predicted}${r.within ? "" : " (out)"}`).join("; "));
}
fs.writeFileSync(R + "tools/validation/yeast_gem_tests.json", JSON.stringify(out, null, 1));
