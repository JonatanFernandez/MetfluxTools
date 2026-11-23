// Validation 4: S. cerevisiae in vitro. Mass isotopomer distributions (MIDs) of proteinogenic amino acid fragments
// measured by GC-MS (TBDMS) in yeast S288C (BY4947) grown on SD medium (0.5 % glucose, yeast nitrogen base without
// amino acids) with [1-13C]glucose; corrected for natural isotopes by the authors (Fujiwara et al. 2025,
// Metab Eng Commun 20:e00260, Supplementary Data S1, CC BY-NC 4.0; tools/validation/fujiwara2025_S1.json),
// compared with MetfluxLive on yeast-GEM 9.1.1 as converted for the app.
//
// Three flux solutions, the same labelling simulation for all:
//   A  FBA as the app does: YNB glucose medium, glucose uptake fixed to the measured 14.84 mmol/gDW/h, maximal
//      growth, pFBA
//   B  A + the measured secretion of ethanol, glycerol and acetate (± 2 SD)
//   C  B + the net fluxes of the authors' 13C-MFA imposed on the central reactions of yeast-GEM (95 % confidence
//      intervals as bounds where given, otherwise ± 10 %; widened step by step if infeasible)
//   D  C + the exchange (reverse) fluxes of the MFA fit, as reverse copies of the reactions, and the unlabelled CO2
//      that their model takes up
// The labelled fraction is taken as 1 (cells grown from a small inoculum; the authors fit no dilution term), and no
// natural isotopes are added, as the measured MIDs are already corrected.
// For reference, the error of the authors' own best fit (their 13C-MFA model) against the same data is printed.
//
//   node tools/validation/mids_vs_fujiwara_yeast.mjs [n molecules = 10000]
import fs from "node:fs";
const R = new URL("../../", import.meta.url).pathname;
const { default: loadHighs } = await import(R + "vendor/highs/highs.mjs");
const { FbaProblem, packNetwork } = await import(R + "js/core/fba.js");
const { presetMedium, findExchanges, applyMedium } = await import(R + "js/core/medium.js");
const { Tracer } = await import(R + "js/core/tracer.js");
const { IsotopologueSampler } = await import(R + "js/core/isotopologues.js");
const { ModelIndex } = await import(R + "js/modelindex.js");
const { matchPositions } = await import(R + "js/core/tracerlib.js");

const N = +(process.argv[2] || 10000);
const D = JSON.parse(fs.readFileSync(R + "tools/validation/fujiwara2025_S1.json"));
const rd = (f) => JSON.parse(fs.readFileSync(R + "data/yeast-gem/" + f));
const model = rd("model.json"), structures = rd("structures.json"), maps = rd("atommaps.json"), layout = rd("layout.json");
const highs = await loadHighs({ print: () => {}, printErr: () => {} });
const P = new FbaProblem(highs, packNetwork(model));
const index = new ModelIndex(model, layout);
const rx = (id) => { const j = model.rxns.findIndex((r) => r.id === id); if (j < 0) throw new Error(id); return j; };
const metByName = (n, c) => model.mets.findIndex((m) => m.c === c && model.bases[m.b].n === n);
const GROWTH = model.rxns.findIndex((r) => r.obj);
const F = D.fluxes;

// ---------------------------------------------------------------- fluxes
const ex = findExchanges(model), med = presetMedium(model, "ynb");
const base = applyMedium(model, ex, med, model.rxns.map((r) => r.lb), model.rxns.map((r) => r.ub));
const GLC = rx("r_1714");
const within = (lb, ub, j, lo, hi) => { lb[j] = Math.min(lo, hi); ub[j] = Math.max(lo, hi); };
const solve = (lb, ub) => P.solve({ lb, ub, objective: GROWTH, maximize: true, parsimonious: true });
const fresh = () => [Float64Array.from(base.lb), Float64Array.from(base.ub)];
const glcUp = F.Glcin.measured;

let [lb, ub] = fresh();
within(lb, ub, GLC, -glcUp, -glcUp);
const A = solve(lb, ub);
const exf = (id, v) => v.v[rx(id)];
const show = (tag, s) => console.log(`${tag}: ${s.status}, growth ${s.objective.toFixed(3)} h-1 (measured ${F.Biomass.measured}), ethanol ${exf("r_1761", s).toFixed(2)} (measured ${F.EtOHex.measured.toFixed(2)}), glycerol ${exf("r_1808", s).toFixed(2)} (${F.Glycerolex.measured.toFixed(2)}), acetate ${exf("r_1634", s).toFixed(2)} (${F.Acetateex.measured.toFixed(2)}), O2 ${(-exf("r_1992", s)).toFixed(2)}`);
show("A  FBA default", A);

const SECR = { r_1761: "EtOHex", r_1808: "Glycerolex", r_1634: "Acetateex" };
const lbB = Float64Array.from(lb), ubB = Float64Array.from(ub);
for (const [id, k] of Object.entries(SECR)) within(lbB, ubB, rx(id), F[k].measured - 2 * F[k].sd, F[k].measured + 2 * F[k].sd);
const B = solve(lbB, ubB);
show("B  + measured secretion", B);

// MFA reaction (or net flux) -> yeast-GEM reaction, same direction
const MAP = { oxPPP: "r_0466", PGI: "r_0467", r4: "r_0450", r13: "r_0962", r15: "r_0958", r16: "r_2034", r19: "r_0959",
  r20: "r_0961", r21: "r_0300", r22: "r_0658", r23: "r_0832", r28: "r_0718", r32: "r_0112", MDH: "r_0713" };
let C = null, widen = 1;
for (; widen <= 64 && !(C && C.status === "optimal"); widen *= 2) {
  const l = Float64Array.from(lbB), u = Float64Array.from(ubB);
  for (const [k, id] of Object.entries(MAP)) {
    const f = F[k], lo0 = f.lo ?? f.best - 0.1 * Math.abs(f.best), hi0 = f.hi ?? f.best + 0.1 * Math.abs(f.best);
    const half = Math.max(f.best - lo0, hi0 - f.best) * widen;
    within(l, u, rx(id), f.best - half, f.best + half);
  }
  within(l, u, rx("r_0719"), 0, 0); // malic enzyme: the MFA has one, mapped to the NAD isoform
  C = solve(l, u);
}
widen /= 2;
show(`C  + MFA central fluxes (95 % intervals widened x${widen})`, C);
for (const [k, id] of Object.entries(MAP)) console.log(`     ${k.padEnd(6)} ${id}  A ${A.v[rx(id)].toFixed(2).padStart(7)}  B ${B.v[rx(id)].toFixed(2).padStart(7)}  C ${C.v[rx(id)].toFixed(2).padStart(7)}  MFA ${F[k].best.toFixed(2).padStart(7)}`);

// D: exchange fluxes of the MFA (min of the two directions) on the matching reversible reactions of yeast-GEM
const xch = (a, b) => Math.min(F[a].best, F[b].best);
const EXCH = { r_0467: xch("r2", "r3"), r_1054: xch("r5", "r6"), r_0486: xch("r7", "r8"), r_0892: xch("r7", "r8"),
  r_0893: xch("r11", "r12"), r_0366: xch("r11", "r12"), r_0451: xch("r24", "r25"), r_1021: xch("r24", "r25"),
  r_0713: xch("r26", "r27"), r_0982: xch("r33", "r34"), r_0984: xch("r35", "r36"), r_1048: xch("r39", "r40"),
  r_1049: xch("r41", "r42"), r_1050: xch("r37", "r38") };
const CAP = 10 * glcUp; // exchange fluxes the MFA leaves practically free (> 1000): 10 x the glucose uptake
function withExchange(v) {
  const m2 = { ...model, rxns: model.rxns.slice() }, maps2 = { ...maps }, v2 = Array.from(v);
  for (const [id, x] of Object.entries(EXCH)) {
    const j = rx(id), r = model.rxns[j], xs = Math.min(x, CAP);
    if (!(xs > 1e-9)) continue;
    const n = v[j];
    v2[j] = Math.max(n, 0) + xs;
    const rr = { ...r, id: r.id + "_rev", c: r.c.map((c) => -c) };
    m2.rxns.push(rr); v2.push(Math.max(-n, 0) + xs);
    const mp = maps[r.id];
    if (mp) maps2[rr.id] = { ...mp, s: mp.s.map((sg) => [sg[1], sg[0], sg[2] === -1 ? -1 : sg[2].map((_, i, a) => a[i % 2 ? i - 1 : i + 1])].concat(sg.slice(3))) };
  }
  return { m2, maps2, v2: Float64Array.from(v2) };
}

// ---------------------------------------------------------------- labelling
const glc = metByName("glucose", "e"), st = structures[model.bases[model.mets[glc].b].id];
const pos = matchPositions({
  atoms: [["C1", "C"], ["C2", "C"], ["C3", "C"], ["C4", "C"], ["C5", "C"], ["C6", "C"], ["O1", "O"], ["O2", "O"], ["O3", "O"], ["O4", "O"], ["O5", "O"], ["O6", "O"]],
  bonds: [["C1", "C2"], ["C2", "C3"], ["C3", "C4"], ["C4", "C5"], ["C5", "O5"], ["O5", "C1"], ["C5", "C6"], ["C1", "O1"], ["C2", "O2"], ["C3", "O3"], ["C4", "O4"], ["C6", "O6"]],
}, st);
const tracer = new Float64Array(st.e.length); tracer[pos.get("C1")] = 1; // [1-13C]glucose

const AA = { Ala: "alanine", Gly: "glycine", Val: "valine", Leu: "leucine", Ile: "isoleucine", Ser: "serine", Phe: "phenylalanine",
  Asp: "aspartate", Glu: "glutamate", Tyr: "tyrosine", Pro: "proline", Thr: "threonine" };
// TBDMS fragments: [M-57]+ keeps every carbon; [M-85]+ and [M-159]+ lose C1; m/z 302 is C1-C2
const KIND = { "M-57": "all", "M-85": "noC1", "M-159": "noC1", 302: "C1C2" };
const parseFrag = (f) => { const [, aa, k] = /^([A-Z][a-z]+)\[(.+)\]$/.exec(f); return { aa, kind: KIND[k] }; };

/** carbon indices C1, C2, other carbons of an amino acid structure */
function aaCarbons(m) {
  const s = structures[model.bases[model.mets[m].b].id], nb = s.e.map(() => []);
  for (let k = 0; k < s.b.length; k += 4) { nb[s.b[k]].push(s.b[k + 1]); nb[s.b[k + 1]].push(s.b[k]); }
  const isC = (i) => s.e[i] === "C", nO = (j) => nb[j].filter((x) => s.e[x] === "O").length;
  const c2 = s.e.findIndex((e, i) => e === "C" && nb[i].some((j) => s.e[j] === "N") && nb[i].some((j) => isC(j) && nO(j) === 2));
  const c1 = nb[c2].find((j) => isC(j) && nO(j) === 2);
  return { all: s.e.flatMap((e, i) => (e === "C" ? [i] : [])), c1, c2 };
}

function simulate(v, model_ = model, maps_ = maps, extra = []) {
  const T = new Tracer(model_, structures, maps_);
  const inputs = [...extra];
  for (const [mi, e] of ex) {
    const up = v[e.rxn] * e.dir;
    if (up > 1e-9) inputs.push({ met: mi, flux: up, label: mi === glc ? tracer : null });
  }
  const prog = T.compile(v, inputs), sol = T.advance(prog, null, Infinity);
  const S = new IsotopologueSampler(T, model_ === model ? index : new ModelIndex(model_, layout), v, inputs, sol.e, new Set(["C"]));
  const out = {};
  for (const [aa, name] of Object.entries(AA)) {
    const m = metByName(name, "c"), { all, c1, c2 } = aaCarbons(m);
    const p = all.map((_, i) => i), o = new Uint8Array(all.length);
    const counts = { all: new Float64Array(all.length + 1), noC1: new Float64Array(all.length), C1C2: new Float64Array(3) };
    for (let k = 0; k < N; k++) {
      S._sample(m, all, p, o, 0, 0);
      let a = 0, n1 = 0, n12 = 0;
      all.forEach((at, i) => { if (o[i]) { a++; if (at === c1) n1++; if (at === c1 || at === c2) n12++; } });
      counts.all[a]++; counts.noC1[a - n1]++; counts.C1C2[n12]++;
    }
    out[aa] = counts;
  }
  return out;
}

function compare(pred) {
  const res = [];
  for (const [frag, meas] of Object.entries(D.mids)) {
    const p = pred(frag).slice(0, meas.length), ps = p.reduce((a, b) => a + b, 0), ms = meas.reduce((a, b) => a + b, 0);
    res.push({ frag, meas: meas.map((x) => (100 * x) / ms), pred: p.map((x) => (100 * x) / ps) });
  }
  return res;
}
const variants = {};
for (const [key, s, exch] of [["A", A], ["B", B], ["C", C], ["D", C, true]]) {
  let sim;
  if (exch) { const { m2, maps2, v2 } = withExchange(s.v); sim = simulate(v2, m2, maps2, [{ met: metByName("CO2", "c"), flux: F.CO2uptake.best, label: null }]); }
  else sim = simulate(s.v);
  variants[key] = compare((frag) => { const { aa, kind } = parseFrag(frag), c = sim[aa][kind]; const out = Array.from(c, (x) => x / N); while (out.length < 12) out.push(0); return out; });
}
variants.MFA = compare((frag) => D.fitted[frag].concat(new Array(12).fill(0)));

// ---------------------------------------------------------------- metrics
const summary = {};
for (const [key, res] of Object.entries(variants)) {
  const d = res.flatMap((r) => r.meas.map((m, k) => r.pred[k] - m));
  const mae = d.reduce((s, x) => s + Math.abs(x), 0) / d.length;
  const xm = res.flatMap((r) => r.meas), xp = res.flatMap((r) => r.pred);
  const mx = xm.reduce((a, b) => a + b) / xm.length, my = xp.reduce((a, b) => a + b) / xp.length;
  let sxy = 0, sxx = 0, syy = 0; xm.forEach((x, i) => { sxy += (x - mx) * (xp[i] - my); sxx += (x - mx) ** 2; syy += (xp[i] - my) ** 2; });
  summary[key] = { n: d.length, mae: +mae.toFixed(3), r: +(sxy / Math.sqrt(sxx * syy)).toFixed(4), within1: +(d.filter((x) => Math.abs(x) <= 1).length / d.length).toFixed(3), within2: +(d.filter((x) => Math.abs(x) <= 2).length / d.length).toFixed(3) };
  console.log(`${key.padEnd(3)}: ${res.length} fragment MIDs, ${d.length} values: r = ${summary[key].r}, mean |error| ${mae.toFixed(2)} mol%, within 1 mol% ${(100 * summary[key].within1).toFixed(0)} %, within 2 mol% ${(100 * summary[key].within2).toFixed(0)} %`);
}
for (const key of ["A", "D"]) {
  const worst = variants[key].map((r) => ({ ...r, e: r.meas.reduce((s, m, k) => s + Math.abs(r.pred[k] - m), 0) / r.meas.length })).sort((a, b) => b.e - a.e).slice(0, 5);
  for (const w of worst) console.log(`   worst ${key}: ${w.frag}  meas ${w.meas.map((x) => x.toFixed(1)).join(" ")} | pred ${w.pred.map((x) => x.toFixed(1)).join(" ")}`);
}
fs.writeFileSync(R + "tools/validation/mids_vs_fujiwara_yeast.json", JSON.stringify({ summary, widen, variants }, null, 0));
