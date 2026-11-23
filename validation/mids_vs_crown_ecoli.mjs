// Validation 2: E. coli in vitro. Mass isotopomer distributions (MIDs) of proteinogenic amino acid fragments
// measured by GC-MS in 14 parallel labelling experiments (Crown, Long & Antoniewicz, Metab Eng 28:151, 2015;
// E. coli K-12 MG1655, M9 glucose, aerobic batch, exponential growth; Table S3) compared with MetfluxLive on iML1515.
//
// Four flux solutions, the same labelling simulation for all:
//   A  FBA as the app does by default: M9 glucose, maximal growth, pFBA
//   B  A + the exchange rates of the experiment (acetate secretion and O2 uptake per glucose, from Table S4)
//   C  the net fluxes of the authors' 13C-MFA (COMPLETE-MFA fit, Table S4) imposed on the central reactions of
//      iML1515 (95 % confidence intervals as bounds), the rest by maximal growth + pFBA
//   D  C + the exchange (reverse) fluxes of the MFA fit and the dilution of CO2 by unlabelled CO2 of their model
// In every case the predicted MID of each amino acid is diluted by the unlabelled fraction estimated in that
// experiment (biomass from the inoculum; Table S4 "fractional labeling"), and natural isotopes of the whole
// TBDMS fragment are added, because Table S3 is not corrected for natural abundance.
//
//   node tools/validation/mids_vs_crown_ecoli.mjs [n molecules = 4000]
import fs from "node:fs";
const R = new URL("../../", import.meta.url).pathname;
const { default: loadHighs } = await import(R + "vendor/highs/highs.mjs");
const { FbaProblem, packNetwork } = await import(R + "js/core/fba.js");
const { presetMedium, findExchanges, applyMedium } = await import(R + "js/core/medium.js");
const { Tracer } = await import(R + "js/core/tracer.js");
const { IsotopologueSampler } = await import(R + "js/core/isotopologues.js");
const { ModelIndex } = await import(R + "js/modelindex.js");
const { matchPositions } = await import(R + "js/core/tracerlib.js");

const N = +(process.argv[2] || 4000);
const D = JSON.parse(fs.readFileSync(R + "tools/validation/crown2015.json"));
const rd = (f) => JSON.parse(fs.readFileSync(R + "data/ecoli-gem/" + f));
const model = rd("model.json"), structures = rd("structures.json"), maps = rd("atommaps.json"), layout = rd("layout.json");
const highs = await loadHighs({ print: () => {}, printErr: () => {} });
const P = new FbaProblem(highs, packNetwork(model));
const index = new ModelIndex(model, layout);
const rx = (id) => { const j = model.rxns.findIndex((r) => r.id === id); if (j < 0) throw new Error(id); return j; };
const metByName = (n, c) => model.mets.findIndex((m) => m.c === c && model.bases[m.b].n === n);
const bio = index.biomassRxn;

// ---------------------------------------------------------------- fluxes
const ex = findExchanges(model), med = presetMedium(model, "m9glc");
const base = applyMedium(model, ex, med, model.rxns.map((r) => r.lb), model.rxns.map((r) => r.ub));
const GLC = rx("EX_glc__D_e");
const solve = (lb, ub) => P.solve({ lb, ub, objective: bio, maximize: true, parsimonious: true });

const A = solve(base.lb, base.ub);
const glcUp = -A.v[GLC];
console.log(`A  FBA default: ${A.status}, growth ${A.objective.toFixed(3)} h-1, glucose ${glcUp.toFixed(2)} mmol/gDW/h, acetate ${A.v[rx("EX_ac_e")].toFixed(2)}`);

const net = D.net, scale = glcUp / 100; // MFA fluxes are per 100 glucose
const within = (lb, ub, j, lo, hi) => { lb[j] = Math.min(lo, hi); ub[j] = Math.max(lo, hi); };
// B: exchange rates of the experiment, glucose uptake as in A
const lbB = Float64Array.from(base.lb), ubB = Float64Array.from(base.ub);
within(lbB, ubB, GLC, -glcUp, -glcUp);
within(lbB, ubB, rx("EX_ac_e"), net[65][2] * scale, net[65][3] * scale);
within(lbB, ubB, rx("EX_o2_e"), -net[67][3] * scale, -net[67][2] * scale);
const B = solve(lbB, ubB);
console.log(`B  + measured acetate/O2: ${B.status}, growth ${B.objective.toFixed(3)}, acetate ${B.v[rx("EX_ac_e")].toFixed(2)}`);

// C: MFA net fluxes on the central reactions ([MFA flux no., sign] per iML1515 reaction)
const MAP = {
  GLCptspp: [[1, 1]], PGI: [[2, 1]], PFK: [[3, 1]], FBA: [[4, 1]], TPI: [[5, 1]], GAPD: [[6, 1]], PGK: [[6, -1]],
  PGM: [[7, -1]], ENO: [[7, 1]], PYK: [[8, 1]], G6PDH2r: [[9, 1]], PGL: [[9, 1]], GND: [[10, 1]], RPE: [[11, 1]],
  RPI: [[12, -1]], TKT1: [[15, -1]], TKT2: [[14, -1]], TALA: [[17, 1]], EDD: [[18, 1]], EDA: [[19, 1]], PDH: [[20, 1]],
  CS: [[21, 1]], ACONTa: [[22, 1]], ACONTb: [[22, 1]], ICDHyr: [[23, 1]], AKGDH: [[24, 1]], SUCOAS: [[25, -1]],
  SUCDi: [[26, 1]], FUM: [[27, 1]], MDH: [[28, 1]], ICL: [[29, 1]], MALS: [[30, 1]], ME2: [[31, 1]], ME1: [[32, 1]],
  PPC: [[33, 1]], PPCK: [[34, 1]], PTAr: [[35, 1]], ACKr: [[35, -1]], EX_ac_e: [[65, 1]], EX_o2_e: [[67, -1]],
  GHMT2r: [[44, 1]],
};
let C = null, widen = 1;
for (; widen <= 64 && !(C && C.status === "optimal"); widen *= 2) {
  const lb = Float64Array.from(base.lb), ub = Float64Array.from(base.ub);
  within(lb, ub, GLC, -glcUp, -glcUp);
  for (const [id, terms] of Object.entries(MAP)) {
    const [[k, s]] = terms, [, best, lo, hi] = net[k];
    const half = Math.max(best - lo, hi - best) * widen;
    within(lb, ub, rx(id), s * (best - half) * scale, s * (best + half) * scale);
  }
  C = solve(lb, ub);
}
widen /= 2;
console.log(`C  MFA net fluxes imposed: ${C.status}, growth ${C.objective.toFixed(3)}, 95 % intervals widened x${widen}`);
for (const id of ["PGI", "G6PDH2r", "EDD", "PDH", "CS", "PPC", "MALS", "EX_ac_e"])
  console.log(`     ${id.padEnd(8)} A ${(A.v[rx(id)] / scale).toFixed(1).padStart(6)}  B ${(B.v[rx(id)] / scale).toFixed(1).padStart(6)}  C ${(C.v[rx(id)] / scale).toFixed(1).padStart(6)}  MFA ${(MAP[id][0][1] * net[MAP[id][0][0]][1]).toFixed(1).padStart(6)}  (per 100 glucose)`);

// D: exchange fluxes as reverse copies of the reactions (tracer only; net fluxes as in C)
const EXCH = { PGI: 2, FBA: 4, TPI: 5, GAPD: 6, PGK: 6, RPE: 11, RPI: 12, ACONTa: 22, ACONTb: 22, ICDHyr: 23, SUCOAS: 25,
  SUCDi: 26, FUM: 27, MDH: 28, ICL: 29, PTAr: 35, ACKr: 35, GHMT2r: 44 };
// the MFA splits transketolase and transaldolase into half-reactions with a C2/C3 pool; the exchange of a full
// iML1515 reaction is limited by its slower half
const EXCH_TK = { TKT1: Math.min(D.exch[13][1], D.exch[15][1]), TKT2: Math.min(D.exch[13][1], D.exch[14][1]), TALA: Math.min(D.exch[16][1], D.exch[17][1]) };
const INF_EXCH = 1000; // "Inf" in the MFA (unidentifiable, fast): 10 x the glucose uptake
function withExchange(v) {
  const m2 = { ...model, rxns: model.rxns.slice() }, maps2 = { ...maps }, v2 = Array.from(v);
  const add = (id, x) => {
    const j = rx(id), r = model.rxns[j], xs = Math.min(x, INF_EXCH) * scale;
    if (!(xs > 1e-9)) return;
    const n = v[j];
    v2[j] = Math.max(n, 0) + xs;
    const rr = { ...r, id: r.id + "_rev", c: r.c.map((c) => -c) };
    m2.rxns.push(rr); v2.push(Math.max(-n, 0) + xs);
    const mp = maps[r.id];
    if (mp) maps2[rr.id] = { ...mp, s: mp.s.map((sg) => [sg[1], sg[0], sg[2] === -1 ? -1 : sg[2].map((_, i, a) => a[i % 2 ? i - 1 : i + 1])].concat(sg.slice(3))) };
  };
  for (const [id, k] of Object.entries(EXCH)) add(id, D.exch[k][1]);
  for (const [id, x] of Object.entries(EXCH_TK)) add(id, x);
  return { m2, maps2, v2: Float64Array.from(v2) };
}

// ---------------------------------------------------------------- labelling
const glc = metByName("glucose", "e"), st = structures[model.bases[model.mets[glc].b].id];
const pos = matchPositions({
  atoms: [["C1", "C"], ["C2", "C"], ["C3", "C"], ["C4", "C"], ["C5", "C"], ["C6", "C"], ["O1", "O"], ["O2", "O"], ["O3", "O"], ["O4", "O"], ["O5", "O"], ["O6", "O"]],
  bonds: [["C1", "C2"], ["C2", "C3"], ["C3", "C4"], ["C4", "C5"], ["C5", "O5"], ["O5", "C1"], ["C5", "C6"], ["C1", "O1"], ["C2", "O2"], ["C3", "O3"], ["C4", "O4"], ["C6", "O6"]],
}, st);
const lab = (...cs) => { const l = new Float64Array(st.e.length); for (const c of cs) l[pos.get("C" + c)] = 1; return l; };
const U = [1, 2, 3, 4, 5, 6];
// molecule-level mixtures: [weight, carbons labelled]
const TRACERS = [[[1, [1, 2]]], [[1, [2, 3]]], [[1, [4, 5, 6]]], [[1, [2, 3, 4, 5, 6]]], [[0.5, [1]], [0.5, [4, 5, 6]]],
  [[0.5, [1]], [0.5, U]], [[0.8, [1]], [0.2, U]], [[0.205, U], [0.795, []]],
  [[1, [1]]], [[1, [2]]], [[1, [3]]], [[1, [4]]], [[1, [5]]], [[1, [6]]]];

// amino acid fragments (TBDMS, Antoniewicz et al. 2007): carbons of the amino acid in the fragment and fragment type
const AA = { Ala: ["alanine", "C3H7NO2", 2], Gly: ["glycine", "C2H5NO2", 2], Val: ["valine", "C5H11NO2", 2],
  Leu: ["leucine", "C6H13NO2", 2], Ile: ["isoleucine", "C6H13NO2", 2], Ser: ["serine", "C3H7NO3", 3],
  Phe: ["phenylalanine", "C9H11NO2", 2], Asp: ["aspartate", "C4H7NO4", 3], Glu: ["glutamate", "C5H9NO4", 3],
  Tyr: ["tyrosine", "C9H11NO3", 3] };
const parse = (f) => { const c = {}; for (const [, e, n] of f.matchAll(/([A-Z][a-z]?)(\d*)/g)) c[e] = (c[e] || 0) + (n ? +n : 1); return c; };
const addF = (a, b, k = 1) => { const c = { ...a }; for (const e in b) c[e] = (c[e] || 0) + k * b[e]; return c; };
const MASS = { C: 12, H: 1, N: 14, O: 16, Si: 28 };
function fragment(name) {
  const [, aa, nom] = /([A-Z][a-z]+)(\d+)/.exec(name), [, formula, nT] = AA[aa];
  const deriv = addF(parse(formula), { C: 6, H: 14, Si: 1 }, nT); // each TBDMS replaces an H by Si(CH3)2C(CH3)3
  const M = Object.entries(deriv).reduce((s, [e, n]) => s + MASS[e] * n, 0);
  let f, carbons;
  if (+nom === M - 57) { f = addF(deriv, { C: 4, H: 9 }, -1); carbons = "all"; }               // [M-57]+: loses tert-butyl
  else if (+nom === M - 85) { f = addF(deriv, { C: 5, H: 9, O: 1 }, -1); carbons = "noC1"; }   // [M-85]+: and CO (C1)
  else if (+nom === 302) { f = { C: 14, H: 32, N: 1, O: 2, Si: 2 }; carbons = "C1C2"; }          // C1-C2 fragment
  else throw new Error("unknown fragment " + name);
  const m = Object.entries(f).reduce((s, [e, n]) => s + MASS[e] * n, 0);
  if (m !== +nom) throw new Error(`${name}: formula mass ${m}`);
  return { aa, carbons, formula: f };
}
// natural isotope abundances (IUPAC): mass shift -> fraction
const NAT = { C: [[1, 0.0107]], H: [[1, 0.000115]], N: [[1, 0.00364]], O: [[1, 0.00038], [2, 0.00205]], Si: [[1, 0.04685], [2, 0.03092]] };
const conv = (a, b) => { const o = new Array(a.length + b.length - 1).fill(0); a.forEach((x, i) => b.forEach((y, j) => { o[i + j] += x * y; })); return o; };
function natural(counts) {
  let d = [1];
  for (const [e, n] of Object.entries(counts)) {
    const iso = NAT[e]; if (!iso || !n) continue;
    const one = new Array(3).fill(0); one[0] = 1 - iso.reduce((s, [, p]) => s + p, 0); for (const [k, p] of iso) one[k] += p;
    for (let i = 0; i < n; i++) d = conv(d, one).slice(0, 14);
  }
  return d;
}
const FR = Object.fromEntries(Object.keys(D.mids).map((f) => [f, fragment(f)]));

/** carbon indices C1, C2, other carbons of an amino acid structure */
function aaCarbons(m) {
  const s = structures[model.bases[model.mets[m].b].id], nb = s.e.map(() => []);
  for (let k = 0; k < s.b.length; k += 4) { nb[s.b[k]].push(s.b[k + 1]); nb[s.b[k + 1]].push(s.b[k]); }
  const isC = (i) => s.e[i] === "C";
  const c2 = s.e.findIndex((e, i) => e === "C" && nb[i].some((j) => s.e[j] === "N") &&
    nb[i].some((j) => isC(j) && nb[j].filter((x) => s.e[x] === "O").length === 2));
  const c1 = nb[c2].find((j) => isC(j) && nb[j].filter((x) => s.e[x] === "O").length === 2);
  return { all: s.e.flatMap((e, i) => (e === "C" ? [i] : [])), c1, c2 };
}

function simulate(model_, maps_, v, extraInputs, tracer) {
  const T = new Tracer(model_, structures, maps_), idx = new ModelIndex(model_, layout);
  const inputs = [];
  for (const [mi, e] of ex) {
    const up = v[e.rxn] * e.dir;
    if (up > 1e-9) {
      if (mi !== glc) { inputs.push({ met: mi, flux: up, label: null }); continue; }
      const batches = tracer.map(([w, cs]) => ({ w, label: cs.length ? lab(...cs) : null, frac: 1 }));
      const avg = new Float64Array(st.e.length); for (const b of batches) if (b.label) b.label.forEach((x, i) => { avg[i] += b.w * x; });
      inputs.push({ met: mi, flux: up, label: avg, batches });
    }
  }
  inputs.push(...extraInputs);
  const prog = T.compile(v, inputs), sol = T.advance(prog, null, Infinity);
  const S = new IsotopologueSampler(T, idx, v, inputs, sol.e, new Set(["C"]));
  const out = {};
  for (const [aa, [name]] of Object.entries(AA)) {
    const m = metByName(name, "c"), { all, c1, c2 } = aaCarbons(m);
    const pos = all.map((_, i) => i), o = new Uint8Array(all.length);
    const counts = { all: new Float64Array(all.length + 1), noC1: new Float64Array(all.length), C1C2: new Float64Array(3) };
    for (let k = 0; k < N; k++) {
      S._sample(m, all, pos, o, 0, 0);
      let a = 0, n1 = 0, n12 = 0;
      all.forEach((at, i) => { if (o[i]) { a++; if (at === c1) n1++; if (at === c1 || at === c2) n12++; } });
      counts.all[a]++; counts.noC1[a - n1]++; counts.C1C2[n12]++;
    }
    out[aa] = { counts, nC: all.length };
  }
  return out;
}

function predictFragment(sim, frag, ds) {
  const { aa, carbons, formula } = FR[frag], { counts, nC } = sim[aa];
  const tracked = carbons === "all" ? nC : carbons === "noC1" ? nC - 1 : 2;
  const lab = Array.from(counts[carbons], (x) => x / N);
  const f = D.fraclab[aa]?.[ds + 1] ?? 1;
  // unlabelled fraction of the amino acid (old biomass): no tracer atoms
  const dil = lab.map((p, k) => f * p + (k === 0 ? 1 - f : 0));
  // tracked carbons: k labelled, the others natural; then the rest of the fragment
  let tot = new Array(16).fill(0);
  dil.forEach((p, k) => { if (!(p > 0)) return; const rest = natural({ ...formula, C: formula.C - k }); rest.forEach((q, j) => { if (k + j < 16) tot[k + j] += p * q; }); });
  return tot;
}

const variants = {};
const co2 = metByName("CO2", "c");
for (const [key, sol, exch] of [["A", A, false], ["B", B, false], ["C", C, false], ["D", C, true]]) {
  const res = [];
  for (let ds = 0; ds < 14; ds++) {
    let sim;
    if (exch) {
      const { m2, maps2, v2 } = withExchange(sol.v);
      sim = simulate(m2, maps2, v2, [{ met: co2, flux: net[71][1] * scale, label: null }], TRACERS[ds]);
    } else sim = simulate(model, maps, sol.v, [], TRACERS[ds]);
    for (const frag of Object.keys(D.mids)) {
      const meas = D.mids[frag].map((row) => row[ds]);
      if (meas.some((x) => x === null)) continue;
      const pred = predictFragment(sim, frag, ds).slice(0, meas.length);
      const s = pred.reduce((a, b) => a + b, 0), ms = meas.reduce((a, b) => a + b, 0);
      res.push({ frag, ds, meas: meas.map((x) => (100 * x) / ms), pred: pred.map((x) => (100 * x) / s) });
    }
  }
  variants[key] = res;
}

// ---------------------------------------------------------------- metrics
const summary = {};
for (const [key, res] of Object.entries(variants)) {
  const d = res.flatMap((r) => r.meas.map((m, k) => r.pred[k] - m));
  const mae = d.reduce((s, x) => s + Math.abs(x), 0) / d.length, rms = Math.sqrt(d.reduce((s, x) => s + x * x, 0) / d.length);
  const w1 = d.filter((x) => Math.abs(x) <= 1).length / d.length, w2 = d.filter((x) => Math.abs(x) <= 2).length / d.length;
  const xm = res.flatMap((r) => r.meas), xp = res.flatMap((r) => r.pred);
  const mx = xm.reduce((a, b) => a + b) / xm.length, my = xp.reduce((a, b) => a + b) / xp.length;
  let sxy = 0, sxx = 0, syy = 0; xm.forEach((x, i) => { sxy += (x - mx) * (xp[i] - my); sxx += (x - mx) ** 2; syy += (xp[i] - my) ** 2; });
  summary[key] = { n: d.length, mae, rms, within1: w1, within2: w2, r: sxy / Math.sqrt(sxx * syy) };
  console.log(`${key}: ${res.length} fragment MIDs, ${d.length} values: r = ${summary[key].r.toFixed(3)}, mean |error| ${mae.toFixed(2)} mol%, RMS ${rms.toFixed(2)}, within 1 mol% ${(100 * w1).toFixed(0)} %, within 2 mol% ${(100 * w2).toFixed(0)} %`);
}
// worst fragments in D and per-tracer error in A and D
for (const key of ["A", "D"]) {
  const per = D.tracers.map((t, ds) => { const rr = variants[key].filter((r) => r.ds === ds); const d = rr.flatMap((r) => r.meas.map((m, k) => Math.abs(r.pred[k] - m))); return `${t} ${(d.reduce((a, b) => a + b, 0) / d.length).toFixed(2)}`; });
  console.log(`${key} mean |error| per tracer: ${per.join("; ")}`);
  const worst = variants[key].map((r) => ({ ...r, e: r.meas.reduce((s, m, k) => s + Math.abs(r.pred[k] - m), 0) / r.meas.length })).sort((a, b) => b.e - a.e).slice(0, 6);
  for (const w of worst) console.log(`   worst ${key}: ${w.frag} ${D.tracers[w.ds]}  meas ${w.meas.map((x) => x.toFixed(1)).join(" ")} | pred ${w.pred.map((x) => x.toFixed(1)).join(" ")}`);
}
fs.writeFileSync(R + "tools/validation/mids_vs_crown_ecoli.json", JSON.stringify({ summary, variants, widen }, null, 0));
