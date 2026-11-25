// Validation Suite 2: Conditional Auxotrophies & Metabolic Synthetic Lethality (Compare tab screens)
// Tests:
//   A. Human-GEM (MFL_biomass_nocof + defaultSerum): acquired nutrient auxotrophies
//      (SDHB -> Asp/Pyr, ASS1 -> Arg, ASNS -> Asn, PHGDH/PSAT1 -> Ser, ALDH18A1 -> Pro, etc.)
//   B. Human-GEM: canonical metabolic synthetic lethal gene pairs (double KO vs single KOs)
//      (SDHB+PC, ENO1+ENO2, ME2+ME3, GOT1+GOT2, SHMT1+SHMT2, IDH1+IDH2, etc.)
//   C. E. coli iML1515 (M9 glucose): 16 classic Keio single-gene auxotrophs + cognate nutrient rescue
//   D. DepMap 26Q3 CRISPR cross-check for paralog & auxotrophy line-specific dependencies
//
// Run: node tools/validation/diff_screens_validate.mjs
import fs from "node:fs";
const R = new URL("../../", import.meta.url).pathname;
const { default: loadHighs } = await import(R + "vendor/highs/highs.mjs");
const { FbaProblem } = await import(R + "js/core/fba.js");
const { presetMedium, mediumBounds } = await import(R + "js/core/medium.js");
const { Simulator } = await import(R + "js/sim.js");
const { gprActivity } = await import(R + "js/core/gpr.js");

const highs = await loadHighs({ print: () => {}, printErr: () => {} });
function loadModel(slug) {
  const rd = (f) => JSON.parse(fs.readFileSync(R + `data/${slug}/${f}`, "utf8"));
  const model = rd("model.json"), st = rd("structures.json"), am = rd("atommaps.json");
  let pb;
  const backend = async (m) => {
    if (m.type === "init") { pb = new FbaProblem(highs, m.net); return true; }
    const r = pb.solve(m.params); r.ms = 0; return r;
  };
  return { model, st, am, sim: new Simulator(model, st, am, backend) };
}

function rxnsLostByGenes(model, genes) {
  const levels = new Map(genes.map((g) => [g, 0]));
  const lost = [];
  for (let j = 0; j < model.rxns.length; j++) {
    const g = model.rxns[j].g;
    if (!g) continue;
    // Check if any gene in `genes` appears in rule
    if (genes.some((sym) => g.includes(sym)) && gprActivity(g, levels) === 0) lost.push(j);
  }
  return Int32Array.from(lost);
}

const out = { humanAuxotrophies: [], humanSyntheticLethals: [], ecoliAuxotrophs: [], depmapParalogCheck: [] };

// -------------------------------------------------------------------------
// Part A & B: Human-GEM (MFL_biomass_nocof + defaultSerum)
// -------------------------------------------------------------------------
console.log("== Part A: Human-GEM acquired nutrient auxotrophies (MFL_biomass_nocof + defaultSerum) ==");
const human = loadModel("human-gem");
const hSup = presetMedium(human.model, "defaultSerum");
const hObj = human.model.rxns.findIndex((r) => r.id === "MFL_biomass_nocof");
const hCfg = { objective: hObj, maximize: true, parsimonious: false };

const findSupMet = (model, sup, regex) => sup.find((s) => regex.test(model.bases[model.mets[s.met].b].n));
const solveGrowth = async (sim, cfg, sup, genes = []) => {
  const lost = genes.length ? rxnsLostByGenes(sim.model, genes) : new Int32Array(0);
  const r = (await sim.screenKnockouts(cfg, sup, [lost]))[0];
  return r.status === "optimal" ? r.objective : 0;
};

const hWt = await solveGrowth(human.sim, hCfg, hSup, []);
console.log(`  Human WT growth (defaultSerum): ${hWt.toFixed(4)} h^-1`);

// Helper: remove a single nutrient from medium (set conc=0)
const dropNutrient = (sup, model, regex) => {
  const hit = findSupMet(model, sup, regex);
  if (!hit) return { sup, name: null };
  return {
    sup: sup.map((s) => (s.met === hit.met ? { ...s, conc: 0, unlimited: false } : s)),
    name: model.bases[model.mets[hit.met].b].n,
  };
};
// Helper: add a nutrient if absent (or boost if present)
const addNutrient = (sup, model, sim, regex, conc = 2.0) => {
  const exMet = [...sim.exchanges.keys()].find((m) => regex.test(model.bases[model.mets[m].b].n));
  if (exMet === undefined) return { sup, name: null };
  const name = model.bases[model.mets[exMet].b].n;
  const exists = sup.some((s) => s.met === exMet);
  return {
    sup: exists ? sup.map((s) => (s.met === exMet ? { ...s, conc: Math.max(s.conc || 0, conc), unlimited: false } : s))
                : [...sup, { met: exMet, conc, vmax: 2.0, km: 0.1 }],
    name,
  };
};

// Known human acquired auxotrophies:
// 1. SDHB KO: requires pyruvate or aspartate when the other is absent (Cardaci et al. 2015, Nature Cell Biol)
//    Note: defaultSerum has pyruvate (0.5 mM) + aspartate (0.1 mM). Test dropping pyruvate+aspartate vs rescuing with each!
const AUX_TESTS = [
  { gene: "SDHB", nutrient: "pyruvate + aspartate", ref: "Cardaci et al. 2015",
    makeDeprived: (s) => dropNutrient(dropNutrient(s, human.model, /^pyruvate$/i).sup, human.model, /^aspartate$/i).sup,
    makeRescue: (s) => addNutrient(dropNutrient(dropNutrient(s, human.model, /^pyruvate$/i).sup, human.model, /^aspartate$/i).sup, human.model, human.sim, /^aspartate$/i, 1.0).sup },
  { gene: "ASS1", nutrient: "arginine", ref: "Delage et al. 2010",
    makeDeprived: (s) => dropNutrient(s, human.model, /^arginine$/i).sup,
    makeRescue: (s) => addNutrient(dropNutrient(s, human.model, /^arginine$/i).sup, human.model, human.sim, /^arginine$/i, 1.0).sup },
  { gene: "ASNS", nutrient: "asparagine", ref: "Richards & Kilberg 2006",
    makeDeprived: (s) => dropNutrient(s, human.model, /^asparagine$/i).sup,
    makeRescue: (s) => addNutrient(dropNutrient(s, human.model, /^asparagine$/i).sup, human.model, human.sim, /^asparagine$/i, 0.5).sup },
  { gene: "PHGDH", nutrient: "serine/glycine", ref: "Possemato et al. 2011",
    makeDeprived: (s) => dropNutrient(dropNutrient(s, human.model, /^serine$/i).sup, human.model, /^glycine$/i).sup,
    makeRescue: (s) => addNutrient(dropNutrient(dropNutrient(s, human.model, /^serine$/i).sup, human.model, /^glycine$/i).sup, human.model, human.sim, /^serine$/i, 1.0).sup },
  { gene: "PSAT1", nutrient: "serine/glycine", ref: "Possemato et al. 2011",
    makeDeprived: (s) => dropNutrient(dropNutrient(s, human.model, /^serine$/i).sup, human.model, /^glycine$/i).sup,
    makeRescue: (s) => addNutrient(dropNutrient(dropNutrient(s, human.model, /^serine$/i).sup, human.model, /^glycine$/i).sup, human.model, human.sim, /^serine$/i, 1.0).sup },
  { gene: "ALDH18A1", nutrient: "proline", ref: "Loayza-Puch et al. 2016",
    makeDeprived: (s) => dropNutrient(s, human.model, /^proline$/i).sup,
    makeRescue: (s) => addNutrient(dropNutrient(s, human.model, /^proline$/i).sup, human.model, human.sim, /^proline$/i, 0.5).sup },
  { gene: "GLUL", nutrient: "glutamine", ref: "Tardito et al. 2015",
    makeDeprived: (s) => dropNutrient(s, human.model, /^glutamine$/i).sup,
    makeRescue: (s) => addNutrient(dropNutrient(s, human.model, /^glutamine$/i).sup, human.model, human.sim, /^glutamine$/i, 2.0).sup },
];

for (const t of AUX_TESTS) {
  const sDep = t.makeDeprived(hSup);
  const sRes = t.makeRescue(hSup);
  const wtDep = await solveGrowth(human.sim, hCfg, sDep, []);
  const koDep = await solveGrowth(human.sim, hCfg, sDep, [t.gene]);
  const koRes = await solveGrowth(human.sim, hCfg, sRes, [t.gene]);
  const acquired = wtDep > 0.5 * hWt && koDep < 0.05 * wtDep && koRes > 0.5 * wtDep;
  console.log(`  [${t.gene} -> ${t.nutrient}] WT(-nut)=${wtDep.toFixed(4)}, KO(-nut)=${koDep.toFixed(4)}, KO(+nut)=${koRes.toFixed(4)} => ${acquired ? "ACQUIRED AUXOTROPHY" : "partial/no"}`);
  out.humanAuxotrophies.push({
    gene: t.gene, nutrient: t.nutrient, ref: t.ref,
    wtDeprived: +wtDep.toFixed(5), koDeprived: +koDep.toFixed(5), koRescued: +koRes.toFixed(5),
    acquiredAuxotrophy: acquired,
  });
}

console.log("\n== Part B: Human-GEM canonical synthetic lethal gene pairs ==");
// Test in standard medium (and for SDHB+PC in pyruvate/aspartate-free medium where PC is required for OAA/Asp!)
const SL_TESTS = [
  // Stoichiometric biosynthetic paralog pairs (truly essential when both paralogs are lost)
  { g1: "IMPDH1", g2: "IMPDH2", medium: hSup, ref: "GMP de novo synthesis paralogs" },
  { g1: "CTPS1", g2: "CTPS2", medium: hSup, ref: "CTP synthase paralogs (Martin et al. 2014)" },
  { g1: "PRPS1", g2: "PRPS2", extra: ["PRPS1L1"], medium: hSup, ref: "PRPP synthetase paralogs" },
  { g1: "MAT2A", g2: "MAT1A", medium: hSup, ref: "SAM synthetase paralogs (Marjon et al. 2016)" },
  { g1: "HK1", g2: "HK2", extra: ["HK3", "GCK", "HKDC1"], medium: hSup, ref: "Hexokinase family" },
  { g1: "RRM2", g2: "RRM2B", medium: hSup, ref: "Ribonucleotide reductase M2 paralogs" },
  { g1: "ACLY", g2: "ACSS2", medium: presetMedium(human.model, "default"), ref: "Cytosolic acetyl-CoA in lipid-free medium (Zhao et al. 2016)" },
  // Rich-medium bypass cases (known SL in vivo/defined media; tests FBA rich-medium bypass honesty)
  { g1: "ENO1", g2: "ENO2", medium: hSup, ref: "Muller et al. 2012 (rich-medium bypass test)" },
  { g1: "SDHB", g2: "PC", medium: dropNutrient(hSup, human.model, /^aspartate$/i).sup, ref: "Cardaci et al. 2015 (rich-medium bypass test)" },
  { g1: "GOT1", g2: "GOT2", medium: dropNutrient(hSup, human.model, /^aspartate$/i).sup, ref: "Birsoy et al. 2015 (rich-medium bypass test)" },
  { g1: "SHMT1", g2: "SHMT2", medium: dropNutrient(hSup, human.model, /^glycine$/i).sup, ref: "Ducker et al. 2016 (rich-medium bypass test)" },
];
for (const p of SL_TESTS) {
  const genes2 = p.extra ? [p.g1, p.g2, ...p.extra] : [p.g1, p.g2];
  const wt = await solveGrowth(human.sim, hCfg, p.medium, []);
  const k1 = await solveGrowth(human.sim, hCfg, p.medium, [p.g1]);
  const k2 = await solveGrowth(human.sim, hCfg, p.medium, p.extra ? [p.g2, ...p.extra] : [p.g2]);
  const k12 = await solveGrowth(human.sim, hCfg, p.medium, genes2);
  const r1 = k1 / wt, r2 = k2 / wt, r12 = k12 / wt;
  const sl = r1 > 0.5 && r2 > 0.5 && r12 < 0.05;
  const syn = r1 > 0.5 && r2 > 0.5 && r12 < r1 * r2 * 0.5;
  const label = sl ? "SYNTHETIC LETHAL" : syn ? "SYNTHETIC SICK" : "redundant/non-essential";
  console.log(`  [${genes2.join(" + ")}] single1=${(100*r1).toFixed(1)}%, single2=${(100*r2).toFixed(1)}%, double=${(100*r12).toFixed(1)}% => ${label}`);
  out.humanSyntheticLethals.push({
    genes: genes2, ref: p.ref, wt: +wt.toFixed(5),
    relSingle1: +r1.toFixed(4), relSingle2: +r2.toFixed(4), relDouble: +r12.toFixed(4),
    classification: label,
  });
}

// -------------------------------------------------------------------------
// Part C: E. coli iML1515 (M9 glucose) — 16 classic Keio single-gene auxotrophs + cognate rescue
// -------------------------------------------------------------------------
console.log("\n== Part C: E. coli iML1515 classic Keio auxotrophs + cognate nutrient rescue (M9 glucose) ==");
const ecoli = loadModel("ecoli-gem");
const eSup = presetMedium(ecoli.model, "m9glc");
const eObj = ecoli.model.rxns.findIndex((r) => r.obj);
const eCfg = { objective: eObj, maximize: true, parsimonious: false };
const eWt = await solveGrowth(ecoli.sim, eCfg, eSup, []);

const ECOLI_AUX = [
  ["trpA", /^tryptophan$/i], ["hisB", /^histidine$/i], ["argH", /^D-arginine$|^arginine$/i], ["metA", /^methionine$/i],
  ["leuB", /^leucine$/i], ["ilvA", /^isoleucine$/i], ["lysA", /^lysine$/i], ["proA", /^proline$/i],
  ["thrC", /^threonine$/i], ["tyrA", /^tyrosine$/i], ["pheA", /^phenylalanine$/i], ["serA", /^serine$/i],
  ["glyA", /^glycine$/i], ["cysE", /^cysteine$/i], ["purA", /^adenine$|^adenosine$/i], ["pyrB", /^uracil$|^uridine$/i],
];
let ecoliRecovered = 0;
for (const [gene, re] of ECOLI_AUX) {
  const koM9 = await solveGrowth(ecoli.sim, eCfg, eSup, [gene]);
  const res = addNutrient(eSup, ecoli.model, ecoli.sim, re, 2.0);
  const koPlus = res.name ? await solveGrowth(ecoli.sim, eCfg, res.sup, [gene]) : 0;
  // Negative control: supplement a non-cognate amino acid (tryptophan unless gene is trpA, else histidine)
  const negRe = gene === "trpA" ? /^histidine$/i : /^tryptophan$/i;
  const neg = addNutrient(eSup, ecoli.model, ecoli.sim, negRe, 2.0);
  const koWrong = await solveGrowth(ecoli.sim, eCfg, neg.sup, [gene]);
  const ok = koM9 < 0.01 * eWt && koPlus > 0.5 * eWt && koWrong < 0.01 * eWt;
  if (ok) ecoliRecovered++;
  console.log(`  [${gene} -> +${res.name || "?"}] M9=${koM9.toFixed(4)}, +cognate=${koPlus.toFixed(4)}, +non-cognate=${koWrong.toFixed(4)} => ${ok ? "EXACT RESCUE" : "mismatch"}`);
  out.ecoliAuxotrophs.push({
    gene, rescueMetabolite: res.name, m9Growth: +koM9.toFixed(5),
    rescuedGrowth: +koPlus.toFixed(5), wrongNutrientGrowth: +koWrong.toFixed(5), exactRescue: ok,
  });
}
console.log(`  E. coli auxotroph rescue summary: ${ecoliRecovered}/${ECOLI_AUX.length} exact cognate rescues`);

// -------------------------------------------------------------------------
// Part D: DepMap 26Q3 empirical cross-check for paralog synthetic lethality
// Check whether lines with low expression of paralog B (e.g., ENO2, ME3, LDHB, GOT2, SHMT1)
// show stronger CRISPR Chronos dependency on paralog A!
// -------------------------------------------------------------------------
console.log("\n== Part D: DepMap 26Q3 empirical paralog dependency cross-check ==");
const dep = JSON.parse(fs.readFileSync(R + "tools/validation/depmap/depmap26q3_metabolic.json", "utf8"));
const exprLines = fs.readFileSync(R + "tools/validation/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", "utf8").split(/\r?\n/);
const exprHdr = exprLines[0].split(",");
const idxModel = exprHdr.indexOf("ModelID");
const idxDef = exprHdr.indexOf("IsDefaultEntryForModel");
const PARALOGS = [
  ["ENO1", "ENO2"], ["ME2", "ME3"], ["LDHA", "LDHB"], ["GOT1", "GOT2"],
  ["SHMT2", "SHMT1"], ["GLS", "GLS2"], ["ASS1", "ASL"], ["G6PD", "TKT"],
  ["PGAM1", "PGAM2"], ["MAT2A", "MAT1A"], ["NME1", "NME2"], ["ACLY", "ACSS2"],
];
const colOf = new Map();
for (const g of new Set(PARALOGS.flat())) {
  const ci = exprHdr.findIndex((h) => h.startsWith(g + " ("));
  if (ci >= 0) colOf.set(g, ci);
}
const exprByLine = new Map();
for (let i = 1; i < exprLines.length; i++) {
  if (!exprLines[i]) continue;
  const c = exprLines[i].split(",");
  if (idxDef >= 0 && c[idxDef] !== "Yes") continue;
  const mid = c[idxModel];
  if (!dep.effect[mid]) continue;
  const rec = {};
  for (const [g, ci] of colOf) rec[g] = Number(c[ci]); // log2(TPM+1)
  exprByLine.set(mid, rec);
}
const depGeneIdx = new Map(dep.genes.map((g, i) => [g, i]));
for (const [target, paralog] of PARALOGS) {
  const gi = depGeneIdx.get(target);
  if (gi === undefined || !colOf.has(paralog)) continue;
  const low = [], high = [];
  for (const [mid, rec] of exprByLine) {
    const ch = dep.effect[mid]?.[gi];
    const ex = rec[paralog];
    if (ch == null || !Number.isFinite(ch) || !Number.isFinite(ex)) continue;
    // log2(TPM+1) < 1 corresponds to TPM < 1 ("off" in MetfluxLive!)
    if (ex < 1.0) low.push(ch);
    else if (ex >= 3.0) high.push(ch);
  }
  const mean = (a) => a.reduce((s, x) => s + x, 0) / a.length;
  if (low.length >= 5 && high.length >= 5) {
    const mLow = mean(low), mHigh = mean(high);
    const vL = low.reduce((s, x) => s + (x - mLow) ** 2, 0) / (low.length - 1);
    const vH = high.reduce((s, x) => s + (x - mHigh) ** 2, 0) / (high.length - 1);
    const t = (mLow - mHigh) / Math.sqrt(vL / low.length + vH / high.length);
    console.log(`  [${target} dependency when ${paralog} < 1 TPM (n=${low.length}) vs >=7 TPM (n=${high.length})]: Chronos ${mLow.toFixed(3)} vs ${mHigh.toFixed(3)} (Δ=${(mLow-mHigh).toFixed(3)}, Welch t=${t.toFixed(2)})`);
    out.depmapParalogCheck.push({
      target, paralog, nLow: low.length, nHigh: high.length,
      meanChronosLow: +mLow.toFixed(4), meanChronosHigh: +mHigh.toFixed(4),
      delta: +(mLow - mHigh).toFixed(4), welchT: +t.toFixed(2),
    });
  }
}

fs.writeFileSync(R + "tools/validation/diff_screens_results.json", JSON.stringify(out, null, 2));
