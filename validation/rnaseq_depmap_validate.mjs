// Validation Suite 3: Browser User RNA-seq Context Builder (js/core/ctxbuild.js + tasks.js)
// validated directly on DepMap 26Q3 matched cell-line RNA-seq (OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv)
// against DepMap 26Q3 CRISPR Chronos gene essentiality (depmap26q3_metabolic.json).
//
// Uses 4-core parallelism (--worker sub-processes) and exact monotonicity pruning (any gene already
// essential in the generic model for the same medium and not switched off in the context remains
// mathematically guaranteed essential when additional reactions are blocked).
//
// Evaluates:
//   1. All 6 culture media (default, defaultSerum, rpmi1640, dmem, plasmax, hplm + line-matched RPMI/DMEM)
//      scored across ALL 730 DepMap 26Q3 cell lines (2,014,472 cell-line x gene evaluations).
//   2. Browser ContextBuilder (57-task rescue) build-time & rescue audit across 50 DepMap lines x 4 TPM cutoffs (0.5, 1, 2, 5).
//   3. Head-to-head CRISPR essentiality on a stratified 12-line panel (RPMI & DMEM lines across lineages):
//      Generic vs. prebuilt HPA contexts.json vs. ContextBuilder(DepMap matched RNA-seq at 0.5 / 1 / 2 / 5 TPM) vs. matched medium.
//
// Run: node tools/validation/rnaseq_depmap_validate.mjs
import fs from "node:fs";
import { spawn } from "node:child_process";
import { performance } from "node:perf_hooks";

const R = new URL("../../", import.meta.url).pathname;
const WORKER_ARG = process.argv.indexOf("--worker");
const WORKER_ROLE = WORKER_ARG >= 0 ? process.argv[WORKER_ARG + 1] : null;
const TMP_DIR = R + "tools/validation/depmap/rnaseq_tmp/";

// -------------------------------------------------------------------------
// Orchestrator mode: spawn 4 parallel worker processes across CPU cores
// -------------------------------------------------------------------------
if (!WORKER_ROLE) {
  fs.mkdirSync(TMP_DIR, { recursive: true });
  const runWorker = (role) => new Promise((resolve, reject) => {
    const p = spawn(process.execPath, [new URL(import.meta.url).pathname, "--worker", role], { stdio: ["ignore", "pipe", "pipe"] });
    let out = "", err = "";
    p.stdout.on("data", (d) => { out += d; process.stdout.write(`[${role}] ${d}`); });
    p.stderr.on("data", (d) => { err += d; process.stderr.write(`[${role}] ${d}`); });
    p.on("close", (code) => (code === 0 ? resolve(out) : reject(new Error(`worker ${role} exited ${code}: ${err}`))));
  });
  const t0 = Date.now();
  // Phase 1 (4 cores in parallel): media screens (2 cores) + 50-line ContextBuilder audit (1 core) + lines0 (1 core)
  await Promise.all([
    runWorker("media1"),   // default, rpmi1640, dmem
    runWorker("media2"),   // plasmax, hplm
    runWorker("audit50"),  // 50 DepMap lines x 4 TPM cutoffs ContextBuilder + task rescue
    runWorker("panel0"),   // panel lines 0..5 (defaultSerum across 0.5/1/2/5 TPM)
  ]);
  // Phase 2 (2 cores in parallel): panel lines 6..11 + line-matched media (needs rpmi1640/dmem generic ess)
  await Promise.all([
    runWorker("panel1"),       // panel lines 6..11 (defaultSerum across 0.5/1/2/5 TPM)
    runWorker("panelMatched"), // all 12 panel lines on line-matched RPMI/DMEM (1 TPM)
  ]);
  // Phase 3: merge & evaluate all metrics cleanly
  await runWorker("merge");
  console.log(`\nSuite 3 completed in ${((Date.now() - t0) / 1000).toFixed(1)} s`);
  process.exit(0);
}

// -------------------------------------------------------------------------
// Shared setup inside worker processes
// -------------------------------------------------------------------------
const { default: loadHighs } = await import(R + "vendor/highs/highs.mjs");
const { FbaProblem } = await import(R + "js/core/fba.js");
const { presetMedium } = await import(R + "js/core/medium.js");
const { measuredPreset, measuredNgam } = await import(R + "js/core/measured.js");
const { Simulator } = await import(R + "js/sim.js");
const { gprActivity } = await import(R + "js/core/gpr.js");
const { blockedReactions, isMitoGene } = await import(R + "js/core/contexts.js");
const { parseTasks, TaskRunner } = await import(R + "js/core/tasks.js");
const { ContextBuilder } = await import(R + "js/core/ctxbuild.js");

const rd = (f) => JSON.parse(fs.readFileSync(R + f, "utf8"));
const model = rd("data/human-gem/model.json");
const dep = rd("tools/validation/depmap/depmap26q3_metabolic.json");
const genericDefSerum = new Set(rd("tools/validation/depmap/pred_generic_nocof+serum_0of1.json").generic.essential);

// Parse Model.csv for RPMI vs DMEM culture medium
const modelCsv = fs.readFileSync(R + "tools/validation/depmap/Model.csv", "utf8").split(/\r?\n/);
const mHdr = modelCsv[0].split(",");
const iMid = mHdr.indexOf("ModelID"), iForm = mHdr.indexOf("FormulationID"), iLin = mHdr.indexOf("OncotreeLineage");
const mediaByLine = new Map(), lineageByLine = new Map();
for (let i = 1; i < modelCsv.length; i++) {
  if (!modelCsv[i]) continue;
  const cols = []; let cur = "", inQ = false;
  for (const ch of modelCsv[i]) {
    if (ch === '"') inQ = !inQ;
    else if (ch === "," && !inQ) { cols.push(cur); cur = ""; }
    else cur += ch;
  }
  cols.push(cur);
  const id = cols[iMid]; if (!id) continue;
  lineageByLine.set(id, cols[iLin] || "");
  if (/\bRPMI\b/i.test(cols[iForm] || "")) mediaByLine.set(id, "rpmi1640");
  else if (/\bDMEM\b/i.test(cols[iForm] || "")) mediaByLine.set(id, "dmem");
}

// Load DepMap matched RNA-seq expression for Human-GEM genes
const modelGeneSyms = new Set(model.genes.map((g) => (Array.isArray(g) ? g[1] || g[0] : g)));
const exprLines = fs.readFileSync(R + "tools/validation/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", "utf8").split(/\r?\n/);
const eHdr = exprLines[0].split(",");
const eMid = eHdr.indexOf("ModelID"), eDef = eHdr.indexOf("IsDefaultEntryForModel");
const geneCols = [];
for (let c = 0; c < eHdr.length; c++) {
  const m = eHdr[c].match(/^([^\s(]+)\s+\(/);
  if (m && modelGeneSyms.has(m[1])) geneCols.push([m[1], c]);
}
const exprByLine = new Map();
for (let i = 1; i < exprLines.length; i++) {
  if (!exprLines[i]) continue;
  const c = exprLines[i].split(",");
  if (eDef >= 0 && c[eDef] !== "Yes") continue;
  const mid = c[eMid];
  if (dep.lines[mid] && dep.effect[mid]) exprByLine.set(mid, c);
}
const candidates = [...exprByLine.keys()].sort((a, b) => {
  const ma = mediaByLine.has(a) ? 0 : 1, mb = mediaByLine.has(b) ? 0 : 1;
  return ma - mb || a.localeCompare(b);
});
// Pick 12 stratified lines (6 RPMI + 6 DMEM across diverse lineages)
const panelLines = [];
const linCnt = new Map();
let nR = 0, nD = 0;
for (const mid of candidates) {
  const m = mediaByLine.get(mid);
  if (!m) continue;
  if (m === "rpmi1640" && nR >= 6) continue;
  if (m === "dmem" && nD >= 6) continue;
  const lin = lineageByLine.get(mid) || "Other";
  if ((linCnt.get(lin) || 0) >= 2) continue;
  panelLines.push(mid);
  linCnt.set(lin, (linCnt.get(lin) || 0) + 1);
  if (m === "rpmi1640") nR++; else nD++;
  if (panelLines.length >= 12) break;
}

function getTpmMap(mid) {
  const row = exprByLine.get(mid);
  const map = new Map();
  for (const [sym, ci] of geneCols) {
    const v = Number(row[ci]);
    if (Number.isFinite(v)) map.set(sym, Math.pow(2, v) - 1);
  }
  return map;
}

// -------------------------------------------------------------------------
// Merge worker (no HiGHS needed)
// -------------------------------------------------------------------------
if (WORKER_ROLE === "merge") {
  const allDepLines = Object.keys(dep.lines).filter((m) => dep.effect[m]);
  const commonEss = new Set(dep.commonEssential || []);
  const m1 = JSON.parse(fs.readFileSync(TMP_DIR + "media1.json", "utf8"));
  const m2 = JSON.parse(fs.readFileSync(TMP_DIR + "media2.json", "utf8"));
  const mediaEss = { defaultSerum: [...genericDefSerum], ...m1, ...m2 };
  const audit50 = JSON.parse(fs.readFileSync(TMP_DIR + "audit50.json", "utf8"));
  const p0 = JSON.parse(fs.readFileSync(TMP_DIR + "panel0.json", "utf8"));
  const p1 = JSON.parse(fs.readFileSync(TMP_DIR + "panel1.json", "utf8"));
  const pMatch = JSON.parse(fs.readFileSync(TMP_DIR + "panelMatched.json", "utf8"));

  // Load cached HPA complex predictions for all 730 lines
  const hpaPredAll = {};
  for (let s = 0; s < 4; s++) {
    Object.assign(hpaPredAll, rd(`tools/validation/depmap/pred_complex_nocof+serum_${s}of4.json`));
  }

  function scoreLineSubset(lineIds, predByLine, chronosCut = -0.5) {
    let tp = 0, fp = 0, tn = 0, fn = 0;
    let tpNC = 0, fpNC = 0, tnNC = 0, fnNC = 0;
    for (const mid of lineIds) {
      const ess = predByLine.get(mid) || new Set();
      const eff = dep.effect[mid];
      if (!eff) continue;
      for (let gi = 0; gi < dep.genes.length; gi++) {
        const ch = eff[gi];
        if (ch == null || !Number.isFinite(ch)) continue;
        const g = dep.genes[gi];
        const exp = ch <= chronosCut, pred = ess.has(g);
        if (exp && pred) tp++; else if (!exp && pred) fp++; else if (!exp && !pred) tn++; else fn++;
        if (!commonEss.has(g)) {
          if (exp && pred) tpNC++; else if (!exp && pred) fpNC++; else if (!exp && !pred) tnNC++; else fnNC++;
        }
      }
    }
    const met = (a, b, c, d) => {
      const prec = a + b ? a / (a + b) : 0, rec = a + d ? a / (a + d) : 0;
      const den = Math.sqrt((a + b) * (a + d) * (c + b) * (c + d));
      return { tp: a, fp: b, tn: c, fn: d, prec: +prec.toFixed(4), rec: +rec.toFixed(4), mcc: +(den ? (a * c - b * d) / den : 0).toFixed(4) };
    };
    return { all: met(tp, fp, tn, fn), nonCommon: met(tpNC, fpNC, tnNC, fnNC) };
  }

  const report = {
    nAllDepLines: allDepLines.length,
    nPanelLines: panelLines.length,
    panelLines: panelLines.map((id) => ({ id, name: dep.lines[id].name, lineage: lineageByLine.get(id), medium: mediaByLine.get(id) })),
    nGenesMatchedInDepMapRna: geneCols.length,
    nModelGenes: modelGeneSyms.size,
    mediumComparison730: [],
    contextBuildAudit50: audit50,
    thresholdSweepPanel12: [],
    hpaVsDepmapPanel12: {},
  };

  console.log(`\n== 1. Medium comparison across ALL ${allDepLines.length} DepMap 26Q3 lines (generic model) ==`);
  for (const med of ["default", "defaultSerum", "rpmi1640", "dmem", "plasmax", "hplm"]) {
    const ess = new Set(mediaEss[med]);
    const pred = new Map(allDepLines.map((m) => [m, ess]));
    const s05 = scoreLineSubset(allDepLines, pred, -0.5), s10 = scoreLineSubset(allDepLines, pred, -1.0);
    console.log(`  [generic | ${med.padEnd(12)}] nEss=${String(ess.size).padStart(3)} | <=-0.5: Prec=${(100*s05.all.prec).toFixed(1)}% Rec=${(100*s05.all.rec).toFixed(1)}% MCC=${s05.all.mcc.toFixed(3)} | <=-1: MCC=${s10.all.mcc.toFixed(3)}`);
    report.mediumComparison730.push({ mode: "generic", medium: med, nEss: ess.size, chronos05: s05, chronos10: s10 });
  }
  {
    const predM = new Map(); let nRpmi = 0, nDmem = 0, nDef = 0;
    for (const mid of allDepLines) {
      const m = mediaByLine.get(mid) || "defaultSerum";
      if (m === "rpmi1640") nRpmi++; else if (m === "dmem") nDmem++; else nDef++;
      predM.set(mid, new Set(mediaEss[m]));
    }
    const s05 = scoreLineSubset(allDepLines, predM, -0.5), s10 = scoreLineSubset(allDepLines, predM, -1.0);
    console.log(`  [generic | matched(${nRpmi} RPMI,${nDmem} DMEM,${nDef} defS)] <=-0.5: Prec=${(100*s05.all.prec).toFixed(1)}% Rec=${(100*s05.all.rec).toFixed(1)}% MCC=${s05.all.mcc.toFixed(3)} | <=-1: MCC=${s10.all.mcc.toFixed(3)}`);
    report.mediumComparison730.push({ mode: "generic_matched_media", nRpmi, nDmem, nDef, chronos05: s05, chronos10: s10 });
  }

  console.log(`\n== 2. Browser ContextBuilder 57-task rescue audit (50 DepMap lines x 4 TPM cutoffs) ==`);
  for (const r of audit50) {
    console.log(`  [thr=${r.thresholdTpm} TPM] avgBuild=${r.avgBuildMs} ms/line, off0=${r.avgOff0}, rescued=${r.avgRescued}, blockedRxns=${r.avgBlocked}, failedTasks=${r.totalFailedTasks}/${r.nLines * 57}`);
  }

  console.log(`\n== 3. Head-to-head DepMap 26Q3 essentiality on ${panelLines.length} stratified lines (6 RPMI + 6 DMEM) ==`);
  const genS05 = scoreLineSubset(panelLines, new Map(panelLines.map((m) => [m, genericDefSerum])), -0.5);
  const genS10 = scoreLineSubset(panelLines, new Map(panelLines.map((m) => [m, genericDefSerum])), -1.0);
  console.log(`  [generic | defaultSerum]              Prec=${(100*genS05.all.prec).toFixed(1)}% Rec=${(100*genS05.all.rec).toFixed(1)}% MCC=${genS05.all.mcc.toFixed(3)} (nonCommon MCC=${genS05.nonCommon.mcc.toFixed(3)}) | <=-1 MCC=${genS10.all.mcc.toFixed(3)}`);
  report.hpaVsDepmapPanel12.generic_defaultSerum = { chronos05: genS05, chronos10: genS10 };

  const hpaMap = new Map(panelLines.map((m) => [m, new Set(hpaPredAll[m]?.essential || [])]));
  const hpaS05 = scoreLineSubset(panelLines, hpaMap, -0.5), hpaS10 = scoreLineSubset(panelLines, hpaMap, -1.0);
  console.log(`  [HPA contexts.json | defaultSerum]    Prec=${(100*hpaS05.all.prec).toFixed(1)}% Rec=${(100*hpaS05.all.rec).toFixed(1)}% MCC=${hpaS05.all.mcc.toFixed(3)} (nonCommon MCC=${hpaS05.nonCommon.mcc.toFixed(3)}) | <=-1 MCC=${hpaS10.all.mcc.toFixed(3)}`);
  report.hpaVsDepmapPanel12.hpa_defaultSerum = { chronos05: hpaS05, chronos10: hpaS10 };

  for (const thr of [0.5, 1.0, 2.0, 5.0]) {
    const merged = { ...(p0[thr] || {}), ...(p1[thr] || {}) };
    const pred = new Map(panelLines.map((m) => [m, new Set(merged[m] || [])]));
    const s05 = scoreLineSubset(panelLines, pred, -0.5), s10 = scoreLineSubset(panelLines, pred, -1.0);
    console.log(`  [DepMap RNA-seq thr=${thr} | defSerum]  Prec=${(100*s05.all.prec).toFixed(1)}% Rec=${(100*s05.all.rec).toFixed(1)}% MCC=${s05.all.mcc.toFixed(3)} (nonCommon MCC=${s05.nonCommon.mcc.toFixed(3)}) | <=-1 MCC=${s10.all.mcc.toFixed(3)}`);
    report.thresholdSweepPanel12.push({ thresholdTpm: thr, chronos05: s05, chronos10: s10 });
  }

  const predMatch = new Map(panelLines.map((m) => [m, new Set(pMatch[m] || [])]));
  const mS05 = scoreLineSubset(panelLines, predMatch, -0.5), mS10 = scoreLineSubset(panelLines, predMatch, -1.0);
  console.log(`  [DepMap RNA-seq 1 TPM | matched med]  Prec=${(100*mS05.all.prec).toFixed(1)}% Rec=${(100*mS05.all.rec).toFixed(1)}% MCC=${mS05.all.mcc.toFixed(3)} (nonCommon MCC=${mS05.nonCommon.mcc.toFixed(3)}) | <=-1 MCC=${mS10.all.mcc.toFixed(3)}`);
  report.hpaVsDepmapPanel12.depmap1Tpm_matchedMedium = { chronos05: mS05, chronos10: mS10 };

  fs.writeFileSync(R + "tools/validation/rnaseq_depmap_results.json", JSON.stringify(report, null, 2));
  process.exit(0);
}

// -------------------------------------------------------------------------
// HiGHS-backed workers
// -------------------------------------------------------------------------
const st = rd("data/human-gem/structures.json");
const am = rd("data/human-gem/atommaps.json");
const highs = await loadHighs({ print: () => {}, printErr: () => {} });
let pb;
const backend = async (m) => {
  if (m.type === "init") { pb = new FbaProblem(highs, m.net); return true; }
  const r = pb.solve(m.params); r.ms = 0; return r;
};
const sim = new Simulator(model, st, am, backend);
const objIdx = model.rxns.findIndex((r) => r.id === "MFL_biomass_nocof");
const meas = { items: measuredPreset(model, sim.exchanges, "nci60").map((x) => ({ met: x.met, rate: x.rate })), tol: 0.2, ngam: measuredNgam("nci60") };
const cfgBase = { objective: objIdx, maximize: true, parsimonious: true, measured: meas };

const geneRxns = new Map();
model.rxns.forEach((r, j) => {
  if (!r.g) return;
  for (const g of new Set(r.g.match(/[^\s()]+/g).filter((x) => !/^(and|or)$/i.test(x)))) {
    if (!geneRxns.has(g)) geneRxns.set(g, []);
    geneRxns.get(g).push(j);
  }
});
const supCache = new Map();
const getSup = (id) => { if (!supCache.has(id)) supCache.set(id, presetMedium(model, id)); return supCache.get(id); };

async function screenWithOff(offSet, mediumId, knownGenericEss = null) {
  const sup = getSup(mediumId);
  const off = new Set([...offSet].filter((g) => !isMitoGene(g)));
  const block = off.size ? blockedReactions(model, off, { complexes: true }) : null;
  const cfg = { ...cfgBase, context: block ? { id: "user", mode: "block", block } : undefined };
  const f = (await sim.runSteady(cfg, sup, st)).frames[0];
  if (f.status !== "optimal" || !(f.objective > 1e-9)) return [];
  const blocked = new Set(block || []);
  const ess = new Set();
  const groups = new Map();
  for (const [sym, rx] of geneRxns) {
    if (off.has(sym)) continue;
    // Exact monotonicity pruning: if sym is already essential in the unconstrained generic model
    // for the SAME medium (and not off in this context), adding blocked reactions cannot rescue it!
    if (knownGenericEss && knownGenericEss.has(sym)) {
      ess.add(sym);
      continue;
    }
    if (!rx.some((j) => Math.abs(f.v[j]) > 1e-9)) continue;
    const levels = new Map([...off].map((g) => [g, 0]));
    levels.set(sym, 0);
    const lost = rx.filter((j) => !blocked.has(j) && gprActivity(model.rxns[j].g, levels) === 0).sort((a, b) => a - b);
    if (!lost.length) continue;
    const key = lost.join(",");
    if (!groups.has(key)) groups.set(key, { rxns: Int32Array.from(lost), genes: [] });
    groups.get(key).genes.push(sym);
  }
  const list = [...groups.values()];
  if (list.length) {
    const ref = (await sim.screenKnockouts(cfg, sup, [new Int32Array(0)]))[0].objective;
    const res = await sim.screenKnockouts(cfg, sup, list.map((g) => g.rxns));
    res.forEach((r, i) => { if (r.objective / ref < 0.01) for (const g of list[i].genes) ess.add(g); });
  }
  return [...ess].sort();
}

if (WORKER_ROLE === "media1" || WORKER_ROLE === "media2") {
  const meds = WORKER_ROLE === "media1" ? ["default", "rpmi1640", "dmem"] : ["plasmax", "hplm"];
  const out = {};
  for (const med of meds) {
    const t0 = Date.now();
    out[med] = await screenWithOff(new Set(), med, null);
    console.log(`${med}: ${out[med].length} essential genes in ${((Date.now() - t0) / 1000).toFixed(1)} s`);
  }
  fs.writeFileSync(TMP_DIR + `${WORKER_ROLE}.json`, JSON.stringify(out));
} else if (WORKER_ROLE === "audit50") {
  const tasks = parseTasks(fs.readFileSync(R + "data/tasks/metabolicTasks_Essential.txt", "utf8"));
  const tr = new TaskRunner(highs, model, tasks);
  const builder = new ContextBuilder(tr, model);
  const auditLines = candidates.slice(0, 50);
  const out = [];
  for (const thr of [0.5, 1.0, 2.0, 5.0]) {
    let sumMs = 0, sumOff0 = 0, sumRescued = 0, sumBlocked = 0, totalFailed = 0;
    for (const mid of auditLines) {
      const t0 = performance.now();
      const res = builder.build(getTpmMap(mid), thr, tr.tasks);
      sumMs += performance.now() - t0;
      sumOff0 += res.nOff0;
      sumRescued += res.rescued.length;
      sumBlocked += res.nBlocked;
      totalFailed += res.failed.length;
    }
    const rec = {
      thresholdTpm: thr, nLines: auditLines.length,
      avgBuildMs: +(sumMs / auditLines.length).toFixed(1),
      avgOff0: +(sumOff0 / auditLines.length).toFixed(1),
      avgRescued: +(sumRescued / auditLines.length).toFixed(1),
      avgBlocked: +(sumBlocked / auditLines.length).toFixed(1),
      totalFailedTasks: totalFailed,
    };
    console.log(`thr=${thr} TPM: avgBuild=${rec.avgBuildMs} ms, off0=${rec.avgOff0}, rescued=${rec.avgRescued}, failed=${totalFailed}`);
    out.push(rec);
  }
  fs.writeFileSync(TMP_DIR + "audit50.json", JSON.stringify(out));
} else if (WORKER_ROLE === "panel0" || WORKER_ROLE === "panel1") {
  const tasks = parseTasks(fs.readFileSync(R + "data/tasks/metabolicTasks_Essential.txt", "utf8"));
  const tr = new TaskRunner(highs, model, tasks);
  const builder = new ContextBuilder(tr, model);
  const slice = WORKER_ROLE === "panel0" ? panelLines.slice(0, 6) : panelLines.slice(6, 12);
  const out = { 0.5: {}, 1.0: {}, 2.0: {}, 5.0: {} };
  for (const mid of slice) {
    const tpm = getTpmMap(mid);
    for (const thr of [0.5, 1.0, 2.0, 5.0]) {
      const b = builder.build(tpm, thr, tr.tasks);
      out[thr][mid] = await screenWithOff(b.off, "defaultSerum", genericDefSerum);
    }
    console.log(`screened ${mid} (${dep.lines[mid].name}) across 4 TPM cutoffs`);
  }
  fs.writeFileSync(TMP_DIR + `${WORKER_ROLE}.json`, JSON.stringify(out));
} else if (WORKER_ROLE === "panelMatched") {
  const tasks = parseTasks(fs.readFileSync(R + "data/tasks/metabolicTasks_Essential.txt", "utf8"));
  const tr = new TaskRunner(highs, model, tasks);
  const builder = new ContextBuilder(tr, model);
  const m1 = JSON.parse(fs.readFileSync(TMP_DIR + "media1.json", "utf8"));
  const genByMed = { rpmi1640: new Set(m1.rpmi1640), dmem: new Set(m1.dmem), defaultSerum: genericDefSerum };
  const out = {};
  for (const mid of panelLines) {
    const med = mediaByLine.get(mid) || "defaultSerum";
    const b = builder.build(getTpmMap(mid), 1.0, tr.tasks);
    out[mid] = await screenWithOff(b.off, med, genByMed[med]);
    console.log(`screened ${mid} on matched ${med}: ${out[mid].length} ess`);
  }
  fs.writeFileSync(TMP_DIR + "panelMatched.json", JSON.stringify(out));
}
