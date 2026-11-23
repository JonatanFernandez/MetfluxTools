// Cell-type contexts: which genes are switched off in each cell type, after rescuing the metabolic tasks.
//
//   node tools/contexts/build_contexts.mjs <slug> [shard/nShards]      -> tools/cache/contexts/<slug>_part<k>.json
//   node tools/contexts/build_contexts.mjs <slug> merge <nShards>       -> data/<slug>/contexts.json
//
// For each context (tools/cache/contexts/<slug>_expr.json, made by build_expression.py):
//  1. A gene is off when it is not detected: below 1 nTPM / nCPM in the Human Protein Atlas (HPA's own
//     detection cut-off), or called "absent" by Bgee (mouse). Genes without data stay on.
//  2. A reaction is blocked when its GPR rule evaluates to 0 with those genes off (and = min, or = max).
//     Transport of freely diffusing molecules (H2O, O2, CO2, NH3, urea, NO, H2S) is never blocked.
//  3. Every essential metabolic task of Human-GEM (data/metabolicTasks/metabolicTasks_Essential.txt, 57 tasks,
//     RAVEN semantics with closed boundaries) must still be possible. For a task that is not, an LP finds the
//     cheapest set of blocked reactions that makes it possible (minimise sum w_j |v_j| over blocked reactions,
//     w_j = 2 - level of the best-expressed gene, so genes just under the cut-off are preferred) and the genes
//     those reactions need are switched back on (for "or", the branch needing fewest/most expressed genes).
//     This is the task-protection idea of tINIT (Agren et al. 2014) applied to a gene-level knockout list.
//  4. Contexts paired with a measured exchange-flux set must also reproduce that phenotype (phenotype task).

const R = new URL("../../", import.meta.url).pathname;
const fs = await import("node:fs");
const { execSync } = await import("node:child_process");
const { default: loadHighs } = await import(R + "vendor/highs/highs.mjs");
const { parseTasks, TaskRunner } = await import(R + "tools/contexts/tasks.mjs");

const slug = process.argv[2] || "human-gem";
const OUT = R + "tools/cache/contexts/";
export const THRESHOLD = { "human-gem": 1, "mouse-gem": 1e-9 }; // HPA: detected >= 1 nTPM; Bgee: absent = 0

// measured sets that describe the same cells as a context ("phenotype tasks")
const PAIRS = {
  "human-gem": {
    nci60: [{ set: "nci60", medium: "default" }],
    "hpa_line:HEK293": [{ set: "hek293", medium: "default" }],
    "hpa_line:BJ [Human fibroblast]": [{ set: "fibroP", medium: "default" }, { set: "fibroQ", medium: "default" }],
    "hpa_tissue:liver": [{ set: "liver", medium: "williamsE" }],
    "hpa_cell:hepatocytes": [{ set: "liver", medium: "williamsE" }],
  },
  "mouse-gem": {},
};

// ---- build ----------------------------------------------------------------------------------------------
async function build(shard, nShards) {
  const highs = await loadHighs({ print: () => {}, printErr: () => {} });
  const model = JSON.parse(fs.readFileSync(R + `data/${slug}/model.json`));
  const expr = JSON.parse(fs.readFileSync(OUT + `${slug}_expr.json`));
  const tasks = parseTasks(execSync(`unzip -p ${R}Human-GEM-2.1.1.zip Human-GEM-2.1.1/data/metabolicTasks/metabolicTasks_Essential.txt`).toString());
  const tr = new TaskRunner(highs, model, tasks);
  const nR = model.rxns.length;
  const { ContextBuilder, active } = await import(R + "js/core/ctxbuild.js");
  const cb = new ContextBuilder(tr, model), asts = cb.asts;
  const thr = THRESHOLD[slug] ?? 1e-9; // every other model: Bgee calls

  // phenotype tasks: bounds as the app sets them (medium + measured set), captured from the real Simulator
  const pairs = PAIRS[slug] ?? {};
  const phenoCache = new Map();
  async function phenotype(p) {
    const key = p.set + "|" + p.medium;
    if (phenoCache.has(key)) return phenoCache.get(key);
    const { Simulator } = await import(R + "js/sim.js");
    const { presetMedium } = await import(R + "js/core/medium.js");
    const { measuredPreset, measuredNgam } = await import(R + "js/core/measured.js");
    const { ModelIndex } = await import(R + "js/modelindex.js");
    const structures = JSON.parse(fs.readFileSync(R + `data/${slug}/structures.json`));
    const atommaps = JSON.parse(fs.readFileSync(R + `data/${slug}/atommaps.json`));
    const layout = JSON.parse(fs.readFileSync(R + `data/${slug}/layout.json`));
    let captured = null;
    const backend = async (msg) => {
      if (msg.type === "init") return true;
      if (!captured) captured = { lb: Float64Array.from(msg.params.lb), ub: Float64Array.from(msg.params.ub) };
      return { status: "optimal", objective: 0, v: new Float64Array(nR), sourceFlux: new Float64Array(msg.params.sources.length) };
    };
    const sim = new Simulator(model, structures, atommaps, backend), index = new ModelIndex(model, layout);
    const meas = { items: measuredPreset(model, sim.exchanges, p.set).map((x) => ({ met: x.met, rate: x.rate })), tol: 0.2, ngam: measuredNgam(p.set) };
    await sim.runSteady({ objective: index.biomassRxn, maximize: true, parsimonious: false, measured: meas }, presetMedium(model, p.medium), structures);
    const t = { id: "PH", desc: `Measured phenotype ${p.set} (${p.medium} medium)`, bounds: tr.modelBounds(captured.lb, captured.ub) };
    phenoCache.set(key, t);
    return t;
  }

  // pool of known feasible flux supports per task, seeded with a sparse (min total flux) full-model solution

  // how often each reaction is blocked over all contexts: a small cost that steers solutions away from them
  const freqW = cb.freqW; freqW.fill(0);
  for (const c of expr.contexts) {
    const off = new Set(expr.genes.filter((g, i) => c.levels[i] !== null && c.levels[i] < thr));
    for (let j = 0; j < nR; j++) if (asts[j] && !active(asts[j], off)) freqW[j] += 1;
  }
  for (let j = 0; j < nR; j++) freqW[j] = 1e-3 + freqW[j] / expr.contexts.length;
  console.log("block frequencies done");
  let all = expr.contexts.filter((_, i) => i % nShards === shard);
  if (process.env.ONLY) { const ids = new Set(process.env.ONLY.split(",")); all = expr.contexts.filter((c) => ids.has(c.id)); }
  // resume: contexts already written to this shard's part file are kept and skipped
  const partFile = OUT + `${slug}_part${shard}.json`;
  const out = !process.env.ONLY && fs.existsSync(partFile) ? JSON.parse(fs.readFileSync(partFile)) : [];
  const done = new Set(out.map((c) => c.id));
  if (done.size) { all = all.filter((c) => !done.has(c.id)); console.log(`resuming: ${done.size} done, ${all.length} to go`); }
  const t0 = Date.now();
  for (const [ci, c] of all.entries()) {
    const level = new Map();
    expr.genes.forEach((g, i) => { const x = c.levels[i]; if (x !== null && x !== undefined) level.set(g, x); });
    const ctxTasks = [...tr.tasks.filter((t) => !t.shouldFail && !t.bad)];
    for (const p of pairs[c.id] || []) ctxTasks.push(await phenotype(p));
    const { off, nOff0, rescued, failed, nBlocked } = cb.build(level, thr, ctxTasks);
    const geneIdx = new Map(model.genes.map((g, i) => [g[1], i]));
    out.push({ id: c.id, name: c.name, group: c.group, src: c.src, info: c.info,
      off: [...off].map((g) => geneIdx.get(g)).sort((a, b) => a - b), nOff0, rescued: rescued.map(([g, t]) => [geneIdx.get(g), t]),
      nBlocked, failed, pairs: (pairs[c.id] || []).map((p) => p.set) });
    if (ci % 10 === 0 || ci === all.length - 1) {
      const el = (Date.now() - t0) / 1000;
      console.log(`[${shard}/${nShards}] ${ci + 1}/${all.length} ${c.id}: off ${nOff0} -> ${off.size}, blocked ${out.at(-1).nBlocked}, rescued ${rescued.length}, failed ${failed.length} | ${el.toFixed(0)} s, LPs ${cb.nLp} + ${cb.nRescueLp} rescue`);
      fs.writeFileSync(OUT + `${slug}_part${shard}.json`, JSON.stringify(out));
    }
  }
  fs.writeFileSync(OUT + `${slug}_part${shard}.json`, JSON.stringify(out));
}

// ---- merge into the file the web app loads --------------------------------------------------------------
const GENERIC_ENTITIES = new Set(["UBERON:0000465", "UBERON:0000468", "UBERON:0001062", "UBERON:0000061", "UBERON:0010000"]);
function merge(nShards) {
  const model = JSON.parse(fs.readFileSync(R + `data/${slug}/model.json`));
  const expr = JSON.parse(fs.readFileSync(OUT + `${slug}_expr.json`));
  const order = new Map(expr.contexts.map((c, i) => [c.id, i]));
  let items = [];
  for (let k = 0; k < nShards; k++) items.push(...JSON.parse(fs.readFileSync(OUT + `${slug}_part${k}.json`)));
  items.sort((a, b) => order.get(a.id) - order.get(b.id));
  if (items.length !== expr.contexts.length) throw new Error(`incomplete: ${items.length} of ${expr.contexts.length}`);
  // Bgee root terms ("material anatomical entity", "multicellular organism"...) pool every sample: not a cell type
  items = items.filter((c) => !GENERIC_ENTITIES.has(c.id.replace(/^bgee:/, "")));
  const nG = model.genes.length;
  const bits = (idx) => { const b = new Uint8Array(Math.ceil(nG / 8)); for (const i of idx) b[i >> 3] |= 1 << (i & 7); return Buffer.from(b).toString("base64"); };
  const human = slug === "human-gem";
  // task names are stored once; rescued genes and failed tasks refer to them by index
  const taskNames = [], taskIdx = new Map();
  const ti = (name) => { if (!taskIdx.has(name)) { taskIdx.set(name, taskNames.length); taskNames.push(name); } return taskIdx.get(name); };
  const doc = {
    meta: human ? {
      sources: ["hpaLine", "hpaCell", "hpaTissue", "cellosaurus"],
      hpaVersion: "26.0", cellosaurusVersion: "57.0", threshold: "1 nTPM (cell lines, tissues) / 1 nCPM (single cell types)",
      license: "CC BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0/). Derived from the Human Protein Atlas (proteinatlas.org, version 26.0, CC BY-SA 4.0): expression levels converted to lists of genes not detected, with metabolic-task rescue. Cell line descriptions from Cellosaurus 57.0 (CC BY 4.0).",
    } : {
      sources: ["bgee"], threshold: "Bgee absent call",
      license: "CC0 1.0. Derived from Bgee (bgee.org, CC0): absent calls converted to lists of genes switched off, with metabolic-task rescue.",
    },
    tasks: "Human-GEM 2.1.1 data/metabolicTasks/metabolicTasks_Essential.txt (57 tasks)",
    genes: nG,
    taskNames,
    contexts: items.map((c) => ({ id: c.id, name: c.name, group: c.group, src: c.src, info: c.info, off: bits(c.off), nOff0: c.nOff0,
      rescued: c.rescued.map(([g, t]) => [g, ti(t)]), nBlocked: c.nBlocked, failed: c.failed.map(ti), pairs: c.pairs })),
  };
  fs.writeFileSync(R + `data/${slug}/contexts.json`, JSON.stringify(doc));
  console.log("wrote", items.length, "contexts,", (fs.statSync(R + `data/${slug}/contexts.json`).size / 1e6).toFixed(2), "MB");
}

if (process.argv[3] === "merge") merge(+process.argv[4] || 1);
else { const [k, n] = (process.argv[3] || "0/1").split("/").map(Number); await build(k, n); }
