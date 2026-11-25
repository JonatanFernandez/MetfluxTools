// Build compact signed TF-target regulon files (data/<slug>/regulons.json) for MetfluxLive:
//   - human-gem: CollecTRI via OmniPath API (organism 9606)
//   - mouse-gem: CollecTRI via OmniPath API (organism 10090)
//   - ecoli-gem: RegulonDB NetworkRegulatorGene via GraphQL API
//
// Format of data/<slug>/regulons.json:
//   {
//     meta: { source, version, organism, built, nTFs, nMetGenes, nGenes, nMetEdges, nAllEdges },
//     nMet: <number>,           // indices 0..nMet-1 are metabolic genes present in model.json
//     genes: [sym0, sym1, ...], // canonical target gene symbols (matching model.json case)
//     tfs: [ [tfName, [geneIdx, sign, ...]], ... ]   // sign is +1 (activates) or -1 (represses)
//   }
//
// Usage: node tools/build_regulons.mjs
import fs from "node:fs";
import path from "node:path";

const ROOT = new URL("../", import.meta.url).pathname;
const CACHE = path.join(ROOT, "tools/cache/regulons");
fs.mkdirSync(CACHE, { recursive: true });

async function cachedFetch(cacheName, fetcher) {
  const fp = path.join(CACHE, cacheName);
  if (fs.existsSync(fp) && fs.statSync(fp).size > 1000) {
    console.log(`  [cache] ${cacheName}`);
    return fs.readFileSync(fp, "utf8");
  }
  console.log(`  [fetch] ${cacheName} ...`);
  const txt = await fetcher();
  fs.writeFileSync(fp, txt, "utf8");
  return txt;
}

/** Map of lowercase gene symbol (and id) -> canonical symbol in model.json */
function modelGeneMap(slug) {
  const m = JSON.parse(fs.readFileSync(path.join(ROOT, `data/${slug}/model.json`), "utf8"));
  const map = new Map();
  for (const g of m.genes) {
    const id = g[0], sym = g[1] || g[0];
    if (sym) map.set(sym.toLowerCase(), sym);
    if (id && !map.has(id.toLowerCase())) map.set(id.toLowerCase(), sym);
  }
  return map;
}

/**
 * Compile raw (tf, target, sign) edges into the compact indexed JSON format.
 * Multiple edges for the same (tf, target) are deduplicated; conflicting signs cancel out and are dropped.
 */
function compileRegulons(rawEdges, slug, metaExtra, { minTargets = 3 } = {}) {
  const metMap = modelGeneMap(slug);
  // Deduplicate (tf, targetLower) -> { tf, targetSym, isMet, pos, neg }
  const tfMap = new Map();
  for (const { tf, target, sign } of rawEdges) {
    if (!tf || !target || (sign !== 1 && sign !== -1)) continue;
    const tLow = target.trim().toLowerCase();
    const isMet = metMap.has(tLow);
    const canon = isMet ? metMap.get(tLow) : target.trim();
    let edges = tfMap.get(tf);
    if (!edges) tfMap.set(tf, (edges = new Map()));
    const prev = edges.get(tLow) || { canon, isMet, pos: 0, neg: 0 };
    if (sign > 0) prev.pos++; else prev.neg++;
    edges.set(tLow, prev);
  }

  // Collect clean signed edges per TF
  const cleanTFs = [];
  const metSyms = new Set(), otherSyms = new Set();
  for (const [tf, edges] of tfMap) {
    const kept = [];
    for (const [, e] of edges) {
      if (e.pos && e.neg) continue; // ambiguous / dual regulation dropped
      const sign = e.pos ? 1 : -1;
      kept.push([e.canon, sign, e.isMet]);
    }
    if (kept.length < minTargets) continue;
    for (const [canon, , isMet] of kept) (isMet ? metSyms : otherSyms).add(canon);
    cleanTFs.push([tf, kept]);
  }

  // Order genes: metabolic genes first (0..nMet-1), then non-metabolic genome-wide targets
  const metList = [...metSyms].sort();
  const otherList = [...otherSyms].sort();
  const genes = [...metList, ...otherList];
  const nMet = metList.length;
  const gIdx = new Map(genes.map((g, i) => [g, i]));

  let nMetEdges = 0, nAllEdges = 0;
  const tfs = cleanTFs
    .sort((a, b) => a[0].localeCompare(b[0]))
    .map(([tf, kept]) => {
      kept.sort((a, b) => gIdx.get(a[0]) - gIdx.get(b[0]));
      const flat = [];
      for (const [canon, sign, isMet] of kept) {
        flat.push(gIdx.get(canon), sign);
        nAllEdges++;
        if (isMet) nMetEdges++;
      }
      return [tf, flat];
    });

  const out = {
    meta: {
      ...metaExtra,
      built: new Date().toISOString().slice(0, 10),
      nTFs: tfs.length,
      nMetGenes: nMet,
      nGenes: genes.length,
      nMetEdges,
      nAllEdges,
    },
    nMet,
    genes,
    tfs,
  };
  const outPath = path.join(ROOT, `data/${slug}/regulons.json`);
  const json = JSON.stringify(out);
  fs.writeFileSync(outPath, json, "utf8");
  console.log(
    `  -> ${slug}/regulons.json: ${tfs.length} TFs, ${nMet}/${genes.length} metabolic/all target genes, `
    + `${nMetEdges}/${nAllEdges} metabolic/all signed edges (${(json.length / 1024).toFixed(1)} KB)`
  );
  return out;
}

/** Parse OmniPath CollecTRI TSV into signed edges. */
function parseCollecTRI(tsv) {
  const lines = tsv.split(/\r?\n/).filter((l) => l.trim());
  const head = lines[0].split("\t");
  const col = (name) => head.indexOf(name);
  const iSrc = col("source_genesymbol"), iTgt = col("target_genesymbol");
  const iCS = col("consensus_stimulation"), iCI = col("consensus_inhibition");
  const iS = col("is_stimulation"), iI = col("is_inhibition");
  const edges = [];
  for (let i = 1; i < lines.length; i++) {
    const c = lines[i].split("\t");
    const tf = c[iSrc]?.trim(), target = c[iTgt]?.trim();
    if (!tf || !target) continue;
    const cs = c[iCS] === "True", ci = c[iCI] === "True";
    const s = c[iS] === "True", inh = c[iI] === "True";
    let sign = 0;
    if (cs && !ci) sign = 1;
    else if (ci && !cs) sign = -1;
    else if (!cs && !ci) {
      if (s && !inh) sign = 1;
      else if (inh && !s) sign = -1;
    }
    if (sign) edges.push({ tf, target, sign });
  }
  return edges;
}

/** Parse RegulonDB NetworkRegulatorGene TSV into signed edges. */
function parseRegulonDB(tsv) {
  const edges = [];
  for (const line of tsv.split(/\r?\n/)) {
    if (!line || line.startsWith("#") || line.startsWith("1)")) continue;
    const c = line.split("\t");
    if (c.length < 6) continue;
    const tf = (c[1] || "").trim();
    const target = (c[4] || "").trim();
    const fn = (c[5] || "").trim();
    if (!tf || !target) continue;
    if (fn === "+") edges.push({ tf, target, sign: 1 });
    else if (fn === "-") edges.push({ tf, target, sign: -1 });
  }
  return edges;
}

async function main() {
  // 1. Human CollecTRI
  console.log("Building human-gem regulons (CollecTRI / OmniPath)...");
  const humanTsv = await cachedFetch("collectri_9606.tsv", async () => {
    const r = await fetch("https://omnipathdb.org/interactions?datasets=collectri&organisms=9606&genesymbols=yes&fields=sources");
    if (!r.ok) throw new Error(`OmniPath human HTTP ${r.status}`);
    return r.text();
  });
  compileRegulons(parseCollecTRI(humanTsv), "human-gem", {
    source: "CollecTRI (OmniPath)",
    organism: "Homo sapiens (9606)",
    license: "Open academic/research via OmniPath (Müller-Dott et al., NAR 2023)",
  });

  // 2. Mouse CollecTRI
  console.log("Building mouse-gem regulons (CollecTRI / OmniPath)...");
  const mouseTsv = await cachedFetch("collectri_10090.tsv", async () => {
    const r = await fetch("https://omnipathdb.org/interactions?datasets=collectri&organisms=10090&genesymbols=yes&fields=sources");
    if (!r.ok) throw new Error(`OmniPath mouse HTTP ${r.status}`);
    return r.text();
  });
  compileRegulons(parseCollecTRI(mouseTsv), "mouse-gem", {
    source: "CollecTRI (OmniPath orthology)",
    organism: "Mus musculus (10090)",
    license: "Open academic/research via OmniPath (Müller-Dott et al., NAR 2023)",
  });

  // 3. E. coli RegulonDB
  console.log("Building ecoli-gem regulons (RegulonDB NetworkRegulatorGene)...");
  let rdbVersion = "14.5.1";
  const ecoliTsv = await cachedFetch("regulondb_NetworkRegulatorGene.tsv", async () => {
    const r = await fetch("https://regulondb.ccg.unam.mx/graphql", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ query: '{ getDataOfFile(fileName: "NetworkRegulatorGene") { fileName rdbVersion content } }' }),
    });
    if (!r.ok) throw new Error(`RegulonDB GraphQL HTTP ${r.status}`);
    const j = await r.json();
    const f = j.data?.getDataOfFile;
    if (!f?.content) throw new Error("RegulonDB returned empty content");
    rdbVersion = f.rdbVersion || rdbVersion;
    return f.content;
  });
  compileRegulons(parseRegulonDB(ecoliTsv), "ecoli-gem", {
    source: `RegulonDB v${rdbVersion} (NetworkRegulatorGene)`,
    organism: "Escherichia coli K-12 MG1655",
    license: "Free for academic/noncommercial use (CCG-UNAM)",
  });
}

main().catch((err) => { console.error(err); process.exit(1); });
