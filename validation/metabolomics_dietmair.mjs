// Validation 6: intracellular metabolomics. The producer / non-producer contrast of HEK293
// (Dietmair et al., PLoS ONE 7:e43394, 2012, Table S3; tools/validation/prepare_dietmair.py) placed on the
// Human-GEM network, with the exchange rates of the same experiment (the "hek293" measured set, their
// Table 2) as the flux state.
//
// This is the test case for the intracellular side: a real two-condition contrast whose cause is known from
// the paper (the producer consumes less glucose and less glutamine at the same growth rate, with oxidative
// phosphorylation genes down and ER-stress genes up), measured on 52 metabolites of central metabolism.
//
// What it reports now:
//   1. mapping    how many measured metabolites exist in the model, and which are ambiguous isomer pairs
//   2. coverage   the "observability" of the panel: reactions with both a substrate and a product measured,
//                 how many measured metabolites have another one within one or two steps, pathways covered
//   3. contrast   log2 of producer / non-producer per metabolite, with Welch's t over the 4 time points,
//                 and the pathways where the panel moves as a whole (Mann-Whitney against the rest)
//
//   node tools/validation/metabolomics_dietmair.mjs [--json out.json]
import fs from "node:fs";

const R = new URL("../../", import.meta.url).pathname;
const model = JSON.parse(fs.readFileSync(R + "data/human-gem/model.json"));
const layout = JSON.parse(fs.readFileSync(R + "data/human-gem/layout.json"));
const data = JSON.parse(fs.readFileSync(R + "tools/validation/dietmair2012_intracellular.json"));

const B = (i) => model.mets[i].b;
const baseByName = new Map(model.bases.map((b, i) => [b.n.toLowerCase(), i]));

// ---------------------------------------------------------------- 1. mapping ----
const rows = [];
for (const m of data.metabolites) {
  const names = (m.metabolite || "").split("|").filter(Boolean);
  const bases = names.map((n) => baseByName.get(n.toLowerCase())).filter((x) => x !== undefined);
  rows.push({ ...m, names, bases, ambiguous: names.length > 1 });
}
const mapped = rows.filter((r) => r.bases.length);
const unmapped = rows.filter((r) => !r.bases.length);
console.log(`== mapping\n${mapped.length} of ${rows.length} measured metabolites are in Human-GEM`
  + (unmapped.length ? `; not found: ${unmapped.map((r) => `${r.label} (${r.metabolite})`).join(", ")}` : ""));
const amb = rows.filter((r) => r.ambiguous);
if (amb.length) console.log(`   ambiguous measurements kept as one node: ${amb.map((r) => `${r.label} = ${r.names.join(" or ")}`).join("; ")}`);

// ---------------------------------------------------------------- 2. coverage ----
// metabolite-metabolite graph over metabolic reactions, without the currency metabolites and the biggest
// hubs (through those, everything is two steps from everything)
const currency = new Set(layout.currency.map(B));
const degree = new Map();
for (const r of model.rxns) { if (r.removed) continue; for (const i of new Set(r.m.map(B))) degree.set(i, (degree.get(i) || 0) + 1); }
for (const [b] of [...degree].sort((a, b) => b[1] - a[1]).slice(0, 40)) currency.add(b);
const reactions = [];
model.rxns.forEach((r) => {
  if (r.removed || r.m.length < 2) return;
  const ss = model.subsystems[r.ss] || "";
  if (/transport|exchange|pool|biomass|artificial/i.test(ss) || /pool|biomass/i.test(r.n)) return;
  const sub = new Set(), prod = new Set();
  r.m.forEach((i, k) => { const b = B(i); if (currency.has(b)) return; (r.c[k] < 0 ? sub : prod).add(b); });
  for (const b of [...sub]) if (prod.has(b)) { sub.delete(b); prod.delete(b); }
  if (sub.size && prod.size) reactions.push({ sub: [...sub], prod: [...prod], ss });
});
const adj = new Map();
const link = (a, b) => { if (a !== b) (adj.get(a) || adj.set(a, new Set()).get(a)).add(b); };
for (const r of reactions) for (const a of r.sub) for (const b of r.prod) { link(a, b); link(b, a); }

const measured = new Set(mapped.flatMap((r) => r.bases));
const bothSides = reactions.filter((r) => r.sub.some((b) => measured.has(b)) && r.prod.some((b) => measured.has(b)));
const hops = (src) => {
  const seen = new Set([src]); let front = [src];
  for (let k = 1; k <= 3; k++) {
    const next = [];
    for (const a of front) for (const b of adj.get(a) || []) { if (seen.has(b)) continue; if (measured.has(b)) return k; seen.add(b); next.push(b); }
    front = next;
  }
  return Infinity;
};
const d = [...measured].map(hops);
const pct = (n) => `${Math.round((100 * n) / d.length)} %`;
const subsystems = new Map();
for (const r of bothSides) subsystems.set(r.ss, (subsystems.get(r.ss) || 0) + 1);
console.log(`\n== coverage of the network\nreactions with substrate AND product measured: ${bothSides.length}`);
console.log(`another measured metabolite 1 step away: ${pct(d.filter((x) => x <= 1).length)}; within 2 steps: ${pct(d.filter((x) => x <= 2).length)}`);
console.log(`pathways reached: ${subsystems.size}. The best covered:`);
for (const [ss, n] of [...subsystems].sort((a, b) => b[1] - a[1]).slice(0, 8)) console.log(`   ${String(n).padStart(3)}  ${ss}`);

// ---------------------------------------------------------------- 3. the contrast ----
// relative concentrations: only producer / non-producer means anything. Welch's t over the 4 time points,
// which are not replicates of one state but a course — so this is a difference between the two lines over
// the culture, not a per-time-point test.
const mean = (a) => a.reduce((s, x) => s + x, 0) / a.length;
const vari = (a) => (a.length < 2 ? 0 : a.reduce((s, x) => s + (x - mean(a)) ** 2, 0) / (a.length - 1));
const erf = (x) => { // Abramowitz & Stegun 7.1.26
  const t = 1 / (1 + 0.3275911 * Math.abs(x));
  const y = 1 - ((((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
  return x >= 0 ? y : -y;
};
const normP = (z) => 1 - erf(Math.abs(z) / Math.SQRT2); // two-sided, normal approximation
const contrast = mapped.map((r) => {
  const p = r.producer.map((x) => x.mean), q = r.nonProducer.map((x) => x.mean);
  const mp = mean(p), mq = mean(q);
  const se = Math.sqrt(vari(p) / p.length + vari(q) / q.length);
  return { label: r.label, bases: r.bases, log2fc: Math.log2(mp / mq), t: se > 0 ? (mp - mq) / se : 0, p: se > 0 ? normP((mp - mq) / se) : 1 };
});
// Benjamini-Hochberg
const byP = [...contrast].sort((a, b) => a.p - b.p);
byP.forEach((x, i) => { x.q = Math.min(1, (x.p * byP.length) / (i + 1)); });
for (let i = byP.length - 2; i >= 0; i--) byP[i].q = Math.min(byP[i].q, byP[i + 1].q);
const sig = contrast.filter((x) => x.q < 0.05).sort((a, b) => b.log2fc - a.log2fc);
console.log(`\n== producer against non-producer\nmedian log2 ratio over the whole panel: ${mean(contrast.map((x) => x.log2fc)).toFixed(3)} (a whole-panel shift would show here, and would be a normalisation artefact rather than biology)`);
console.log(`${sig.length} of ${contrast.length} metabolites differ (Benjamini-Hochberg q < 0.05):`);
for (const x of sig) console.log(`   ${x.label.padEnd(12)} ${x.log2fc > 0 ? "+" : ""}${x.log2fc.toFixed(2)}  q = ${x.q.toExponential(1)}`);

// Classic ratios. With RELATIVE data the value of a ratio between two different metabolites is
// meaningless — each one has its own response factor — but the CHANGE of that ratio between the two lines
// is valid, because the factors cancel: Δlog2(A/B) = log2FC(A) − log2FC(B).
const fc = new Map(contrast.map((x) => [x.label, x.log2fc]));
const RATIOS = [["ATP", "AMP", "energy charge: ATP against AMP"], ["ATP", "ADP", "ATP against ADP"],
  ["LAC", "PYR", "lactate / pyruvate, the cytosolic NADH/NAD+ state"], ["GLU", "KGA", "glutamate / 2-oxoglutarate"],
  ["CIT", "MAL", "citrate / malate"]];
console.log(`\n== ratios, as the change between the lines (the only form a relative measurement supports)`);
for (const [a, b, what] of RATIOS) {
  if (!fc.has(a) || !fc.has(b)) continue;
  const d2 = fc.get(a) - fc.get(b);
  console.log(`   ${(a + "/" + b).padEnd(9)} ${d2 > 0 ? "+" : ""}${d2.toFixed(2)} log2 (x${(2 ** d2).toFixed(2)}) in the producer — ${what}`);
}

// pathway-level shift: Mann-Whitney of each subsystem's metabolites against the rest of the panel
const inSub = new Map();
for (const r of reactions) for (const b of [...r.sub, ...r.prod]) {
  if (!measured.has(b)) continue;
  (inSub.get(r.ss) || inSub.set(r.ss, new Set()).get(r.ss)).add(b);
}
const byBase = new Map();
for (const x of contrast) for (const b of x.bases) byBase.set(b, x.log2fc);
const U = (a, b) => { // Mann-Whitney U with a normal approximation
  const all = [...a.map((v) => [v, 0]), ...b.map((v) => [v, 1])].sort((p, q) => p[0] - q[0]);
  let rank = 0, ra = 0;
  for (let i = 0; i < all.length; i++) { rank = i + 1; if (all[i][1] === 0) ra += rank; }
  const n1 = a.length, n2 = b.length;
  const u = ra - (n1 * (n1 + 1)) / 2;
  const mu = (n1 * n2) / 2, sd = Math.sqrt((n1 * n2 * (n1 + n2 + 1)) / 12);
  return { u, p: sd > 0 ? normP((u - mu) / sd) : 1 };
};
const shifts = [];
for (const [ss, set] of inSub) {
  const inside = [...set].map((b) => byBase.get(b)).filter(Number.isFinite);
  if (inside.length < 4) continue;
  const outside = contrast.filter((x) => !x.bases.some((b) => set.has(b))).map((x) => x.log2fc);
  if (outside.length < 4) continue;
  const { p } = U(inside, outside);
  shifts.push({ ss, n: inside.length, median: [...inside].sort((a, b) => a - b)[inside.length >> 1], p });
}
shifts.sort((a, b) => a.p - b.p);
console.log(`\n== pathways that move as a whole (Mann-Whitney against the rest of the panel, uncorrected)`);
for (const s of shifts.slice(0, 6)) console.log(`   p = ${s.p.toFixed(4)}  median log2 ${s.median > 0 ? "+" : ""}${s.median.toFixed(2)}  (${s.n} measured)  ${s.ss}`);
console.log("\nThe statistics of different metabolites are not independent (shared pathways, shared normalisation),\nso these p-values are optimistic: with the sample matrix the right null is a permutation of the sample labels.");

const out = process.argv.indexOf("--json");
if (out > 0 && process.argv[out + 1]) {
  fs.writeFileSync(process.argv[out + 1], JSON.stringify({ mapped: mapped.length, total: rows.length,
    coverage: { reactions: bothSides.length, within1: d.filter((x) => x <= 1).length / d.length, within2: d.filter((x) => x <= 2).length / d.length, pathways: subsystems.size },
    contrast, shifts }, null, 1));
  console.log(`\nwritten to ${process.argv[out + 1]}`);
}
