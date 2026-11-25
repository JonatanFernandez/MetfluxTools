// The worked-example library (js/ui/demos.js) must stay honest: every step needs its texts and an element
// to point at, every action has to be one the application can replay, and the data a demo ships has to
// reproduce the published numbers it claims to come from.
//
// This test exists because the first version of the intracellular table was transcribed by hand and 7 of
// its 9 metabolites were wrong. A worked example with wrong data teaches the wrong thing, and nothing else
// in the suite would have caught it.
//   node tools/test_demos.mjs
import fs from "node:fs";
import { DEMOS, ACTIONS, EXO_TABLE, MTB_TABLE, RNA_TABLE, REG_TABLE } from "../js/ui/demos.js";
import { TOUR_STEPS, TOUR_STEPS_MOBILE } from "../js/ui/tour.js";
import { DICT } from "../js/i18n.js";
import { parseExoTable, groupForRates, ratesFromReplicates, growthRate } from "../js/core/exodata.js";
import { parseWideTable, splitByPattern, contrast, ratioChanges } from "../js/core/metcontrast.js";
import { parseRnaDiffTable, inferRegulatorsMultiOmic } from "../js/core/reginfer.js";

const R = new URL("../", import.meta.url).pathname;
let bad = 0;
const ok = (c, msg) => { console.log((c ? "ok  " : "FAIL") + " " + msg); if (!c) bad++; };
const near = (a, b, tol, msg) => ok(Math.abs(a - b) <= tol, `${msg} (${a.toPrecision(4)} vs ${b.toPrecision(4)}, tol ${tol})`);

const html = fs.readFileSync(R + "index.html", "utf8");
const mobileJs = fs.readFileSync(R + "js/ui/mobile.js", "utf8");
const docsHtml = fs.readFileSync(R + "docs.html", "utf8");
const extraLangs = ["de", "el", "fr", "it", "pt"];
const extraDicts = Object.fromEntries(await Promise.all(
  extraLangs.map(async (l) => [l, (await import(`../js/i18n/${l}.js`)).default])
));
const extraDocs = Object.fromEntries(
  extraLangs.map((l) => [l, fs.readFileSync(`${R}docs/${l}.html`, "utf8")])
);
const allDicts = { en: DICT.en, es: DICT.es, ...extraDicts };

function hasSelectorId(sel, allowMobile = false) {
  if (!sel) return true;
  // Check every comma-separated alternative; at least one leading #id must exist in index.html (or mobile.js).
  const alts = sel.split(",").map((s) => s.trim());
  return alts.some((a) => {
    const m = a.match(/^#([^\s.>+~:[\]]+)/);
    if (!m) return true;
    const id = m[1];
    if (html.includes(`id="${id}"`)) return true;
    if (allowMobile && (mobileJs.includes(`id="${id}"`) || mobileJs.includes(`id = "${id}"`) || mobileJs.includes(`"${id}"`))) return true;
    return false;
  });
}

console.log("== general walkthrough (desktop & mobile) is complete in all 7 languages");
{
  const missDeskSel = TOUR_STEPS.map((s) => s.sel).filter((s) => !hasSelectorId(s, false));
  ok(!missDeskSel.length, `desktop walkthrough: every selector exists in index.html${missDeskSel.length ? " — missing " + missDeskSel.join(", ") : ""}`);
  const missMobSel = TOUR_STEPS_MOBILE.map((s) => s.sel).filter((s) => !hasSelectorId(s, true));
  ok(!missMobSel.length, `mobile walkthrough: every selector exists${missMobSel.length ? " — missing " + missMobSel.join(", ") : ""}`);
  const tourKeys = [...new Set([...TOUR_STEPS, ...TOUR_STEPS_MOBILE].flatMap((s) => [`tour.${s.id}.t`, `tour.${s.id}.b`]))];
  for (const [lang, dict] of Object.entries(allDicts)) {
    const miss = tourKeys.filter((k) => !(k in dict));
    ok(!miss.length, `walkthrough texts in ${lang}${miss.length ? " — missing " + miss.join(", ") : ""}`);
  }
}

console.log("\n== every worked example is complete across all 7 languages and linked to Docs");
for (const d of DEMOS) {
  const keys = [`demo.${d.id}.title`, `demo.${d.id}.what`]
    .concat(d.steps.flatMap((s) => [`tour.${s.id}.t`, `tour.${s.id}.b`]));
  const missByLang = Object.entries(allDicts)
    .map(([l, dict]) => [l, keys.filter((k) => !(k in dict))])
    .filter(([, m]) => m.length);
  ok(!missByLang.length, `${d.id}: texts in all 7 languages${missByLang.length ? " — missing " + missByLang.map(([l, m]) => `${l}:${m.join(",")}`).join("; ") : ""}`);
  const missSel = d.steps.map((s) => s.sel).filter((s) => !hasSelectorId(s, false));
  ok(!missSel.length, `${d.id}: every step points at an element that exists${missSel.length ? " — missing " + missSel.join(", ") : ""}`);
  const badAct = d.actions.map(([a]) => a).filter((a) => !ACTIONS.includes(a));
  ok(!badAct.length, `${d.id}: every action is one the application can replay${badAct.length ? " — unknown " + badAct.join(", ") : ""}`);
  const docsOk = !!d.docs
    && docsHtml.includes(`id="en-${d.docs}"`)
    && docsHtml.includes(`id="es-${d.docs}"`)
    && extraLangs.every((l) => extraDocs[l].includes(`id="${l}-${d.docs}"`));
  ok(docsOk, `${d.id}: linked Docs anchor (#${d.docs}) exists in all 7 languages`);
  ok(!!d.cite, `${d.id}: names the work its data comes from`);
}

console.log("\n== the medium table gives back the rates it was built from");
{
  // Dietmair et al., PLoS ONE 7:e43394 (2012): Table 2 glucose 528, lactate -695, glutamine 29 µmol/gDW/h;
  // Table S1 dry mass 514 pg/cell; growth 0.028 1/h
  const p = parseExoTable(EXO_TABLE, { gdwPerCell: 514e-12 });
  ok(!p.errors.length && !p.warnings.length, `parses clean (${p.errors.length} errors, ${p.warnings.length} warnings)`);
  const rate = new Map(groupForRates(p).map((g) => [g.metabolite, ratesFromReplicates(g).rate]));
  near(rate.get("glucose"), 0.528, 0.004, "glucose uptake is the published 0.528");
  near(rate.get("L-lactate"), -0.695, 0.005, "lactate secretion is the published -0.695");
  near(rate.get("glutamine"), 0.029, 0.002, "glutamine uptake is the published 0.029");
  near(growthRate(p.series[0].times, p.series[0].biomass).mu, 0.028, 0.0005, "growth is the published 0.028 1/h");
}

console.log("\n== the intracellular table is the published one, value by value");
{
  const src = JSON.parse(fs.readFileSync(R + "tools/validation/dietmair2012_intracellular.json"));
  const by = new Map(src.metabolites.map((m) => [m.label, m]));
  const lines = MTB_TABLE.split("\n");
  ok(lines[1].startsWith("group\t"), "a group row names the condition of each sample");
  let wrong = [];
  for (const line of lines.slice(2)) {
    const c = line.split("\t"), m = by.get(c[0]);
    if (!m) { wrong.push(`${c[0]} is not in the source`); continue; }
    const same = (vals, ref) => vals.every((x, i) => Math.abs(x - ref[i]) / Math.max(1e-9, Math.abs(ref[i])) < 0.005);
    if (!same(c.slice(1, 5).map(Number), m.producer.map((x) => x.mean))) wrong.push(`${c[0]} producer`);
    if (!same(c.slice(5, 9).map(Number), m.nonProducer.map((x) => x.mean))) wrong.push(`${c[0]} non-producer`);
  }
  ok(!wrong.length, `all ${lines.length - 2} metabolites match tools/validation/dietmair2012_intracellular.json${wrong.length ? " — wrong: " + wrong.join(", ") : ""}`);

  // and it still shows what the demo says it shows
  const t = parseWideTable(MTB_TABLE);
  const g = splitByPattern(t, "producer", "parent");
  ok(g.nA === 4 && g.nB === 4 && !g.overlap, `the group row splits it ${g.nA} against ${g.nB}, with no sample on both sides`);
  const c = contrast(g.rows);
  const byLabel = new Map(c.rows.map((r) => [r.label, r.log2fc]));
  const atpAmp = ratioChanges([["ATP", "AMP", ""]], byLabel)[0];
  ok(atpAmp && atpAmp.dlog2 < -0.9, `the energy charge falls, as the demo tells the reader: ATP/AMP ${atpAmp?.dlog2.toFixed(2)} log2 (x${atpAmp?.fold.toFixed(2)})`);
  const lacPyr = ratioChanges([["LAC", "PYR", ""]], byLabel)[0];
  ok(lacPyr && Math.abs(lacPyr.dlog2 + 0.3) < 0.15, `and lactate/pyruvate matches the whole 52-metabolite panel (-0.30): ${lacPyr?.dlog2.toFixed(2)}`);
}

console.log("\n== the RNA-seq TPM table contains real Human-GEM HCT116 genes (33 expressed + 10 silenced)");
{
  const humanModel = JSON.parse(fs.readFileSync(R + "data/human-gem/model.json", "utf8"));
  const validGenes = new Set(humanModel.genes.map((g) => g[1]));
  const rows = RNA_TABLE.trim().split("\n").slice(1).map((l) => {
    const [gene, tpm] = l.split("\t");
    return { gene, tpm: Number(tpm) };
  });
  const missing = rows.filter((r) => !validGenes.has(r.gene)).map((r) => r.gene);
  const on = rows.filter((r) => r.tpm >= 1).length;
  const off = rows.filter((r) => r.tpm < 1).length;
  ok(!missing.length, `all ${rows.length} HCT116 genes exist in Human-GEM${missing.length ? " — missing: " + missing.join(", ") : ""}`);
  ok(on === 33 && off === 10, `exact split at 1 TPM cutoff (${on} expressed >= 1 TPM, ${off} silenced < 1 TPM)`);
}

console.log("\n== the regulatory RNA-seq table ranks ArcA #1 on E. coli RegulonDB");
{
  const ecoliModel = JSON.parse(fs.readFileSync(R + "data/ecoli-gem/model.json", "utf8"));
  const ecoliReg = JSON.parse(fs.readFileSync(R + "data/ecoli-gem/regulons.json", "utf8"));
  const p = parseRnaDiffTable(REG_TABLE);
  ok(p.map.size >= 400 && !p.warnings.length, `parses clean (${p.map.size} genes, ${p.warnings.length} warnings)`);
  const res = inferRegulatorsMultiOmic(ecoliModel, null, ecoliReg, null, null, p.map, { minTargets: 5 });
  const top = res.tfs[0];
  ok(top && top.tf === "ArcA" && top.zRna < -4 && top.qRna < 1e-3,
    `ArcA ranks #1 with negative activity (${top?.tf}: t=${top?.zRna?.toFixed(2)}, q=${top?.qRna?.toExponential(1)}, n=${top?.nObsRna})`);
}

console.log(bad ? `\n${bad} check(s) FAILED` : "\nall passed");
process.exit(bad ? 1 : 0);
