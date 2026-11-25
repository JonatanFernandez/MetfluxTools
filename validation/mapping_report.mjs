// How well do the metabolite labels of real published tables land on the models? One number per dataset,
// so that "should we invest in synonyms?" is a measurement and not an opinion (js/core/metnames.js).
//
// For each dataset it reports how many labels resolve, split by the layer that answered — an identifier,
// the model's own name, one of the model's abbreviations, or a naming convention — and lists what is left.
// The unresolved labels are split in two, because they call for different things:
//   · compounds the model simply does not contain (the lipids and plant compounds of an untargeted panel):
//     nothing can be done, and the interface should say so rather than imply a gap in our mapping;
//   · compounds that look like central metabolism and should have matched: those are the ones worth curating.
//
//   node tools/validation/mapping_report.mjs
// Datasets that are not downloaded are skipped with a note (tools/validation/fetch_metabolomics.sh,
// prepare_dietmair.py).
import fs from "node:fs";
import zlib from "node:zlib";
import { buildMetIndex, resolveMetabolite, suggestMetabolites } from "../../js/core/metnames.js";

const R = new URL("../../", import.meta.url).pathname;
const models = {};
const modelOf = (slug) => (models[slug] ??= JSON.parse(fs.readFileSync(`${R}data/${slug}/model.json`)));
const indexOf = (slug) => ((indexOf.c ??= {})[slug] ??= buildMetIndex(modelOf(slug)));

// labels that are plainly not metabolism as the models know it: lipid species, plant and drug compounds,
// internal standards. Counted apart so they do not look like a mapping failure.
const NOT_IN_MODELS = [
  /^\[(?:FA|GL|GP|SP|ST|PR|SL|PK|CE|PC|PE|PI|PS|PG|TG|DG|MG|Fv|Fl)[\s(\]]/i, // LIPID MAPS style annotation
  /ISTD|internal standard/i,
  /^\[?[A-Z]{2}\s?\(\d+:\d+/,                     // "PC (16:0/18:2)"
  /\d+[ZE][,)]|\(\d+[ZE]/,                           // an acyl chain with double-bond geometry
  /-sn-glycero|sphing|ceramide|cholesteryl|-oyl\)/i,   // lipid species written out in full
  /\(\d+:\d+\)/,                                     // "(17:0/20:5)"
];
const notInModels = (l) => NOT_IN_MODELS.some((rx) => rx.test(l));
// rows of the sheet that are headers or sample labels, not metabolites
const NOT_A_METABOLITE = /^(METABOLITE|Metabolites?\b|Samples?|Factors?|#|Log10|Sdhb|BY4742|YGR023)|^\s*$/i;

/** Read the metabolite labels of a Metabolomics Workbench mwTab file. */
function mwtabLabels(path) {
  const txt = fs.readFileSync(path, "utf8");
  const block = txt.split(/MS_METABOLITE_DATA_START|NMR_BINNED_DATA_START/)[1];
  if (!block) return null;
  const lines = block.split(/\r?\n/).slice(1);
  const out = [];
  for (const line of lines) {
    if (/^(MS_METABOLITE_DATA_END|NMR_BINNED_DATA_END|#)/.test(line)) break;
    const first = line.split("\t")[0]?.trim();
    if (first && !NOT_A_METABOLITE.test(first)) out.push(first);
  }
  return out;
}

const datasets = [
  // NOTE: this uses the sheet's own three-letter codes (ASP, GLU, SER…) on purpose, to measure what the
  // automatic cascade achieves alone. Those codes are genuinely ambiguous — ASP is aspartate or asparagine —
  // and prepare_dietmair.py carries a hand-written map, which brings that dataset to 51 of 52. A low number
  // here therefore means "the codes need a per-dataset map", not "the normaliser failed".
  { name: "Dietmair 2012, HEK293 producer (targeted, raw codes)", slug: "human-gem", labels: () => {
      const p = `${R}tools/validation/dietmair2012_intracellular.json`;
      if (!fs.existsSync(p)) return null;
      // the sheet's own abbreviations, not the names prepare_dietmair.py already mapped by hand
      return JSON.parse(fs.readFileSync(p)).metabolites.map((m) => m.label);
    } },
  { name: "Cardaci 2015, Sdhb knockout (untargeted)", slug: "mouse-gem", labels: () => {
      const p = `${R}tools/validation/metabolomics/cardaci2015_labels.txt`;
      if (!fs.existsSync(p)) return null;
      return fs.readFileSync(p, "utf8").split("\n").map((s) => s.trim()).filter((s) => s && !NOT_A_METABOLITE.test(s));
    } },
  { name: "ST002076, E. coli regulator knockouts", slug: "ecoli-gem", mw: "ST002076" },
  { name: "ST004848, mouse methionine restriction", slug: "mouse-gem", mw: "ST004848" },
  { name: "ST004339, rat INS-1 + BSO", slug: "rat-gem", mw: "ST004339" },
  { name: "ST001609, yeast mtl1 deletion", slug: "yeast-gem", mw: "ST001609" },
  { name: "ST003960, HEK293 spent medium", slug: "human-gem", mw: "ST003960" },
];

const total = { n: 0, id: 0, name: 0, alias: 0, rule: 0, absent: 0, missed: 0 };
for (const d of datasets) {
  let labels = null;
  if (d.mw) {
    const p = `${R}tools/validation/metabolomics/${d.mw}.mwtab`;
    labels = fs.existsSync(p) ? mwtabLabels(p) : null;
  } else labels = d.labels();
  if (!labels || !labels.length) { console.log(`\n== ${d.name}\n   not downloaded, skipped`); continue; }
  const index = indexOf(d.slug);
  const by = { id: [], name: [], alias: [], rule: [], none: [] };
  for (const l of labels) {
    const { via } = resolveMetabolite(l, index);
    by[via || "none"].push(l);
  }
  const absent = by.none.filter(notInModels);
  const missed = by.none.filter((l) => !notInModels(l));
  const resolved = labels.length - by.none.length;
  console.log(`\n== ${d.name}  (${d.slug})`);
  console.log(`   ${resolved} of ${labels.length} resolved (${Math.round((100 * resolved) / labels.length)} %)`
    + `: ${by.id.length} by identifier, ${by.name.length} by name, ${by.alias.length} by an abbreviation of the model, ${by.rule.length} by a naming convention`);
  if (by.rule.length) console.log(`   conventions: ${by.rule.slice(0, 6).join(", ")}${by.rule.length > 6 ? ` … and ${by.rule.length - 6} more` : ""}`);
  console.log(`   unresolved: ${absent.length} not in the model at all (lipid species, plant compounds, standards), ${missed.length} worth looking at`);
  if (missed.length) {
    for (const l of missed.slice(0, 8)) {
      const s = suggestMetabolites(l, index, 2).map((x) => `${x.name} (${x.score.toFixed(2)})`).join(", ");
      console.log(`      ${l.slice(0, 44).padEnd(44)} closest: ${s || "nothing similar"}`);
    }
    if (missed.length > 8) console.log(`      … and ${missed.length - 8} more`);
  }
  total.n += labels.length; total.id += by.id.length; total.name += by.name.length;
  total.alias += by.alias.length; total.rule += by.rule.length; total.absent += absent.length; total.missed += missed.length;
}

if (total.n) {
  const resolved = total.id + total.name + total.alias + total.rule;
  console.log(`\n== all datasets together\n   ${resolved} of ${total.n} labels resolved (${Math.round((100 * resolved) / total.n)} %)`);
  console.log(`   by identifier ${total.id}, by name ${total.name}, by an abbreviation ${total.alias}, by a convention ${total.rule}`);
  console.log(`   ${total.absent} are compounds the models do not contain; ${total.missed} are candidates for curation`);
  console.log(`   Of what could possibly be resolved (${resolved + total.missed}), ${Math.round((100 * resolved) / (resolved + total.missed))} % is.`);
  console.log(`   Read the Dietmair line with its note above: it is measured on three-letter codes with no map.`);
}
