// Growth with the generic biomass (MAR13082, the app's default objective) in the standard medium, per model.
import fs from "node:fs";
import loadHighs from "../vendor/highs/highs.mjs";
import { FbaProblem, packNetwork } from "../js/core/fba.js";
import { presetMedium, findExchanges, applyMedium } from "../js/core/medium.js";
const highs = await loadHighs();
for (const slug of process.argv.slice(2)) {
  const model = JSON.parse(fs.readFileSync(new URL(`../data/${slug}/model.json`, import.meta.url)));
  const P = new FbaProblem(highs, packNetwork(model));
  const ex = findExchanges(model), med = presetMedium(model, process.env.MEDIUM || "default");
  const { lb, ub } = applyMedium(model, ex, med, model.rxns.map(r => r.lb), model.rxns.map(r => r.ub));
  let obj = model.rxns.findIndex(r => r.id === "MAR13082"); if (obj < 0) obj = model.rxns.findIndex(r => r.obj); else ub[obj] = 1000;
  const r = P.solve({ lb, ub, objective: obj, maximize: true, parsimonious: true });
  const top = [...ex].map(([mi, e]) => [model.bases[model.mets[mi].b].n, r.v[e.rxn] * e.dir]).filter(x => Math.abs(x[1]) > 0.05).sort((a, b) => b[1] - a[1]);
  console.log(slug, r.status, r.objective.toFixed(4), top.map(([n, v]) => `${n} ${v.toFixed(2)}`).join(", "));
}
