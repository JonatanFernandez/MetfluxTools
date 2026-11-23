// Optimum of each candidate objective in the standard medium (no measured fluxes), per model. An exchange
// objective has its own uptake closed, as the simulator does. Usage: node tools/test_objectives.mjs <slug>... [-- RXN...]
import fs from "node:fs";
import loadHighs from "../vendor/highs/highs.mjs";
import { FbaProblem, packNetwork } from "../js/core/fba.js";
import { presetMedium, findExchanges, applyMedium } from "../js/core/medium.js";
const args = process.argv.slice(2), cut = args.indexOf("--");
const slugs = cut < 0 ? args : args.slice(0, cut);
const RX = cut < 0 ? ["MAR11420", "MAR09852", "MAR09066", "MAR11965", "MAR09063", "MAR09061", "MAR06916", "MAR09451", "MAR09397", "MAR09365", "MAR09433", "MAR00619", "MAR11427", "MAR09354", "MAR09857"] : args.slice(cut + 1);
const highs = await loadHighs();
for (const slug of slugs) {
  const model = JSON.parse(fs.readFileSync(new URL(`../data/${slug}/model.json`, import.meta.url)));
  const P = new FbaProblem(highs, packNetwork(model));
  const ex = findExchanges(model), med = presetMedium(model, "default");
  const base = applyMedium(model, ex, med, model.rxns.map(r => r.lb), model.rxns.map(r => r.ub));
  const out = [];
  for (const id of RX) {
    const j = model.rxns.findIndex(r => r.id === id);
    if (j < 0) { out.push(`${id} -`); continue; }
    const lb = base.lb.slice(), ub = base.ub.slice();
    if (model.rxns[j].m.length === 1) lb[j] = 0;
    if (id === "MAR13082" || id.startsWith("MFL_")) ub[j] = 1000; // closed in the model, opened as objective
    const r = P.solve({ lb, ub, objective: j, maximize: true, parsimonious: false });
    out.push(`${id} ${r.status === "optimal" ? r.objective.toPrecision(3) : r.status}`);
  }
  console.log(slug.padEnd(14), out.join(" | "));
}
