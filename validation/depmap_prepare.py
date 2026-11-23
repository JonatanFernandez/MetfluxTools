"""DepMap Public 26Q3 -> tools/validation/depmap/depmap26q3_metabolic.json, for the metabolic genes of Human-GEM and the
cell lines of the app's catalogue (HPA lines matched by RRID = Cellosaurus accession in Model.csv).
Files in tools/validation/depmap/ (downloaded by hand from depmap.org): Model.csv, CRISPRGeneEffect.csv (Chronos),
CRISPRGeneDependency.csv (probability of dependency), CRISPRInferredCommonEssentials.csv.
    python3 tools/validation/depmap_prepare.py"""
import csv, json, os
R = os.path.join(os.path.dirname(__file__), "..", "..")
D = os.path.join(R, "tools", "validation", "depmap")
model = json.load(open(os.path.join(R, "data", "human-gem", "model.json")))
syms = {g[1] for g in model["genes"] if g[1]}
ctx = json.load(open(os.path.join(R, "data", "human-gem", "contexts.json")))["contexts"]
by_rrid = {c["info"]["cellosaurus"]: c for c in ctx if (c.get("info") or {}).get("cellosaurus")}
meta = {}
for r in csv.DictReader(open(os.path.join(D, "Model.csv"), encoding="utf-8")):
    c = by_rrid.get(r["RRID"])
    if c: meta[r["ModelID"]] = {"context": c["id"], "name": r["CellLineName"], "rrid": r["RRID"], "lineage": r["OncotreeLineage"], "disease": r["OncotreePrimaryDisease"], "media": r["FormulationID"]}
def matrix(fname):
    f = open(os.path.join(D, fname)); rd = csv.reader(f); head = next(rd)
    cols = [(i, h.split(" (")[0]) for i, h in enumerate(head) if i and h.split(" (")[0] in syms]
    out = {}
    for row in rd:
        if row[0] in meta: out[row[0]] = [None if row[i] in ("", "NA") else round(float(row[i]), 3) for i, _ in cols]
    return [s for _, s in cols], out
genes, effect = matrix("CRISPRGeneEffect.csv")
genes2, prob = matrix("CRISPRGeneDependency.csv")
assert genes == genes2
common = sorted({l.strip().split(" (")[0] for l in open(os.path.join(D, "CRISPRInferredCommonEssentials.csv")) if l.strip() and l.strip() != "Essentials"} & syms)
lines = {a: meta[a] for a in effect}
json.dump({"source": "DepMap Public 26Q3: CRISPRGeneEffect (Chronos), CRISPRGeneDependency, CRISPRInferredCommonEssentials, Model.csv; lines matched to the app's HPA cell lines by RRID",
           "genes": genes, "commonEssential": common, "lines": lines, "effect": effect, "prob": prob}, open(os.path.join(D, "depmap26q3_metabolic.json"), "w"))
print(f"catalogue lines in Model.csv: {len(meta)}; with CRISPR screen: {len(lines)}; metabolic genes screened: {len(genes)} of {len(syms)}; common essentials among them: {len(common)}")
