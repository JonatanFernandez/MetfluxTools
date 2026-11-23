"""Match the Human Protein Atlas cell lines to Cellosaurus and keep what the app needs to describe them.

Input : tools/cache/expression/rna_celline.tsv.zip (HPA), tools/cache/expression/cellosaurus.txt (Cellosaurus)
Output: tools/cache/expression/cellline_meta.json  {hpa name: {ac, id, category, disease, site, sex, nci60}}
"""
import json, re, sys, zipfile, io
D = "/home/ouuiea/Projects/MetfluxLive/tools/cache/expression/"
norm = lambda s: re.sub(r"[^a-z0-9]", "", s.lower())

names = set()
with zipfile.ZipFile(D + "rna_celline.tsv.zip") as z:
    f = io.TextIOWrapper(z.open(z.namelist()[0]), encoding="utf-8"); next(f)
    for line in f: names.add(line.split("\t", 3)[2])

entries, cur = [], None
for line in open(D + "cellosaurus.txt", encoding="utf-8"):
    tag, val = line[:2], line[5:].rstrip("\n")
    if tag == "ID": cur = {"id": val, "sy": [], "di": [], "ox": [], "cc": [], "ca": "", "sx": ""}; entries.append(cur)
    elif cur is None: continue
    elif tag == "AC": cur["ac"] = val
    elif tag == "SY": cur["sy"] += [s.strip() for s in val.split(";") if s.strip()]
    elif tag == "DI": cur["di"].append(val.split(";")[-1].strip())
    elif tag == "OX": cur["ox"].append(val)
    elif tag == "CA": cur["ca"] = val
    elif tag == "SX": cur["sx"] = val
    elif tag == "CC": cur["cc"].append(val)
human = [e for e in entries if any("9606" in o for o in e["ox"])]
by_id, by_sy = {}, {}
for e in human:
    by_id.setdefault(norm(e["id"]), []).append(e)
    for s in e["sy"]: by_sy.setdefault(norm(s), []).append(e)

def site(e):
    for c in e["cc"]:
        if c.startswith("Derived from site:"):
            parts = [p.strip() for p in c[len("Derived from site:"):].split(";")]
            return parts[1] if len(parts) > 1 else parts[0]
    return ""

out, miss, amb = {}, [], []
for n in sorted(names):
    k = norm(n); cand = by_id.get(k) or by_sy.get(k) or []
    # prefer a parental line (no problematic/contaminated note) when several share the name
    cand = sorted(cand, key=lambda e: (n.lower() not in [x.lower() for x in [e["id"]] + e["sy"]], any("Problematic cell line" in c for c in e["cc"]), len(e["id"])))
    if not cand: miss.append(n); continue
    if len({e["ac"] for e in cand}) > 1: amb.append((n, [e["id"] for e in cand][:3]))
    e = cand[0]
    out[n] = {"ac": e["ac"], "id": e["id"], "category": e["ca"], "disease": e["di"][0] if e["di"] else "",
              "site": site(e), "sex": e["sx"], "nci60": any("NCI-60" in c for c in e["cc"]),
              "problematic": next((c for c in e["cc"] if c.startswith("Problematic cell line")), "")}
json.dump(out, open(D + "cellline_meta.json", "w"), indent=0)
from collections import Counter
print(len(names), "HPA lines;", len(out), "matched;", len(miss), "missing:", miss[:20])
print("ambiguous:", len(amb), amb[:8])
print(Counter(v["category"] for v in out.values()).most_common())
print("NCI-60 in HPA:", sum(v["nci60"] for v in out.values()), "problematic:", sum(bool(v["problematic"]) for v in out.values()))
