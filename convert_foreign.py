#!/usr/bin/env python3
"""Convert a model that does not use Human-GEM identifiers (yeast-GEM, BiGG models such as iML1515) into a
standard-GEM folder that build_model.py reads: tools/cache/src/<slug>/model/<slug>.yml (+ metabolites.tsv).

Metabolites are matched to Human-GEM by shared identifiers (MetaNetX, then BiGG, ChEBI, KEGG, then exact name), one
to one. A matched metabolite takes the Human-GEM name and identifier of its compartment-free form, so that the
media, tracers, cofactor lines and structures of the app (which go by Human-GEM names) work unchanged; the original
name is kept in the identifier table (metsOriginal). Unmatched metabolites keep their own name; their structure comes
from the model's own SMILES annotation when it has one, otherwise from MetaNetX (tools/fetch_mnx_structures.py, only
structures that MNXref takes from ChEBI, SwissLipids or LIPID MAPS).

Usage:
  convert_foreign.py yeast-GEM-9.1.1.zip yeast-gem
  convert_foreign.py ecoli-src/iML1515.json ecoli-gem
"""
import csv, io, json, re, sys, zipfile
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent


def load(src):
    p = Path(src)
    if p.suffix == ".json":
        return json.loads(p.read_text())
    with zipfile.ZipFile(p) as z:
        name = next(n for n in z.namelist() if re.search(r"model/[^/]+\.ya?ml$", n))
        text = z.read(name).decode("utf-8")

    class loader(getattr(yaml, "CSafeLoader", yaml.SafeLoader)):
        pass
    loader.yaml_implicit_resolvers = {ch: [(t, rx) for t, rx in lst if t != "tag:yaml.org,2002:bool"]
                                      for ch, lst in loader.yaml_implicit_resolvers.items()}

    def omap(n):
        if isinstance(n, list):
            if n and all(isinstance(x, tuple) and len(x) == 2 for x in n):
                return {k: omap(v) for k, v in n}
            return [omap(x) for x in n]
        if isinstance(n, dict):
            return {k: omap(v) for k, v in n.items()}
        return n
    m = omap(yaml.load(text, Loader=loader))
    if isinstance(m, list):
        m = {k: v for d in m for k, v in d.items()}
    return m


def ids(ann, key):
    v = (ann or {}).get(key)
    if v is None:
        return []
    v = v if isinstance(v, list) else [v]
    out = []
    for x in v:
        for y in str(x).split(";"):
            y = y.strip()
            if y:
                out.append(y.replace("CHEBI:", "").replace("chebi:", ""))
    return out


def human_index():
    h = json.loads((ROOT / "data/human-gem/model.json").read_text())
    idx = {k: {} for k in ("mnx", "bigg", "chebi", "kegg", "name")}
    for b in h["bases"]:
        for k in ("mnx", "bigg", "chebi", "kegg"):
            for x in str(b.get(k, "")).split(";"):
                x = x.strip().replace("CHEBI:", "")
                if x:
                    idx[k].setdefault(x, b["id"])
        idx["name"].setdefault(b["n"].lower(), b["id"])
    names = {b["id"]: b["n"] for b in h["bases"]}
    return idx, names


def main():
    src, slug = sys.argv[1], sys.argv[2]
    m = load(src)
    idx, hname = human_index()
    meta = m.get("metaData") or {"id": m.get("id", slug), "version": m.get("version", "")}
    comps = m["compartments"]
    if isinstance(comps, list):
        comps = {k: v for d in comps for k, v in d.items()}

    # compartment-free form of each metabolite: by its id stem (BiGG: glc__D_c -> glc__D) or, when ids carry no
    # compartment (yeast-GEM s_0001), by name + identifiers
    def stem(met):
        c = met.get("compartment", "")
        if c and met["id"].endswith("_" + c):
            return met["id"][: -len(c) - 1]
        return "N:" + met.get("name", met["id"])
    groups = {}
    for met in m["metabolites"]:
        groups.setdefault(stem(met), []).append(met)

    mnx_file = ROOT / "tools/cache/mnx_structures.json"
    mnx_smiles = json.loads(mnx_file.read_text()).get(slug, {}) if mnx_file.exists() else {}
    claimed, base_of, rows, matched, used = set(), {}, [], 0, set()
    for key, mets in groups.items():
        ann = mets[0].get("annotation") or {}
        cand = None
        for k, akey in (("mnx", "metanetx.chemical"), ("bigg", "bigg.metabolite"), ("chebi", "chebi"), ("kegg", "kegg.compound")):
            for x in ids(ann, akey):
                b = idx[k].get(x)
                if b and b not in claimed:
                    cand = b
                    break
            if cand:
                break
        if not cand:
            b = idx["name"].get(mets[0].get("name", "").lower())
            if b and b not in claimed:
                cand = b
        if cand:
            claimed.add(cand)
            matched += 1
            base, name = cand, hname[cand]
        else:
            base = re.sub(r"[^A-Za-z0-9_]", "_", key[2:] if key.startswith("N:") else key)
            base = f"{slug.split('-')[0].upper()}_{base}"[:60]
            while base in used:  # two names that differ only in punctuation
                base += "_"
            name = mets[0].get("name", key)
        used.add(base)
        smi = ids(ann, "smiles")
        for met in mets:
            base_of[met["id"]] = (base, name)
            rows.append({"mets": met["id"], "metsNoComp": base, "metsOriginal": met.get("name", ""),
                         "metSmiles": smi[0] if smi and not cand else "",
                         "metSmilesMNX": mnx_smiles[base]["smiles"] if not smi and not cand and base in mnx_smiles else "",
                         "metKEGGID": ";".join(ids(ann, "kegg.compound")), "metChEBIID": ";".join("CHEBI:" + x for x in ids(ann, "chebi")),
                         "metBiGGID": ";".join(ids(ann, "bigg.metabolite")), "metMetaNetXID": ";".join(ids(ann, "metanetx.chemical"))})
    print(f"{slug}: {len(groups)} metabolites, {matched} matched to Human-GEM")

    # several compartment-free forms can still map to one base only once (claimed): ids stay unique per compartment
    mets_out = []
    for met in m["metabolites"]:
        base, name = base_of[met["id"]]
        mets_out.append({"id": met["id"], "name": name, "compartment": met.get("compartment"), "charge": met.get("charge", 0),
                         "formula": met.get("formula", "")})
    rxns_out = []
    for r in m["reactions"]:
        ss = r.get("subsystem") or "Unassigned"
        ss = ss[0] if isinstance(ss, list) and ss else ss
        ss = re.sub(r"^sce\d+\s+", "", str(ss)).strip() or "Unassigned"
        ann = r.get("annotation") or {}
        rr = {"id": r["id"], "name": r.get("name") or r["id"], "metabolites": r["metabolites"],
              "lower_bound": r.get("lower_bound", -1000), "upper_bound": r.get("upper_bound", 1000),
              "gene_reaction_rule": r.get("gene_reaction_rule", ""), "subsystem": ss}
        ec = ids(ann, "ec-code")
        if ec:
            rr["eccodes"] = ec
        if r.get("objective_coefficient"):
            rr["objective_coefficient"] = r["objective_coefficient"]
        rxns_out.append(rr)
    genes_out = [{"id": g["id"], "name": g.get("name") or g["id"]} for g in m.get("genes", [])]
    doc = {"metaData": {"id": meta.get("id", slug), "name": meta.get("name", ""), "version": str(meta.get("version", "")),
                        "taxonomy": meta.get("taxonomy", "")},
           "metabolites": mets_out, "reactions": rxns_out, "genes": genes_out, "compartments": comps}
    out = ROOT / "tools/cache/src" / slug / "model"
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{slug}.yml").write_text(json.dumps(doc))
    with open(out / "metabolites.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print("wrote", out)


if __name__ == "__main__":
    main()
