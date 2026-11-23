#!/usr/bin/env python3
"""Structures for metabolites of yeast-GEM and iML1515 that have no Human-GEM counterpart, from MetaNetX/MNXref.

Reads tools/cache/src/<slug>/model/metabolites.tsv (written by convert_foreign.py), resolves each metabolite's MetaNetX
ids (deprecated ids followed through chem_depr.tsv) and, failing those, its ChEBI ids (chem_xref.tsv), and takes the
SMILES of chem_prop.tsv. Only structures whose MNXref reference is ChEBI, SwissLipids or LIPID MAPS (all CC BY 4.0) are
kept, so that no structure of a source with redistribution limits (KEGG, MetaCyc, BiGG, HMDB) enters the app.
The large MetaNetX files are streamed and filtered, not stored.

Output tools/cache/mnx_structures.json: {slug: {base: {mnx, ref, smiles}}}; convert_foreign.py reads it.
Usage: fetch_mnx_structures.py   (then rerun convert_foreign.py and build_model.py for yeast-gem and ecoli-gem)
"""
import csv, json, subprocess
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
URL = "https://www.metanetx.org/ftp/latest/"   # MetaNetX/MNXref 4.5 when this was written (CC BY 4.0)
OPEN = ("chebi:", "slm:", "lipidmapsM:")
SLUGS = ["ecoli-gem", "yeast-gem"]


def stream(name, keys):
    """Rows of a MetaNetX file whose first column is in keys (curl + awk, the file is not kept)."""
    keyfile = ROOT / "tools/cache/_mnx_keys.txt"
    keyfile.write_text("\n".join(sorted(keys)))
    awk = "awk -F'\\t' 'NR==FNR{k[$1]=1;next} /^#/{next} ($1 in k){print}' " + str(keyfile) + " -"
    out = subprocess.run(f"curl -sf {URL}{name} | {awk}", shell=True, capture_output=True, text=True, check=True).stdout
    keyfile.unlink()
    return [l.split("\t") for l in out.splitlines()]


def main():
    hs = json.loads((ROOT / "tools/cache/human-gem/smiles.json").read_text())
    todo = {}
    for slug in SLUGS:
        rows = csv.DictReader(open(ROOT / f"tools/cache/src/{slug}/model/metabolites.tsv"), delimiter="\t")
        todo[slug] = {r["metsNoComp"]: r for r in rows if not r["metSmiles"] and r["metsNoComp"] not in hs}
    split = lambda s: [x.strip() for x in s.split(";") if x.strip()]
    mnx = {x for t in todo.values() for r in t.values() for x in split(r["metMetaNetXID"])}
    chebi = {x.replace("CHEBI:", "chebi:") for t in todo.values() for r in t.values() for x in split(r["metChEBIID"])}

    depr = defaultdict(set)
    for l in subprocess.run(f"curl -sf {URL}chem_depr.tsv", shell=True, capture_output=True, text=True, check=True).stdout.splitlines():
        if not l.startswith("#"):
            a, b = l.split("\t")[:2]
            depr[a].add(b)

    def current(x, seen=frozenset()):
        if x in seen or x not in depr:
            return {x}
        return set().union(*(current(y, seen | {x}) for y in depr[x]))

    c2m = defaultdict(set)
    for f in stream("chem_xref.tsv", chebi):
        c2m[f[0]].add(f[1])
    cand = set().union(*(current(x) for x in mnx)) | set().union(*c2m.values())
    prop = {f[0]: {"ref": f[2], "smiles": f[8] if len(f) > 8 else ""} for f in stream("chem_prop.tsv", cand)}

    out = {}
    for slug, bases in todo.items():
        res, why = {}, defaultdict(int)
        for b, r in bases.items():
            ids = [y for x in split(r["metMetaNetXID"]) for y in sorted(current(x))]
            ids += [y for x in split(r["metChEBIID"]) for y in sorted(c2m.get(x.replace("CHEBI:", "chebi:"), ()))]
            reason = "no MetaNetX or ChEBI id"
            for c in dict.fromkeys(ids):
                p = prop.get(c)
                if not p or not p["smiles"]:
                    reason = "no structure in MetaNetX"
                elif not p["ref"].startswith(OPEN):
                    reason = "structure from a source with redistribution limits"
                else:
                    res[b] = {"mnx": c, "ref": p["ref"], "smiles": p["smiles"]}
                    break
            else:
                why[reason] += 1
        out[slug] = res
        print(f"{slug}: {len(bases)} metabolites without a Human-GEM structure, {len(res)} found; missing: {dict(why)}")
    (ROOT / "tools/cache/mnx_structures.json").write_text(json.dumps(out, indent=0))


if __name__ == "__main__":
    main()
