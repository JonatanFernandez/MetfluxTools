"""Expression levels of the model genes in each cell type / tissue, from the downloaded atlases.

Human-GEM (Ensembl gene IDs ENSG...):
  * Human Protein Atlas (HPA) cell line atlas, rna_celline.tsv.zip: nTPM per cell line (1206 lines).
    Cell line metadata (cancer or not, disease, site, NCI-60 membership) from Cellosaurus (cellline_meta.json,
    made by cellline_meta.py), plus 5 HPA derivatives that Cellosaurus does not list (LINE_OVERRIDES).
  * HPA single cell type atlas, rna_single_cell_type.tsv.zip: nCPM per cell type (pseudobulk of scRNA-seq).
  * HPA tissue consensus, rna_tissue_consensus.tsv.zip: nTPM per tissue (HPA + GTEx).
  * NCI-60 consensus: per gene, the median nTPM over the NCI-60 lines present in the HPA cell line atlas
    (pairs with the NCI-60 measured exchange fluxes, which are also a median over the panel).
Mouse-GEM (gene symbols; mapped to Ensembl ENSMUSG with Mouse-GEM's data/MGI_gene_ID_mapping.tsv):
  * Bgee, Mus_musculus_expr_simple.tsv.gz: present/absent calls per anatomical entity or cell type. Present
    -> the Bgee expression score (0-100), absent -> 0, no call -> unknown. Only entities with a call for at
    least 80 % of the model genes, and at least one of them absent, are kept.

Output: tools/cache/contexts/<slug>_expr.json
  { genes: [symbol...], unit, contexts: [{ id, name, group, src, info, levels: [x | null per gene] }] }
"""
import csv, gzip, io, json, statistics, zipfile
from collections import defaultdict

ROOT = "/home/ouuiea/Projects/MetfluxLive/"
EX = ROOT + "tools/cache/expression/"
OUT = ROOT + "tools/cache/contexts/"

LINE_OVERRIDES = {  # HPA cell lines with no Cellosaurus entry of their own (described on their HPA pages)
    "ASC2telo differentiated": {"category": "Telomerase immortalized cell line", "site": "Adipose tissue", "note": "ASC52telo differentiated to adipocytes"},
    "BJ hTERT+ SV40 Large T+": {"category": "Transformed cell line", "site": "Foreskin, skin", "note": "BJ fibroblasts with hTERT and SV40 large T"},
    "BJ hTERT+ SV40 Large T+ RasG12V": {"category": "Transformed cell line", "site": "Foreskin, skin", "note": "BJ fibroblasts with hTERT, SV40 large T and HRAS G12V (tumorigenic)"},
    "HHSteC": {"category": "Primary cell", "site": "Liver", "note": "Primary hepatic stellate cells"},
    "HSkMC": {"category": "Primary cell", "site": "Skeletal muscle", "note": "Primary skeletal muscle cells"},
}


def read_zip_tsv(path):
    with zipfile.ZipFile(path) as z:
        f = io.TextIOWrapper(z.open(z.namelist()[0]), encoding="utf-8")
        r = csv.reader(f, delimiter="\t")
        head = next(r)
        for row in r: yield dict(zip(head, row))


def slugify(s):
    return "".join(c if c.isalnum() else "_" for c in s.lower()).strip("_")


def human():
    model = json.load(open(ROOT + "data/human-gem/model.json"))
    ens = [g[0] for g in model["genes"]]; pos = {e: i for i, e in enumerate(ens)}
    n = len(ens)
    meta = json.load(open(EX + "cellline_meta.json"))
    contexts = []

    def collect(path, ctxcol, valcol):
        data = defaultdict(lambda: [None] * n)
        for row in read_zip_tsv(path):
            i = pos.get(row["Gene"])
            if i is not None: data[row[ctxcol]][i] = float(row[valcol])
        return data

    lines = collect(EX + "rna_celline.tsv.zip", "Cell line", "nTPM")
    for name, lv in sorted(lines.items(), key=lambda kv: kv[0].lower()):
        m = meta.get(name) or {}
        o = LINE_OVERRIDES.get(name, {})
        cat = m.get("category") or o.get("category", "")
        cancer = cat == "Cancer cell line"
        info = {"category": cat, "disease": m.get("disease", ""), "site": m.get("site") or o.get("site", ""),
                "sex": m.get("sex", ""), "cellosaurus": m.get("ac", ""), "nci60": bool(m.get("nci60")),
                "note": o.get("note", ""), "problematic": m.get("problematic", "")}
        contexts.append({"id": "hpa_line:" + name, "name": name, "group": "lineCancer" if cancer else "lineOther",
                         "src": "hpaLine", "info": info, "levels": lv})
    # NCI-60 consensus: median over the NCI-60 lines present in HPA
    nci = [name for name in lines if meta.get(name, {}).get("nci60")]
    med = [None] * n
    for i in range(n):
        xs = [lines[name][i] for name in nci if lines[name][i] is not None]
        med[i] = statistics.median(xs) if xs else None
    contexts.insert(0, {"id": "nci60", "name": "NCI-60", "group": "consensus", "src": "hpaLine",
                        "info": {"lines": sorted(nci), "category": "Cancer cell line"}, "levels": med})
    for name, lv in sorted(collect(EX + "rna_single_cell_type.tsv.zip", "Cell type", "nCPM").items()):
        contexts.append({"id": "hpa_cell:" + name, "name": name, "group": "cellType", "src": "hpaCell", "info": {}, "levels": lv})
    for name, lv in sorted(collect(EX + "rna_tissue_consensus.tsv.zip", "Tissue", "nTPM").items()):
        contexts.append({"id": "hpa_tissue:" + name, "name": name, "group": "tissue", "src": "hpaTissue", "info": {}, "levels": lv})
    return {"genes": [g[1] for g in model["genes"]], "unit": "nTPM", "contexts": contexts}


def mouse():
    model = json.load(open(ROOT + "data/mouse-gem/model.json"))
    syms = [g[1] for g in model["genes"]]; pos = {s: i for i, s in enumerate(syms)}
    n = len(syms)
    with zipfile.ZipFile(ROOT + "Mouse-GEM-1.8.0.zip") as z:
        f = io.TextIOWrapper(z.open("Mouse-GEM-1.8.0/data/MGI_gene_ID_mapping.tsv"), encoding="utf-8-sig")
        r = csv.reader(f, delimiter="\t"); next(r)
        ens = {row[5]: pos[row[2]] for row in r if row[2] in pos and row[5]}
    data, names, calls = defaultdict(lambda: [None] * n), {}, defaultdict(int)
    for line in gzip.open(EX + "Mus_musculus_expr_simple.tsv.gz", "rt"):
        p = line.rstrip("\n").split("\t")
        i = ens.get(p[0])
        if i is None: continue
        ent = p[2]; names[ent] = p[3].strip('"')
        data[ent][i] = float(p[7]) if p[4] == "present" else 0.0
        calls[ent] += 1
    mapped = len(set(ens.values()))
    contexts = []
    for ent in sorted(data, key=lambda e: names[e].lower()):
        if calls[ent] < 0.8 * mapped: continue
        if not any(x == 0.0 for x in data[ent]): continue  # no gene called absent: same as the generic model
        group = "cellType" if ent.startswith("CL:") else "tissue"
        contexts.append({"id": "bgee:" + ent, "name": names[ent], "group": group, "src": "bgee",
                         "info": {"entity": ent}, "levels": data[ent]})
    return {"genes": syms, "unit": "Bgee score", "contexts": contexts}


def bgee(slug, species, zipname, idmap):
    """Any other model, from the Bgee calls of its species. Model genes are matched by the gene name Bgee lists
    (the model uses the same official symbols) and, failing that, by the IDs that idmap(zip) gives (symbol -> Bgee ID)."""
    model = json.load(open(ROOT + f"data/{slug}/model.json"))
    syms = [g[1] for g in model["genes"]]; pos = {s: i for i, s in enumerate(syms)}
    n = len(syms)
    with zipfile.ZipFile(ROOT + zipname) as z:
        byid = {gid: pos[sym] for sym, gid in idmap(z) if sym in pos}
    data, names, calls = defaultdict(lambda: [None] * n), {}, defaultdict(int)
    for line in gzip.open(EX + f"{species}_expr_simple.tsv.gz", "rt"):
        p = line.rstrip("\n").split("\t")
        i = byid.get(p[0])
        if i is None: i = pos.get(p[1].strip('"'))
        if i is None: continue
        ent = p[2]; names[ent] = p[3].strip('"')
        data[ent][i] = float(p[7]) if p[4] == "present" else 0.0
        calls[ent] += 1
    mapped = len({i for v in data.values() for i, x in enumerate(v) if x is not None})
    print(slug, "model genes with Bgee calls:", mapped, "of", n)
    contexts = []
    for ent in sorted(data, key=lambda e: names[e].lower()):
        if calls[ent] < 0.8 * mapped: continue
        if not any(x == 0.0 for x in data[ent]): continue
        group = "cellType" if ent.startswith("CL:") else "tissue"
        contexts.append({"id": "bgee:" + ent, "name": names[ent], "group": group, "src": "bgee",
                         "info": {"entity": ent}, "levels": data[ent]})
    return {"genes": syms, "unit": "Bgee score", "contexts": contexts}


def tsv_rows(z, suffix):
    name = next(x for x in z.namelist() if x.endswith(suffix))
    f = io.TextIOWrapper(z.open(name), encoding="utf-8-sig")
    r = csv.DictReader(f, delimiter="\t")
    return list(r)


def ortho_ids(prefix):  # human2<Species>Orthologs.tsv: toSymbol -> toGeneId without its database prefix
    return lambda z: [(r["toSymbol"], r["toGeneId"].split(":", 1)[1]) for r in tsv_rows(z, "Orthologs.tsv") if r["toGeneId"].startswith(prefix)]


OTHERS = {
    "rat-gem": ("Rattus_norvegicus", "Rat-GEM-1.8.0.zip", lambda z: [(r["SYMBOL"], r["ENSEMBL_ID"]) for r in tsv_rows(z, "GENES_RAT.tsv") if r["ENSEMBL_ID"]]),
    "zebrafish-gem": ("Danio_rerio", "Zebrafish-GEM-1.6.0.zip", lambda z: []),
    "fruitfly-gem": ("Drosophila_melanogaster", "Fruitfly-GEM-1.7.0.zip", ortho_ids("FB:")),
    "worm-gem": ("Caenorhabditis_elegans", "Worm-GEM-1.8.0.zip", ortho_ids("WB:")),
}


if __name__ == "__main__":
    import os, sys
    os.makedirs(OUT, exist_ok=True)
    jobs = [(s, (lambda s=s: bgee(s, *OTHERS[s]))) for s in OTHERS]
    jobs = [("human-gem", human), ("mouse-gem", mouse)] + jobs
    if len(sys.argv) > 1: jobs = [j for j in jobs if j[0] in sys.argv[1:]]
    for slug, fn in jobs:
        d = fn()
        json.dump(d, open(OUT + f"{slug}_expr.json", "w"))
        from collections import Counter
        print(slug, len(d["contexts"]), "contexts", Counter(c["group"] for c in d["contexts"]))
