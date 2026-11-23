#!/usr/bin/env python3
"""Write data/models.json listing every built model (read by the web app at start-up)."""
import json
from pathlib import Path

data = Path(__file__).resolve().parent.parent / "data"
order = {"human-gem": 0, "mouse-gem": 1, "rat-gem": 2, "zebrafish-gem": 3, "fruitfly-gem": 4, "worm-gem": 5, "yeast-gem": 6, "ecoli-gem": 7}
# organism and per-model defaults read by the app: the medium of the usual cell culture of that organism (see docs,
# section "Models") and no measured exchange fluxes (the measured sets are human cells)
ORGANISM = {"human-gem": "Homo sapiens", "mouse-gem": "Mus musculus", "rat-gem": "Rattus norvegicus",
            "zebrafish-gem": "Danio rerio", "fruitfly-gem": "Drosophila melanogaster", "worm-gem": "Caenorhabditis elegans",
            "yeast-gem": "Saccharomyces cerevisiae S288C", "ecoli-gem": "Escherichia coli K-12 MG1655"}
DEFAULTS = {"zebrafish-gem": {"medium": "default", "measured": "none"},
            "fruitfly-gem": {"medium": "schneider", "measured": "none"},
            "worm-gem": {"medium": "default", "measured": "none"},
            "yeast-gem": {"medium": "verduyn", "measured": "none", "kind": "microbe"},
            "ecoli-gem": {"medium": "m9glc", "measured": "none", "kind": "microbe"}}
NAMES = {"yeast-gem": "Yeast-GEM", "ecoli-gem": "E. coli iML1515"}
VERSIONS = {"ecoli-gem": "(BiGG)"}  # the BiGG file's "version" is that of the JSON format  # display names when the folder name is not the model's
models = []
for d in sorted(data.iterdir(), key=lambda p: (order.get(p.name, 9), p.name)):
    need = ["model.json", "structures.json", "atommaps.json", "layout.json"]
    if d.is_dir() and all((d / f).exists() for f in need):
        meta = json.loads((d / "model.json").read_text())["meta"]
        pretty = "-".join(w.upper() if w == "gem" else w.capitalize() for w in d.name.split("-"))
        entry = {"slug": d.name, "name": NAMES.get(d.name, pretty), "version": VERSIONS.get(d.name, meta.get("version", "")),
                 "description": meta.get("name", ""), "organism": ORGANISM.get(d.name, "")}
        if d.name in DEFAULTS:
            entry["defaults"] = DEFAULTS[d.name]
        models.append(entry)
(data / "models.json").write_text(json.dumps(models, indent=1))
print("models:", [m["slug"] for m in models])
