#!/usr/bin/env bash
# Rebuild every data/<slug>/ folder from the model archives in the project root.
#   tools/build_all.sh            -> Human, Mouse, Rat, Zebrafish, Fruitfly, Worm-GEM, yeast-GEM and E. coli iML1515
# Requires the Python environment in tools/.venv (see README).
set -euo pipefail
cd "$(dirname "$0")/.."
PY=tools/.venv/bin/python

build() {  # build <archive-or-folder> <slug> [extra build_model args]
  local src=$1 slug=$2; shift 2
  $PY tools/build_model.py "$src" "$slug" "$@"
  $PY tools/build_atommaps.py "$slug" 2> >(grep -v Warning >&2)
  $PY tools/build_layout.py "$slug" | tail -1
}

DEMAND="--add-demand MFL_DM_protein:protein_pool_biomass:c"   # protein secretion objective (plasma cells)
# growth without the cofactor pool, for gene essentiality (docs, "Validation": DepMap)
DEMAND="$DEMAND --add-biomass-variant MFL_biomass_nocof:MAR13082:cofactor_pool_biomass"
build "$(ls -d Human-GEM-*.zip | tail -1)" human-gem $DEMAND
build "$(ls -d Mouse-GEM-*.zip | tail -1)" mouse-gem --structures-from human-gem $DEMAND
build "$(ls -d Rat-GEM-*.zip | tail -1)" rat-gem --structures-from human-gem $DEMAND
build "$(ls -d Zebrafish-GEM-*.zip | tail -1)" zebrafish-gem --structures-from human-gem $DEMAND
build "$(ls -d Fruitfly-GEM-*.zip | tail -1)" fruitfly-gem --structures-from human-gem $DEMAND
# C. elegans has no bisphosphoglycerate mutase: 2,3-BPG leaves the biomass metabolite pool (docs, "Models")
build "$(ls -d Worm-GEM-*.zip | tail -1)" worm-gem --structures-from human-gem $DEMAND --drop-from "MAR10064:2,3-bisphospho-D-glycerate"
# models with other identifiers: converted first, metabolites matched to Human-GEM names (tools/convert_foreign.py)
# structures of the unmatched metabolites from MetaNetX (needs a first conversion to know which; streamed, ~1.5 GB)
if [ ! -f tools/cache/mnx_structures.json ]; then
  $PY tools/convert_foreign.py "$(ls -d yeast-GEM-*.zip | tail -1)" yeast-gem
  $PY tools/convert_foreign.py ecoli-src/iML1515.json ecoli-gem
  $PY tools/fetch_mnx_structures.py
fi
$PY tools/convert_foreign.py "$(ls -d yeast-GEM-*.zip | tail -1)" yeast-gem
build tools/cache/src/yeast-gem yeast-gem --structures-from human-gem
$PY tools/convert_foreign.py ecoli-src/iML1515.json ecoli-gem   # http://bigg.ucsd.edu/static/models/iML1515.json
build tools/cache/src/ecoli-gem ecoli-gem --structures-from human-gem
node tools/contexts/build_strains.mjs   # laboratory strains of yeast and E. coli
# free energies of reaction (optional: needs tools/.venv-thermo, see README)
if [ -x tools/.venv-thermo/bin/python ]; then
  for slug in human-gem mouse-gem rat-gem zebrafish-gem fruitfly-gem worm-gem yeast-gem ecoli-gem; do tools/.venv-thermo/bin/python -W ignore tools/build_thermo.py "$slug" 2>&1 | tr '\r' '\n' | grep -v '%|' | tail -2; done
fi
$PY tools/write_manifest.py
