# MetfluxTools

**[English](#english)** · **[Español](#español)**

---

## English

Offline preprocessing, testing and validation tools for **[MetfluxLive](https://github.com/JonatanFernandez/MetfluxLive)** ([metflux.live](https://metflux.live)).

**This repository is not part of the web application.** MetfluxLive is a static browser app that runs directly from the preprocessed JSON files shipped in its `data/` folder and needs neither Python, Node.js nor this repository at run time. `MetfluxTools` holds the scripts used to generate those data assets from public sources, verify the solver and tracer modules, and run the validation benchmarks reported in the paper.

### Setup

All scripts expect to live inside the `tools/` folder of a `MetfluxLive` checkout (they import the browser's ES modules from `../js/core/` and read/write `../data/`):

```bash
git clone https://github.com/JonatanFernandez/MetfluxLive.git
git clone https://github.com/JonatanFernandez/MetfluxTools.git MetfluxLive/tools
```

To run the Node.js tests and validation scripts (Node 18+), make sure `MetfluxLive/package.json` exists with:

```json
{ "type": "module" }
```

To rebuild models, atom maps or layouts, create the Python environment (CPU PyTorch is enough for RXNMapper):

```bash
cd MetfluxLive
python3 -m venv --without-pip tools/.venv
python3 -m pip --python tools/.venv/bin/python install rdkit pyyaml igraph rxnmapper \
    torch --extra-index-url https://download.pytorch.org/whl/cpu
```

Gibbs free energies (`build_thermo.py`) use a separate environment because of `equilibrator-api`'s dependencies:

```bash
python3 -m venv --without-pip tools/.venv-thermo
python3 -m pip --python tools/.venv-thermo/bin/python install equilibrator-api==0.8.1 rdkit
```

### Contents

```
DATA_SOURCES.md              provenance, URLs, licences and SHA-256 checksums of all external downloads
build_all.sh                 rebuilds data/<slug>/ for all eight GEMs
build_model.py               converts standard-GEM releases into model.json and structures.json
convert_foreign.py           maps yeast-GEM and E. coli iML1515 to Human-GEM metabolite names
fetch_mnx_structures.py      fetches openly licensed SMILES from MetaNetX/MNXref 4.5
build_atommaps.py            computes atom-to-atom maps (identity, cofactors, RXNMapper, MCS, overrides)
atommap_overrides.json       curated atom-mapping rules and manual corrections
build_layout.py              precomputes pathway and reaction layouts (layout.json)
build_thermo.py              estimates transformed Gibbs energies with eQuilibrator (thermo.json)
write_manifest.py            writes data/models.json
inspect_map.py, map_reaction.py   utilities to inspect maps or run RXNMapper on a single reaction SMILES
contexts/                    cell-type, tissue (HPA, Bgee, Cellosaurus) and microbial strain builders
test_*.mjs, diag_*.mjs       functional tests and diagnostics running the browser's js/core/ modules in Node
validation/                  validation scripts and generated outputs (Recon 3D, E. coli, liver, yeast, DepMap)
upstream/                    proposed upstream fix for yeast-GEM reaction r_0953 and verification scripts
```

Large third-party raw downloads (>100 MB, such as DepMap CSVs, Preciat et al. 2017 archive, HPA/Bgee expression dumps and build caches) and publisher supplementary files are excluded via `.gitignore`; see [`DATA_SOURCES.md`](DATA_SOURCES.md) for where to download each file and verify its checksum.

### Usage

Run commands from the `MetfluxLive` root directory:

```bash
# Rebuild preprocessed assets (requires source model archives in the MetfluxLive root; see DATA_SOURCES.md)
tools/build_all.sh

# Functional tests against MetfluxLive's js/core/ modules
node tools/test_fba.mjs           [human-gem|mouse-gem] [default|plasmax|hplm]
node tools/test_tracer.mjs        [human-gem|mouse-gem]
node tools/test_isotopologues.mjs [human-gem|mouse-gem]
node tools/test_dynamic.mjs       [human-gem|mouse-gem]
node tools/test_edits.mjs         [human-gem|mouse-gem]
node tools/test_contexts.mjs
node tools/test_objectives.mjs

# Validation benchmarks
node tools/validation/mids_vs_crown_ecoli.mjs 4000
node tools/validation/mids_vs_fujiwara_yeast.mjs
node tools/validation/mids_vs_grankvist.mjs tools/validation/grankvist_mids.json
node tools/validation/yeast_gem_tests.mjs
node tools/validation/depmap_evaluate.mjs
```

### Licence

Copyright © 2026 Jonatan Fernández-García.

The scripts written for MetfluxTools are free software, released under the GNU General Public License version 3 or (at your option) any later version (GPL-3.0-or-later). Derived or redistributed reference files under `validation/` (such as `validation/yeast/` from yeast-GEM and `validation/depmap/depmap26q*_metabolic.json` from DepMap) keep the licences of their original sources; see [`DATA_SOURCES.md`](DATA_SOURCES.md) for details.

---

## Español

Herramientas locales de preprocesado, pruebas y validación de **[MetfluxLive](https://github.com/JonatanFernandez/MetfluxLive)** ([metflux.live](https://metflux.live)).

**Este repositorio no forma parte de la aplicación web.** MetfluxLive es una aplicación estática que se ejecuta íntegramente en el navegador a partir de los ficheros JSON ya preprocesados de su carpeta `data/` y no necesita este repositorio, Python ni Node.js en tiempo de ejecución. `MetfluxTools` reúne el código usado para generar esos ficheros a partir de fuentes públicas, probar los módulos del solver y del trazador, y ejecutar las validaciones del artículo.

### Instalación

Todos los scripts están pensados para vivir en la carpeta `tools/` dentro de una copia de `MetfluxLive` (importan los módulos ES del navegador desde `../js/core/` y leen/escriben en `../data/`):

```bash
git clone https://github.com/JonatanFernandez/MetfluxLive.git
git clone https://github.com/JonatanFernandez/MetfluxTools.git MetfluxLive/tools
```

Para ejecutar las pruebas y validaciones con Node.js (Node 18+), asegúrate de que exista `MetfluxLive/package.json` con:

```json
{ "type": "module" }
```

Para regenerar modelos, mapas atómicos o disposiciones, crea el entorno de Python (PyTorch en CPU basta para RXNMapper):

```bash
cd MetfluxLive
python3 -m venv --without-pip tools/.venv
python3 -m pip --python tools/.venv/bin/python install rdkit pyyaml igraph rxnmapper \
    torch --extra-index-url https://download.pytorch.org/whl/cpu
```

Las energías libres de Gibbs (`build_thermo.py`) usan un entorno aparte por las dependencias de `equilibrator-api`:

```bash
python3 -m venv --without-pip tools/.venv-thermo
python3 -m pip --python tools/.venv-thermo/bin/python install equilibrator-api==0.8.1 rdkit
```

### Contenido

```
DATA_SOURCES.md              procedencia, URLs, licencias y sumas SHA-256 de todas las descargas externas
build_all.sh                 regenera data/<slug>/ para los ocho GEM
build_model.py               convierte publicaciones standard-GEM en model.json y structures.json
convert_foreign.py           empareja yeast-GEM e iML1515 (E. coli) con nombres de metabolitos de Human-GEM
fetch_mnx_structures.py      obtiene SMILES con licencia abierta desde MetaNetX/MNXref 4.5
build_atommaps.py            calcula los mapas átomo a átomo (identidad, cofactores, RXNMapper, MCS, overrides)
atommap_overrides.json       reglas curadas de mapeo atómico y correcciones manuales
build_layout.py              precalcula disposiciones de rutas y reacciones (layout.json)
build_thermo.py              estima energías de Gibbs transformadas con eQuilibrator (thermo.json)
write_manifest.py            genera data/models.json
inspect_map.py, map_reaction.py   utilidades para inspeccionar mapas o ejecutar RXNMapper en un SMILES de reacción
contexts/                    generación de contextos de tipos celulares/tejidos (HPA, Bgee, Cellosaurus) y cepas
test_*.mjs, diag_*.mjs       pruebas funcionales y diagnósticos que ejecutan los módulos js/core/ en Node
validation/                  scripts de validación y resultados generados (Recon 3D, E. coli, hígado, levadura, DepMap)
upstream/                    propuesta de corrección para la reacción r_0953 de yeast-GEM y scripts de verificación
```

Las descargas externas grandes (>100 MB, como los CSV de DepMap, el archivo de Preciat et al. 2017, los ficheros de expresión de HPA/Bgee y las cachés) y los suplementos originales de editoriales están excluidos en `.gitignore`; consulta [`DATA_SOURCES.md`](DATA_SOURCES.md) para ver de dónde descargar cada fichero y comprobar su SHA-256.

### Uso

Ejecuta los comandos desde la raíz de `MetfluxLive`:

```bash
# Regenerar los ficheros preprocesados (requiere los archivos fuente de los modelos en la raíz; ver DATA_SOURCES.md)
tools/build_all.sh

# Pruebas funcionales sobre los módulos js/core/ de MetfluxLive
node tools/test_fba.mjs           [human-gem|mouse-gem] [default|plasmax|hplm]
node tools/test_tracer.mjs        [human-gem|mouse-gem]
node tools/test_isotopologues.mjs [human-gem|mouse-gem]
node tools/test_dynamic.mjs       [human-gem|mouse-gem]
node tools/test_edits.mjs         [human-gem|mouse-gem]
node tools/test_contexts.mjs
node tools/test_objectives.mjs

# Validaciones
node tools/validation/mids_vs_crown_ecoli.mjs 4000
node tools/validation/mids_vs_fujiwara_yeast.mjs
node tools/validation/mids_vs_grankvist.mjs tools/validation/grankvist_mids.json
node tools/validation/yeast_gem_tests.mjs
node tools/validation/depmap_evaluate.mjs
```

### Licencia

Copyright © 2026 Jonatan Fernández-García.

El código de MetfluxTools es software libre, publicado bajo la Licencia Pública General de GNU versión 3 o, a tu elección, cualquier versión posterior (GPL-3.0-or-later). Los ficheros de referencia derivados o copiados bajo `validation/` (como `validation/yeast/` de yeast-GEM y `validation/depmap/depmap26q*_metabolic.json` de DepMap) conservan las licencias de sus fuentes originales; ver [`DATA_SOURCES.md`](DATA_SOURCES.md).
