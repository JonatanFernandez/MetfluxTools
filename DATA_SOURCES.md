# Data sources of the tools

Every file that was downloaded to build or validate MetfluxLive, where it came from, which release it is and which
script uses it. Everything else under `tools/` is code, or a file computed from these (listed at the end, so that
nobody mistakes it for a download).

How to read it:

- **Stored** files are in this folder tree (or at the project root, for the model archives) and have a SHA-256 below.
- **Streamed** files were read from the web by a script and filtered; they are not kept.
- **Not kept** files were read once to take numbers or tables out of them; only the result is in the repository.
- Download dates are the dates of the files (all in November 2026), taken from the file system; the sessions that
  made the downloads are not otherwise dated. Where a URL could not be recovered it says so instead of guessing.
- Paths are relative to the project root. Licences are those declared by each source; check them before you
  redistribute a file. Files marked *by hand* were saved from a browser, not by a script.

## 1. Metabolic models (inputs of `tools/build_all.sh`)

| Model | Release | Source | Stored as |
|---|---|---|---|
| Human-GEM | 2.1.1 | https://github.com/SysBioChalmers/Human-GEM (release archive) | `Human-GEM-2.1.1.zip` |
| Mouse-GEM | 1.8.0 | https://github.com/SysBioChalmers/Mouse-GEM | `Mouse-GEM-1.8.0.zip` |
| Rat-GEM | 1.8.0 | https://github.com/SysBioChalmers/Rat-GEM/archive/refs/tags/v1.8.0.zip | `Rat-GEM-1.8.0.zip` |
| Zebrafish-GEM | 1.6.0 | https://github.com/SysBioChalmers/Zebrafish-GEM | `Zebrafish-GEM-1.6.0.zip` |
| Fruitfly-GEM | 1.7.0 | https://github.com/SysBioChalmers/Fruitfly-GEM | `Fruitfly-GEM-1.7.0.zip` |
| Worm-GEM | 1.8.0 | https://github.com/SysBioChalmers/Worm-GEM | `Worm-GEM-1.8.0.zip` |
| yeast-GEM | 9.1.1 | https://github.com/SysBioChalmers/yeast-GEM/archive/refs/tags/v9.1.1.zip | `yeast-GEM-9.1.1.zip` |
| iML1515 (*E. coli*) | BiGG | http://bigg.ucsd.edu/static/models/iML1515.json | `ecoli-src/iML1515.json` |

The SysBioChalmers archives follow the pattern `https://github.com/SysBioChalmers/<Name>/archive/refs/tags/v<release>.zip`.
The URL of the Human-GEM and Mouse-GEM archives (saved on 12 November, at the start of the project) was not logged;
the pattern above is the one used for the others. Mouse-GEM's `data/MGI_gene_ID_mapping.tsv`, which
`tools/contexts/build_expression.py` needs, is inside its archive. Papers: Robinson et al., *Sci. Signal.* 2020
(Human-GEM); Wang et al., *PNAS* 2021 (the animal models); Lu et al., *Nat. Commun.* 2019 (yeast-GEM); Monk et al.,
*Nat. Biotechnol.* 2017 (iML1515).

| File | Size | SHA-256 |
|---|---|---|
| `Human-GEM-2.1.1.zip` | 9.2 MB | `62a99070bc017a0b46cdb34a1d1cb9fc3391be9f3171cf0756ebd9c7f789ee9f` |
| `Mouse-GEM-1.8.0.zip` | 7.1 MB | `30a2bc55fcd142349279fdf0e615041a7bdb566b5ad0af12c1748689ab07daf8` |
| `Rat-GEM-1.8.0.zip` | 4.5 MB | `d385530c0a28e012f30e2e30797196dc7d27b5be9a8fbf8bcddb1606b1b4e6d5` |
| `Zebrafish-GEM-1.6.0.zip` | 3.3 MB | `c782a54769007d3ffe63a12ec70ebc85c25f7c2789432ae9ba62f34ff9518f3a` |
| `Fruitfly-GEM-1.7.0.zip` | 2.9 MB | `f5e0f016c5c2c4f690251afb2a2fc58e11959437360cfdf71fbd8061cae8d05e` |
| `Worm-GEM-1.8.0.zip` | 2.9 MB | `f2c39e7bd9053cd258a33d497193d7d227afc0aa3606a4c74f7516f02702f037` |
| `yeast-GEM-9.1.1.zip` | 8.5 MB | `f9cfb4e05f578d48b97919df498f32dd8990f882a554217cebfc66d522d21289` |
| `ecoli-src/iML1515.json` | 3.1 MB | `b0f9199f048779bb08a14dfa6c09ec56d35b8750d2f99681980d0f098355fbf5` |

`tools/convert_foreign.py` turns the yeast-GEM archive and `iML1515.json` into the folders `tools/cache/src/yeast-gem/`
and `tools/cache/src/ecoli-gem/` (computed, not downloaded).

## 2. Gene expression and cell lines (`tools/contexts/`, `tools/cache/expression/`)

| File | Source | Release | Used by |
|---|---|---|---|
| `rna_celline.tsv.zip` | https://www.proteinatlas.org/download/tsv/rna_celline.tsv.zip | Human Protein Atlas 26.0 | `contexts/build_expression.py`, `contexts/cellline_meta.py` |
| `rna_single_cell_type.tsv.zip` | https://www.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip | HPA 26.0 | `contexts/build_expression.py` |
| `rna_tissue_consensus.tsv.zip` | https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip | HPA 26.0 | `contexts/build_expression.py` |
| `cellosaurus.txt` | https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt | Cellosaurus 57.0 (last update 22 October 2026) | `contexts/cellline_meta.py` (cell line categories, diseases, DepMap and NCI-60 cross-references) |
| `Mus_musculus_expr_simple.tsv.gz` | https://www.bgee.org/ftp/bgee_v15_2/download/calls/expr_calls/Mus_musculus_expr_simple.tsv.gz | Bgee 15.2 | `contexts/build_expression.py` (Mouse-GEM) |
| `Rattus_norvegicus_expr_simple.tsv.gz` | same folder, `Rattus_norvegicus_expr_simple.tsv.gz` | Bgee 15.2 | Rat-GEM |
| `Danio_rerio_expr_simple.tsv.gz` | same folder | Bgee 15.2 | Zebrafish-GEM |
| `Drosophila_melanogaster_expr_simple.tsv.gz` | same folder | Bgee 15.2 | Fruitfly-GEM |
| `Caenorhabditis_elegans_expr_simple.tsv.gz` | same folder | Bgee 15.2 | Worm-GEM |

All in `tools/cache/expression/`. The first Mouse download came from the `current` Bgee folder
(`https://www.bgee.org/ftp/current/...`); the other species were fetched with the explicit `bgee_v15_2` path. Licences:
HPA CC BY-SA 4.0, Cellosaurus CC BY 4.0; Bgee: see https://www.bgee.org/about.

| File | Size | SHA-256 |
|---|---|---|
| `tools/cache/expression/rna_celline.tsv.zip` | 205.9 MB | `e5b11779308c2720cfe89d9cb19932ad4c0c3fea0fa37e8c61a14d501995852a` |
| `tools/cache/expression/rna_single_cell_type.tsv.zip` | 16.4 MB | `2fdc5094f0182400e6a611ddf66b1251a4c703e45ae0ef72196fae62452160de` |
| `tools/cache/expression/rna_tissue_consensus.tsv.zip` | 5.3 MB | `1238e89fb1094c031db0b44cba9ed97e69f2465bdcc9065d13771fd05414108b` |
| `tools/cache/expression/cellosaurus.txt` | 122.4 MB | `3b5a9fc56b49854230ae90fb7cee57a3bf9f961df98be215a634203526f476ef` |
| `tools/cache/expression/Mus_musculus_expr_simple.tsv.gz` | 160.6 MB | `a8a001cf604106f4b0b98228487508a0c98a9bc1c22745b957d9a51ae3976890` |
| `tools/cache/expression/Rattus_norvegicus_expr_simple.tsv.gz` | 9.4 MB | `393262828d5937f3d71bfbf6db260dd7e0fa809a6d783e8c605f17dd04e3e069` |
| `tools/cache/expression/Danio_rerio_expr_simple.tsv.gz` | 15.8 MB | `5136bdacdcfce7c0e9f33b28556668152488522b2247c4de9435356a2482c025` |
| `tools/cache/expression/Drosophila_melanogaster_expr_simple.tsv.gz` | 43.6 MB | `d6e8de57a6405a0fdad1d878386029ba922f32448f58b2f669573b0cce068066` |
| `tools/cache/expression/Caenorhabditis_elegans_expr_simple.tsv.gz` | 2.6 MB | `c59b4519943644e45718e391e8a270fa49495d99c84f792de2f8e0ca70e185c5` |

Computed from them (not downloads): `cache/expression/cellline_meta.json`, `cache/contexts/*` and, in the app, `data/<slug>/contexts.json`.

## 3. Metabolite structures (`tools/fetch_mnx_structures.py`)

**Streamed**, not kept: `chem_prop.tsv`, `chem_xref.tsv` and `chem_depr.tsv` of MetaNetX/MNXref 4.5, from
https://www.metanetx.org/ftp/latest/ (CC BY 4.0; Moretti et al., *Nucleic Acids Res.* 2021). Only structures whose
reference is ChEBI, SwissLipids or LIPID MAPS are kept. Result, computed: `tools/cache/mnx_structures.json`.
`latest` moves: re-running the script later may use a newer MNXref than 4.5.

## 4. Thermodynamics (`tools/build_thermo.py`)

The eQuilibrator compound database (Zenodo, https://doi.org/10.5281/zenodo.4128543) is downloaded by `equilibrator-api`
the first time `build_thermo.py` runs, into that user's cache (outside this repository). Midpoint potentials of the
cytochrome c and ubiquinone couples are literature values written in `COUPLES` inside the script, with their references.
Packages: `tools/.venv-thermo` (see section 8).

## 5. Atom maps (`tools/build_atommaps.py`, `tools/validation/atommaps_vs_recon3d.py`)

| File | Source | Stored as |
|---|---|---|
| Additional file 3 of Preciat Gonzalez et al., *J. Cheminform.* 2017, 9:39 (curated Recon 3D mappings and the published algorithms' predictions) | https://pmc.ncbi.nlm.nih.gov/articles/instance/5471289/bin/13321_2017_223_MOESM3_ESM.zip (doi 10.1186/s13321-017-0223-1) | `tools/validation/preciat2017_additional_file3.zip` |

| File | Size | SHA-256 |
|---|---|---|
| `tools/validation/preciat2017_additional_file3.zip` | 108.4 MB | `8acf91706f8f77da8d471f162f7e7afe0c48973168e2b98c1000c52d7d0e80e7` |

RXNMapper (Schwaller et al., *Sci. Adv.* 2021, https://github.com/rxn4chemistry/rxnmapper) is installed with `pip` in
`tools/.venv`; its pretrained model comes with the package. Computed or curated here, not downloaded:
`tools/atommap_overrides.json`, `tools/validation/curated_rules.json`, `tools/cache/rxnmapper_cache.json`,
`tools/cache/mcs_cache.json`, `tools/cache/<slug>/canon.json` and `smiles.json`.

## 6. Validation data

### Labelling experiments

| File | Source | Used by |
|---|---|---|
| `tools/validation/crown2015_mmc1.xlsx`, `_mmc2.docx`, `_mmc3.docx`, `_mmc4.xlsx`, `_mmc5.xlsx` | Supplementary files 1–5 of Crown, Long & Antoniewicz, *Metab. Eng.* 2015, 28:151 (doi 10.1016/j.ymben.2015.01.001; PMC5845449). Pattern `https://ars.els-cdn.com/content/image/1-s2.0-S1096717615000038-mmc<n>.<ext>` | `validation/mids_vs_crown_ecoli.mjs` (through `crown2015.json`) |
| `tools/validation/fujiwara2025_mmc1.xlsx` | Supplementary Data S1 of Fujiwara et al., *Metab. Eng. Commun.* 2025, 20:e00260 (doi 10.1016/j.mec.2025.e00260; PMC12008597, via Europe PMC supplementary files). CC BY-NC 4.0 | `validation/mids_vs_fujiwara_yeast.mjs` (through `fujiwara2025_S1.json`) |
| Supplementary Data 3, sheet "MIDs", of Grankvist et al., *Nat. Metab.* 2024, 6:1963 (doi 10.1038/s42255-024-01119-3) | `https://static-content.springer.com/esm/art%3A10.1038%2Fs42255-024-01119-3/MediaObjects/42255_2024_1119_MOESM<n>_ESM.xlsx` (the number <n> of the file was not recorded) | **Not kept**: only `tools/validation/grankvist_mids.json`, extracted from it; `validation/mids_vs_grankvist.mjs` |

| File | Size | SHA-256 |
|---|---|---|
| `tools/validation/crown2015_mmc1.xlsx` | 95 kB | `ee7531e261bcc247b0e48e4e09f4201d211bc677b12d7968d766504ec80d5058` |
| `tools/validation/crown2015_mmc2.docx` | 21 kB | `27b0fa49b40cb5162ae437aacd796d9534c7017ea91340e9b2b87f42f48dee1e` |
| `tools/validation/crown2015_mmc3.docx` | 25 kB | `031eb31390ba1af09b0d4fe5db3b396933ebcd39bdf56a281a7272f924e29b7d` |
| `tools/validation/crown2015_mmc4.xlsx` | 21 kB | `7941616ae0d77cf937f78a177912cd606eb1c2f182f175b52cb3d7cace9af503` |
| `tools/validation/crown2015_mmc5.xlsx` | 27 kB | `d6fd75bc0c0e60343ecbe93a6aa83a3d57ffe4694026536ccdc22fe8e164fa92` |
| `tools/validation/fujiwara2025_mmc1.xlsx` | 57 kB | `ed2643f5ae93ab85020033d49f57b5bda9c95c8034fc059bafdce98a25afdbf3` |

### Intracellular metabolomics (HEK293 producer vs non-producer)

| File | Source | Used by |
|---|---|---|
| `tools/validation/dietmair2012_intracellular.json` | Table S3 of Dietmair et al., *PLoS ONE* 7:e43394 (2012), doi 10.1371/journal.pone.0043394, file `pone.0043394.s003.xlsx`, CC BY 4.0. The whole supplement is one zip at `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3427347/supplementaryFiles`; only this sheet is used | `validation/metabolomics_dietmair.mjs`, through `prepare_dietmair.py` |

52 metabolites of central metabolism, relative intracellular concentrations, two cell lines x 4 time points,
mean ± standard error. The xlsx itself is **not kept**: `prepare_dietmair.py` turns it into the JSON above.
The same paper's Table 2 is the `hek293` measured exchange set in `js/core/measured.js`, so the exometabolome
and the intracellular metabolome of one experiment are both in the repository.

### Metabolomics studies for the worked examples

`tools/validation/fetch_metabolomics.sh` downloads five studies from the
[Metabolomics Workbench](https://www.metabolomicsworkbench.org) REST API (NIH Common Fund) into
`tools/validation/metabolomics/`, as mwTab text files that carry the metabolite table, the sample factors
and the units together. `tools/validation/metabolomics_catalogue.json` describes each one and why it was
chosen. The files are **not kept** in the repository, and each study has its own DOI and terms: check them
before redistributing any of the data.

| Study | Analysis | Organism | Metabolites | Samples | Units | Design |
|---|---|---|---|---|---|---|
| ST002076 | AN003383 | *E. coli* K-12 MG1655 | 22 | 24 | **µmol/gDCW** | wild type and three regulator knockouts, anaerobic |
| ST004848 | AN008205 | *M. musculus* MPC11 | 118 | 90 | AU | control vs methionine restriction, 3 time points |
| ST004339 | AN007246 | *R. norvegicus* INS-1 | 53 | 19 | abundance | buthionine sulfoximine (inhibits GCLC) vs control |
| ST001609 | AN002643 | *S. cerevisiae* BY4742 | 1 | 26 | **mM** | wild type vs *mtl1Δ*, with OD600 per culture |
| ST003960 | AN006508 | *H. sapiens* HEK293 | 2 | 12 | ratio | spent medium: fresh, WT, WT + antimycin A, SDHAF2 + antimycin A |

The REST path that works is `/rest/study/analysis_id/<AN>/mwtab/txt`; the `/datatable` endpoint named in the
API documentation answers 200 with "This function does not exist".

### Yeast benchmarks

`tools/validation/yeast/` holds the files that yeast-GEM 9.1.1 uses for its own release tests, copied unchanged from
`yeast-GEM-9.1.1.zip` (section 1): `data/essentialGenes/{inviable_orfs,verified_orfs}.txt`,
`data/testResults/essentialGenes.tsv`, `data/conditions/anaerobic.yml`,
`data/physiology/{Biolog_Substrate, aminoAcid_Bjorkeroth2020, biomassComposition_*, chemostatData_Tobias2013,
exchange_data_anaerobic, flux_data_anaerobic, chemostatData_VanHoek1998}.tsv`, plus the repository's `LICENSE.md` and the
`README.md` of its `data/physiology` folder (CC BY 4.0). Used by `validation/yeast_gem_tests.mjs`.

`tools/upstream/yeast-GEM-r_0953/` is a proposed correction to yeast-GEM (reaction r_0953), written during this work. Its
`files/` were taken from the yeast-GEM repository (https://github.com/SysBioChalmers/yeast-GEM, `develop` branch, 19 November)
and `verification/*.py` are scripts of ours.

### DepMap (gene essentiality; `validation/depmap_prepare.py`, `depmap_essentiality.mjs`, `depmap_evaluate.mjs`)

DepMap Public **26Q3**, downloaded *by hand* on 21 November from https://depmap.org/portal/data_page/ (Broad Institute;
CC BY 4.0). The portal's download API did not give a stable URL, so the files are identified by name and checksum:

| File | Content |
|---|---|
| `Model.csv` | cell line metadata (`ModelID`, `RRID` = Cellosaurus accession, lineage, growth medium `FormulationID`) |
| `CRISPRGeneEffect.csv` | Chronos gene effect |
| `CRISPRGeneDependency.csv` | probability of dependency |
| `CRISPRInferredCommonEssentials.csv` | genes essential in nearly all lines |
| `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | log2(TPM+1) expression; downloaded for future cell-type contexts, **not used yet** |

| File | Size | SHA-256 |
|---|---|---|
| `tools/validation/depmap/Model.csv` | 0.8 MB | `65eb1895800e71a963a3ab56e5d08c8d006cb0d9e7aa68f7e34a3810a2a3cbc7` |
| `tools/validation/depmap/CRISPRGeneEffect.csv` | 497.1 MB | `1e0c492944a7297c0eda179081a1e81a1cfe41a97628facd7eddbdda1a7504ee` |
| `tools/validation/depmap/CRISPRGeneDependency.csv` | 488.2 MB | `c1177633cba85e460d8c4dad8dd1388b5744629e1e267bc354cb7dfc1794a718` |
| `tools/validation/depmap/CRISPRInferredCommonEssentials.csv` | 31 kB | `149e173a83545f296a933ecaab337c2842b59900d4dc32e84f08569112b43613` |
| `tools/validation/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | 395.0 MB | `252d1fedc7014953b19ffe5d68151b4df973c8078b12c9e5ccd97e292548d177` |

Computed: `tools/validation/depmap/depmap26q3_metabolic.json` (`depmap_prepare.py`), the predictions `pred_*`, `evaluation*.json`
and the logs.

The first screen, on 122 lines, used **26Q1** instead: `gene_effect.csv` from the DepMap Public 26Q1 release on figshare
(article 31660582, https://ndownloader.figshare.com/files/67214582), with the lines matched to the catalogue through the
`DR DepMap` cross-references of Cellosaurus. That file was not kept, and the one-off script that made
`tools/validation/depmap/depmap26q1_metabolic.json` was not saved (its first lines are in the session log). The 26Q1 file
is only needed to reproduce the 122-line subset, not the 730-line run.

## 6b. Signed regulatory networks (`tools/build_regulons.mjs`)

Downloaded on **24 November** (rat on **25 November**) and cached under `tools/cache/regulons/`. The build script keeps the cache, so a rerun does
not download again; delete a cached file to refresh it.

**CollecTRI for human, mouse and rat**, through the OmniPath web service (Muller-Dott *et al.*, *NAR* 2023; Turei *et al.*).
CollecTRI is **CC BY-SA 4.0**, which is why the human, mouse and rat regulon files are distributed under that licence. The mouse and rat
networks are CollecTRI mapped by orthology, not curated in those species; rat coverage is thinner (425 of 2,804
metabolic genes). CollecTRI serves no further organism: the API rejects 7227 (fly), 6239 (worm) and 7955 (zebrafish).

| File | Query | SHA-256 |
|---|---|---|
| `collectri_9606.tsv` | `https://omnipathdb.org/interactions?datasets=collectri&organisms=9606&genesymbols=yes&fields=sources` | `67e89b12ab66b59acbcd4dc3c10d7fb7045a2fce17f85d3a8fceaa51dab984ec` |
| `collectri_10090.tsv` | the same with `organisms=10090` | `c7873e461a92da93c9946d5539b83b40dd88e0e1cbd00321ec939383d6719932` |
| `collectri_10116.tsv` | the same with `organisms=10116` (rat, by orthology) | `2bfeadac4f3719620c4aa63e01d9c135aded129dcef6972f33744005dd9f1045` |

**RegulonDB v14.5.1 for *E. coli*** (CCG-UNAM, Salgado *et al.* 2024), through the GraphQL endpoint at
`https://regulondb.ccg.unam.mx/graphql`, dataset `NetworkRegulatorGene`. Its plain download URLs no longer resolve; the
GraphQL service is the working route. **Its licence is academic / non-commercial**, which is narrower than everything else
here, so the *E. coli* regulon file cannot be redistributed commercially. This is stated in the paper.

| File | SHA-256 |
|---|---|
| `regulondb_NetworkRegulatorGene.tsv` | `99b2477bcc8cb11e8d5b9a7fd9e16f33d93311c5dbdabaa07d304652259e3483` |

The build keeps only edges whose target is a metabolic gene of the corresponding model and whose sign is unambiguous;
edges annotated as both activating and repressing are discarded. The result is `data/<slug>/regulons.json`
(human 917 TFs / 1,327 metabolic genes; mouse 838 / 1,069; rat 837 / 425; *E. coli* 239 / 1,129), used by `js/core/reginfer.js`.

## 6c. The paired multi-omic worked example (`data/demos/`)

The only experiment we found where a whole transcriptome and an intracellular metabolome of the **same cultures**
are both public: anaerobic double mutants of *E. coli* K-12 MG1655 from the Systems Biology and Metabolic
Engineering Laboratory, IIT Bombay. Four strains (WT, Δ*fnr*Δ*arcA*, Δ*arcA*Δ*ihf*, Δ*fnr*Δ*ihf*) grown together
in glucose minimal medium under anaerobic fermentation, with matching sample names in both repositories.

| Source | What | Where |
|---|---|---|
| **GEO GSE195954** | RNA-seq raw counts, 4,466 genes, 2 replicates per strain | `tools/cache/GSE195954_Raw_counts_msi.txt.gz` |
| **Metabolomics Workbench ST002076** | 21 intracellular metabolites, 6 samples per strain | `tools/validation/metabolomics/ST002076.mwtab` |

`tools/` is not shipped, so the two tables the worked example loads are derived into `data/demos/` by
`scratchpad` generation (CPM + moderated *t* for the transcriptome, as every other RNA-seq benchmark here;
the metabolite matrix is the mwtab block for the two conditions with a `group` row added):

| File | SHA-256 |
|---|---|
| `data/demos/ecoli_arcA_ihf_rnaseq.tsv` | `be64868409299246587b695ed4c2361b750118f98f50d8f5f60f730b844989bf` |
| `data/demos/ecoli_arcA_ihf_metabolomics.tsv` | `9ff851213f273cbc55cdea2a95184ad028a3c3a914fba576deb00e6f3c196de7` |

The design carries its own specificity control: each strain lacks two of the three global anaerobic regulators,
so the expected answer differs per contrast. Measured with `inferRegulatorsMultiOmic`, the two deleted factors
outrank the intact one every time -- Δ*arcA*Δ*ihf*: ArcA 3rd, IHF 8th, FNR 130th; Δ*fnr*Δ*arcA*: ArcA 3rd,
FNR 4th, IHF 150th; Δ*fnr*Δ*ihf*: IHF 2nd, FNR 33rd, ArcA 44th.

## 7. Values taken from papers into the application (the files are not kept)

These were read once; the numbers are in `js/core/measured.js`, `js/core/medium.js` and `docs.html`, each with its
reference. Listed so that the downloads are traceable.

| Used for | Source | Where it was read |
|---|---|---|
| NCI-60 exchange fluxes | Jain et al., *Science* 2012 (doi 10.1126/science.1218595); per-line values from S1 Table Q of Aurich et al., *PLoS Comput. Biol.* 2017 (doi 10.1371/journal.pcbi.1005698) | Database S1: https://pmc.ncbi.nlm.nih.gov/articles/instance/3526189/bin/NIHMS422185-supplement-Database_S1.xls |
| Plasmax and HPLM | Table S1 of Vande Voorde et al., *Sci. Adv.* 2019 (doi 10.1126/sciadv.aau7314) | https://pmc.ncbi.nlm.nih.gov/articles/instance/6314821/bin/aau7314_SM.pdf |
| Fibroblast fluxes | Lemons et al., *PLoS Biol.* 2010 (doi 10.1371/journal.pbio.1000514) | supplementary files, PMC2958657 |
| HEK293 fluxes | Dietmair et al., *PLoS ONE* 2012 | the article and its supplement (URL not recorded) |
| Liver slices, Williams E | Supplementary Table 2 and Data of Grankvist et al., *Nat. Metab.* 2024 | static-content.springer.com, as above |
| Yeast batch exchange rates (`yeastSD` measured set) | Fujiwara et al. 2025 (section 6) | PMC12008597 |
| RPMI 1640, DMEM, Schneider's | Gibco formulation tables 11875, 11965 and 21720 (Thermo Fisher) | the manufacturer's product pages (URLs not recorded) |

## 8. Software environments (not data)

- `tools/.venv`: model conversion, atom maps (RXNMapper with PyTorch; the CPU wheels come from
  https://download.pytorch.org/whl/cpu), RDKit.
- `tools/.venv-thermo`: `equilibrator-api` and its dependencies.
- `vendor/highs/`: HiGHS compiled to WebAssembly (https://highs.dev), part of the application.

No requirements file has been written for them yet; `pip freeze` in each environment gives the exact versions.

## 9. Computed files that are not downloads

`tools/cache/{mcs_cache,rxnmapper_cache,mnx_structures}.json`, `tools/cache/<slug>/*.json`, `tools/cache/contexts/*`,
`tools/cache/src/*`, `tools/cache/*.log`, `tools/atommap_overrides.json`, `tools/validation/*.json` made by the scripts next to
them, `tools/validation/*.tsv` and everything in `tools/validation/depmap/` except the five DepMap files above. They can
be regenerated from the files in sections 1–6, except the curated ones (`atommap_overrides.json`, `curated_rules.json`),
which hold manual decisions.
