# Contribución a yeast-GEM: reacción r_0953 (piridoxal oxidasa)

Este documento tiene tres partes:

1. **Para ti (en español):** qué se propone, cómo se ha verificado y los pasos exactos para enviarlo.
2. **Textos para pegar en GitHub (en inglés):** la issue (A) y el pull request (B).
3. **Apéndice (en inglés):** los cálculos completos y cómo reproducirlos. La issue enlaza o resume este apéndice.

Todo lo necesario está en esta carpeta (`tools/upstream/yeast-GEM-r_0953/`):

| Fichero | Qué es |
|---|---|
| `r_0953.patch` | El cambio completo, listo para `git apply` sobre `develop` |
| `files/` | Las versiones modificadas de cada fichero, por si prefieres copiarlas a mano |
| `verification/` | Los scripts con los que se ha comprobado todo |

---

## 1. Resumen para ti

### Por qué yeast-GEM y no Human-GEM

En Human-GEM ya estaba todo corregido (PR #1132) o abierto (#1128, #1129). En yeast-GEM, en cambio, la issue #428, abierta por un mantenedor, lista 22 reacciones con "alarma" termodinámica y dice explícitamente que **ninguna se ha revisado todavía**. r_0953 es la más extrema de las 22 (+375 kJ/mol). Revisar una de ellas a fondo es una aportación útil y no duplica trabajo de nadie.

### Qué está mal en r_0953

La reacción está escrita así:

`piridoxal + NH4+ + 2 H2O + 0,5 O2 → piridoxamina + 2 H2O2`  (gen PDX3/YBR035C, irreversible)

- **Está equilibrada en átomos y en carga.** Lo he comprobado: no es un error de fórmula.
- **El problema es químico.** Reduce piridoxal a piridoxamina usando como donante de electrones el agua, que se oxida a agua oxigenada. Es como fabricar H2O2 a partir de agua.
- **Termodinámicamente imposible:** ΔG'm = +375 ± 7 kJ/mol. Para compensarlo harían falta concentraciones unas 10⁶⁶ veces más favorables que 1 mM.
- **Origen:** viene de la reacción `PYDXO` de BiGG (modelos iND750 e iMM904), que pasó también a ModelSEED (`rxn09566`) y MetaNetX (`MNXR103365`). No existe en KEGG ni en Rhea, y la nomenclatura IUBMB no la recoge.

### Qué se propone: eliminarla

Se había pensado corregirla a su forma química razonable (piridoxamina + O2 + H2O → piridoxal + NH4+ + H2O2), pero no hay pruebas de que Pdx3 actúe sobre piridoxamina **no fosforilada**:

- EC 1.4.3.5 (IUBMB), Rhea y UniProt (P38075) solo recogen los sustratos fosforilados (PMP, PNP).
- La enzima equivalente de *E. coli* (PdxH) "no aceptó sustratos no fosforilados" (Notheis et al. 1995, comprobado experimentalmente).

Corregirla sería, por tanto, inventar una actividad. Eliminarla sigue el precedente de los propios mantenedores, que en la versión 9.1.0 eliminaron r_0887 con el comentario "no hay evidencia de esta reacción en levadura". La alternativa (corregirla) se deja como pregunta abierta para los mantenedores.

**Efecto en el modelo:** **ningún cambio en crecimiento, esencialidad génica ni en el banco de pruebas de crecimiento.** La única consecuencia es que r_0954 (PMP oxidasa) se queda sin flujo posible: en el modelo solo le llegaba sustrato a través de la reacción imposible. Lo mismo ocurre si se corrige en lugar de eliminarse. Se explica en la issue y se sugiere la alternativa honesta: añadir la entrada de piridoxamina, porque la levadura la usa como fuente de B6.

### Cómo se ha verificado (más de tres veces)

| Comprobación | Resultado |
|---|---|
| ΔG con eQuilibrator 0.8.1 (component contribution) | +374,6 ± 7,1 kJ/mol (1 mM); igual a 30 °C y entre pH 6,7 y 7,2 |
| ΔG con la tabla propia de yeast-GEM (yETFL, otro método de contribución de grupos) | +84,02, que solo es coherente si está en kcal/mol: +351,5 kJ/mol (estándar 1 M). Las otras tres reacciones de Pdx3 cuadran con eQuilibrator con el mismo factor |
| Cálculo a mano de la parte "agua → H2O2" con las tablas NBS | +259 kJ/mol por sí sola, sin depender de ninguna estimación sobre el piridoxal |
| Comprobación con potenciales redox (Wood 1988) | Igual que las tablas NBS (+94,9 kJ/mol por H2O2 en estándar 1 M) |
| Identidad algebraica: reacción actual = −(reacción corregida) + 3 × (H2O + ½ O2 → H2O2) | Cuadra al décimo de kJ con eQuilibrator |
| Balance de masa y carga (cobrapy) | Actual, propuesta y modelo final: equilibrados |
| Modelo editado frente al original menos r_0953 | Idénticos en reacciones, metabolitos, genes y crecimiento |
| Ciclos que generan ATP sin sustrato | Ninguno, en ninguna de las tres variantes |
| Pruebas oficiales de yeast-GEM (esencialidad génica y crecimiento) | Idénticas en las tres variantes (original, eliminada, corregida): crecimiento 0,0811; esencialidad: exactitud 0,9015, MCC 0,5323 y las mismas predicciones gen a gen; R² de crecimiento 0,9019 |

### Lo que no he podido comprobar (y se dice en la issue)

- **Loubbardi et al. 1995**, el artículo del gen PDX3: no pude leer el texto completo (el PDF está detrás de un captcha). Solo cito su resumen y lo indico.
- **El fichero `model/yeast-GEM.yml`:** lo genera su herramienta de MATLAB, que no tenemos. He quitado el bloque de r_0953 a mano y he verificado con cobrapy que el resultado es exactamente el modelo original sin esa reacción. Las herramientas de Python de yeast-GEM no reproducen el fichero byte a byte (añaden cientos de campos vacíos), así que no se pueden usar para esto. Lo explico en el PR y les pido que, si quieren, ejecuten `v9_1_2.m`.

### Pasos para enviarlo

Su guía pide **abrir primero una issue** y después el pull request. Necesitas una cuenta de GitHub y `git`; no hace falta MATLAB.

**Paso 1 — Abrir la issue** (antes de tocar código)
1. Ve a <https://github.com/SysBioChalmers/yeast-GEM/issues/new/choose> y elige **"✏️ Propose curation"**.
2. Título: `r_0953 (pyridoxal oxidase) is thermodynamically infeasible and has no enzymatic support: proposal to remove it`
3. Borra el contenido de la plantilla y pega el texto de la **parte A**.
4. Revisa el párrafo de presentación del principio y ajústalo a tu gusto.
5. Envíala y anota su número (por ejemplo, #430).
6. Opcional: deja un comentario breve en la issue #428 enlazando la nueva ("Reviewed `r_0953` from the alarm list in #NNN").

**Paso 2 — Esperar respuesta.** Puede que prefieran corregirla en lugar de eliminarla (está como pregunta en la issue). Si es así, avísame y preparo esa variante: ya está probada, es la que aparece como "corrected" en las pruebas.

**Paso 3 — Preparar el pull request** (en una terminal)
```bash
# 1. En GitHub: botón "Fork" en https://github.com/SysBioChalmers/yeast-GEM
# 2. Clona TU fork (cambia TU_USUARIO) y sitúate en develop
git clone https://github.com/TU_USUARIO/yeast-GEM.git
cd yeast-GEM
git checkout develop
git checkout -b fix/r_0953-pyridoxal-oxidase

# 3. Aplica el cambio (comprueba primero que encaja)
git apply --check ~/Projects/MetfluxLive/tools/upstream/yeast-GEM-r_0953/r_0953.patch
git apply ~/Projects/MetfluxLive/tools/upstream/yeast-GEM-r_0953/r_0953.patch

# 4. Sustituye #NNN por el número real de la issue (ejemplo con 430)
sed -i 's/#NNN/#430/g' code/modelCuration/v9_1_2.m data/deprecatedIdentifiers/deprecatedReactions.tsv
git diff --stat      # deben salir 4 ficheros: +63 / -19 líneas aprox.

# 5. Commit y push
git add -A
git commit -m "fix: remove thermodynamically infeasible pyridoxal oxidase r_0953" -m "Closes #430"
git push -u origin fix/r_0953-pyridoxal-oxidase
```
Si `git apply --check` falla, es que `develop` ha cambiado desde el 19 de noviembre de 2026 (commit `b54e27f`). Avísame y lo rehago sobre la versión nueva.

**Paso 4 — Abrir el pull request**
1. GitHub te mostrará el botón **"Compare & pull request"**. Comprueba que el destino es `SysBioChalmers/yeast-GEM`, rama **`develop`** (no `main`).
2. Título: `fix: remove thermodynamically infeasible pyridoxal oxidase r_0953`
3. Pega el texto de la **parte B** (cambia #NNN por el número real).
4. Marca **"Allow edits from maintainers"** y la casilla de la plantilla que confirma que el destino es `develop`.
5. Su CI comentará automáticamente con los resultados de calidad. Si algo falla, pásame el enlace.

---

## A. Issue text (paste into "Propose curation")

> I am not an expert in metabolic modelling. While studying the free-energy estimates of yeast-GEM reactions, I asked an AI assistant (Claude, by Anthropic) to review this reaction and to check the calculations thoroughly and independently several times. Everything below can be reproduced with the scripts listed at the end. I apologise in advance if I have missed something obvious.

#### Brief description:

[r_0953](https://metabolicatlas.org/explore/Yeast-GEM/gem-browser/reaction/r_0953) (pyridoxal oxidase, PDX3) is one of the 22 alarms listed in #428 that have not been reviewed yet. As written, it reduces pyridoxal to pyridoxamine with electrons taken from water, which is oxidised to H2O2. The reaction is mass- and charge-balanced, but its ΔG'm is **+375 ± 7 kJ/mol**, so it cannot run in the only direction its bounds allow. Correcting it to the non-phosphorylated analogue of r_0954 is not supported by any enzymatic evidence, so I propose to remove it.

#### Current status:

Identifier | Reaction equation | Genes | Bounds
---|---|---|---
r_0953 | ammonium[c] + 2 H2O[c] + 0.5 oxygen[c] + pyridoxal[c] => 2 hydrogen peroxide[c] + pyridoxamine[c] | YBR035C | [0, 1000]

Annotations: BiGG `PYDXO`, MetaNetX `MNXR103365`, EC 1.4.3.5, no KEGG reaction. The same stoichiometry is in BiGG `PYDXO` (iND750, iMM904), ModelSEED `rxn09566` and MetaNetX `MNXR103365`. There is no matching Rhea or KEGG reaction.

The other PDX3 reactions are thermodynamically sound and are not affected:
- [r_0954](https://metabolicatlas.org/explore/Yeast-GEM/gem-browser/reaction/r_0954): pyridoxamine 5'-phosphate + O2 + H2O => PLP + NH4+ + H2O2, ΔG'm = −114 ± 5 kJ/mol
- [r_0955](https://metabolicatlas.org/explore/Yeast-GEM/gem-browser/reaction/r_0955): pyridoxine 5'-phosphate + O2 => PLP + H2O2, −113 ± 4 kJ/mol
- [r_0956](https://metabolicatlas.org/explore/Yeast-GEM/gem-browser/reaction/r_0956): pyridoxine + O2 => pyridoxal + H2O2, −106 ± 5 kJ/mol

**Evidence**

- **The reaction cannot run forward.** With eQuilibrator 0.8.1 (component contribution) at the conditions proposed in #428 (pH 7.2, I = 0.25 M, pMg 3, 25 °C), ΔG'° = +383.2 ± 7.1 and ΔG'm = +374.6 ± 7.1 kJ/mol. At 30 °C and at pH 6.7 the values change by less than 0.5 kJ/mol. Each tenfold change in the reaction quotient shifts ΔG' by 5.7 kJ/mol, so concentrations would have to be about 10^66-fold more favourable than 1 mM.
- **The conclusion does not depend on estimates for the B6 vitamers.** The reaction is exactly −(pyridoxamine + O2 + H2O → pyridoxal + NH4+ + H2O2) + 3 × (H2O + ½ O2 → H2O2). The second term uses only water, O2 and H2O2, whose formation energies are measured (NBS tables, Wagman et al. 1982: H2O(l) −237.129, H2O2(aq) −134.03, O2(aq) +16.4 kJ/mol). This gives ΔG° = +94.9 kJ/mol per H2O2 (eQuilibrator: +95.9 ± 3.4), and +86.3 kJ/mol at 1 mM. Three of them already give **+259 kJ/mol**, whatever the vitamer chemistry. The same value follows from the pH-7 potentials of the O2/H2O2 and H2O2/H2O couples ([Wood 1988](https://doi.org/10.1042/bj2530287)).
- **An independent estimate agrees.** The repository's own table `data/databases/model_rxnDeltaG.tsv` (from yETFL, a different group-contribution method) lists r_0953 at 84.02. Its units are not stated (as #428 notes), but the values are consistent with kcal/mol: 84.02 kcal/mol = +351.5 kJ/mol, against eQuilibrator's +383 ± 7 kJ/mol for ΔG'°. With the same conversion, r_0954, r_0955 and r_0956 (−105.8, −118.4 and −117.8 kJ/mol) are within 5–12 kJ/mol of eQuilibrator.
- **Enzyme specificity: PDX3 acts on the 5'-phosphates.** The IUBMB entry for [EC 1.4.3.5](https://www.enzyme-database.org/query.php?ec=1.4.3.5) lists only pyridoxamine 5'-phosphate and pyridoxine 5'-phosphate as substrates. So do [Rhea](https://www.rhea-db.org/rhea/15817) ([RHEA:15817](https://www.rhea-db.org/rhea/15817), [RHEA:15149](https://www.rhea-db.org/rhea/15149)) and UniProt [P38075](https://www.uniprot.org/uniprotkb/P38075/entry) (PDX3_YEAST). Every one of the ~360 reviewed UniProt entries with EC 1.4.3.5 is annotated only with the phosphorylated substrates.
- **The E. coli orthologue does not oxidise non-phosphorylated vitamers.** Purified PdxH oxidised PNP and PMP and "did not accept non-phosphorylated substrates" ([Notheis et al. 1995](https://doi.org/10.1016/0167-4838(94)00235-9); tested directly). Nothing equivalent has been measured for the yeast enzyme, as far as I could find. The inference that Pdx3 behaves the same is based on homology (Loubbardi et al. report high sequence homology with the E. coli enzyme).
- **Yeast uses pyridoxamine without PDX3.** pdx3 (aux30) mutants lack P(N/M)P oxidase activity. They fail to grow in YNB with pyridoxine, but grow normally when pyridoxine is replaced by pyridoxal or pyridoxamine ([Loubbardi et al. 1995](https://doi.org/10.1128/jb.177.7.1817-1823.1995); I could only read the abstract). The direct observation is growth. My inference is that turning pyridoxamine into PLP in vivo does not require a Pdx3 activity on free pyridoxamine.

#### Proposed status/to-do list:

Identifier | Reaction equation | Genes
---|---|---
r_0953 | *(removed)* |

- [ ] Remove r_0953 (`removeReactions(model,'r_0953',true,true,true)`, as was done for r_0887 in v9_1_0.m). No metabolite or gene is removed with it: pyridoxamine[c] remains in r_2026 and r_2027, and YBR035C in r_0954–r_0956.
- [ ] Remove its row from `model/reactions.tsv` and add it to `data/deprecatedIdentifiers/deprecatedReactions.tsv`.
- [ ] Leave `data/databases/BiGGrxnDictionary.csv` and `model_rxnDeltaG.tsv` untouched, as for r_0887 (happy to prune them if you prefer).

**Other reactions with the same error:** among the vitamin B6 reactions, none is thermodynamically infeasible. r_0956 (pyridoxine oxidase, PDX3) has the same lack of evidence for a non-phosphorylated substrate, but it is thermodynamically sound, so I left it alone (see the question below).

**Counter-arguments considered**

- **"Correct it instead of removing it: pyridoxamine + O2 + H2O => pyridoxal + NH4+ + H2O2."** That reaction is chemically sound (ΔG'm = −112 ± 5 kJ/mol), mass- and charge-balanced, and I tested it in the model too (results below). But it would give PDX3 an activity that, as far as I can find, has not been shown for any PNPO, and that was ruled out for E. coli PdxH. I am happy to change the proposal if you prefer to keep the vitamer conversion.
- **"Removing it loses coverage of B6 salvage."** It does, in one respect: r_0954 (PMP oxidase) becomes blocked, because in the model r_0953 was the only source of pyridoxamine and PMP. Pyridoxamine kinase (r_2026) and PMP phosphatase (r_2027) form a closed loop with no other entry, and there is no pyridoxamine uptake. The corrected variant blocks r_0954 too, since it consumes pyridoxamine instead of producing it. In cells, PMP is formed from PLP in the half-reactions of PLP-dependent aminotransferases, which the model does not represent. So r_0954 was reachable only through an impossible reaction. If you want it to stay active, a pyridoxamine uptake reaction would be the honest way, since yeast grows on pyridoxamine as its B6 source (Loubbardi et al. 1995).
- **"Maybe the intended reaction was a transamination (pyridoxal + amino acid → pyridoxamine + 2-oxo acid)."** Pyridoxamine-pyruvate transaminase (EC 2.6.1.30) is known from bacteria, not yeast, and r_0953 carries an oxidase EC number and the oxidase gene. So I did not consider it the intended reaction.

#### Expected effect

No change in growth, gene essentiality or the growth benchmark (tables in the appendix, A5). The canonical medium contains no B6 vitamer, so r_0953 carries no flux there. Even with pyridoxine supplied, the vitamer conversions are not needed for growth, because PLP is not in the biomass. The only structural change is that r_0954 becomes blocked: it was reachable only through r_0953 (see the second counter-argument). The number of blocked reactions with every boundary open goes from 991 to 992; nothing becomes unblocked.

#### Question for the maintainers

Would you prefer removal (proposed) or correction to the non-phosphorylated pyridoxamine oxidase? The same question applies to r_0956 (pyridoxine oxidase, PDX3), which also assumes a non-phosphorylated substrate.

#### Related

- #428: reversibility check; r_0953 is in its alarm list. While testing I also noticed something about the PLP synthase reactions r_4211 and r_4212, also alarms in #428 (+305 and +307 kJ/mol backwards). With every boundary reaction open, they cannot run in their biosynthetic direction: r_4211 carries no flux and r_4212 only runs backwards. Running backwards, r_4212 breaks PLP down into ribulose 5-phosphate, glyceraldehyde 3-phosphate and glutamine. The model can therefore use pyridoxine as a carbon and nitrogen source: allowing 1 mmol/gDW/h of pyridoxine uptake raises growth on the canonical medium from 0.0811 to 0.2057 h⁻¹, through r_4330 → r_2025 → r_4212 (backwards). This is independent of r_0953: it is the same with r_0953 removed, and with r_0953 and r_0956 both switched off. I have not looked further into it and mention it as a possible next item from the alarm list.

#### Reproducibility

eQuilibrator-api 0.8.1, equilibrator-cache 0.8.1, component-contribution 0.8.0, cobrapy 0.32.1 (GLPK), yeast-GEM `develop` at b54e27f. The scripts (ΔG estimates, FVA, energy-generating cycle test, `essential_genes` and `growth` benchmarks for the three variants, and the equality check of the edited YAML) can be attached or shared on request.

---

## B. Pull request text (target: `develop`)

```
Closes #NNN

Removes r_0953 (pyridoxal oxidase, PYDXO). As written it oxidises water to H2O2 to reduce pyridoxal to pyridoxamine (ΔG'm = +375 ± 7 kJ/mol), and there is no evidence that PDX3 acts on non-phosphorylated pyridoxamine. Arguments and evidence are in #NNN.

| Reaction | Name | Current | New |
| --- | --- | --- | --- |
| `r_0953` | pyridoxal oxidase | `ammonium[c] + 2 H2O[c] + 0.5 oxygen[c] + pyridoxal[c] => 2 hydrogen peroxide[c] + pyridoxamine[c]` (YBR035C) | removed |

- `code/modelCuration/v9_1_2.m`: new curation script on top of 9.1.1, which removes r_0953 (`removeReactions(model,'r_0953',true,true,true)`, as for r_0887) and saves with `saveYeastYaml`.
- `model/yeast-GEM.yml`: r_0953 removed.
- `model/reactions.tsv`: row of r_0953 removed.
- `data/deprecatedIdentifiers/deprecatedReactions.tsv`: r_0953 added.

Note on how the YAML was made: I do not have MATLAB, so I removed the r_0953 block from `model/yeast-GEM.yml` by hand instead of running `v9_1_2.m`. I checked with cobrapy that the result equals the `develop` model with r_0953 removed (identical reactions, metabolites, genes, bounds, GPRs and growth). Saving through `yeastgem.save_yeast_yaml` (raven-toolbox 0.4.0) would have rewritten ~1400 unrelated lines, so I did not use it. Please re-run `v9_1_2.m` if you prefer the file to come from the script.

Benchmarks (Python `yeastgem.model_tests`, GLPK), develop → this PR:
- growth on the canonical medium: 0.081100 → 0.081100
- gene essentiality: accuracy 0.9015 → 0.9015, MCC 0.5323 → 0.5323; per-gene calls identical
- growth benchmark (Tobias 2013) R²: 0.9019 → 0.9019
```

---

## Appendix: calculations in full

### A1. eQuilibrator 0.8.1 (component contribution)

Conditions: I = 0.25 M, pMg 3. The compounds are matched by KEGG identifiers: pyridoxal C00250, pyridoxamine C00534, PMP C00647, PLP C00018, pyridoxine C00314, PNP C00627, NH3 C00014, H2O C00001, O2 C00007, H2O2 C00027. Script: `verification/thermo_equilibrator.py`. Values in kJ/mol (± 1σ); ΔG'm has every reactant at 1 mM except water.

| Reaction | ΔG'° (pH 7.2, 25 °C) | ΔG'm (pH 7.2, 25 °C) | ΔG'm (pH 7.2, 30 °C) | ΔG'm (pH 6.7, 25 °C) |
|---|---|---|---|---|
| r_0953 as written | +383.2 ± 7.1 | **+374.6 ± 7.1** | +374.5 | +374.5 |
| r_0953 corrected (PM + O2 + H2O → PL + NH4+ + H2O2) | −95.3 ± 5.1 | −112.4 ± 5.1 | −112.8 | −112.3 |
| r_0954 (PMP oxidase) | −97.0 ± 4.6 | −114.1 ± 4.6 | −114.6 | −113.5 |
| r_0955 (PNP oxidase) | −113.4 ± 4.2 | −113.4 ± 4.2 | −113.8 | −112.2 |
| r_0956 (PN oxidase) | −106.3 ± 4.7 | −106.3 ± 4.7 | −106.5 | −105.7 |
| H2O + ½ O2 → H2O2 | +95.9 ± 3.4 | +87.4 ± 3.4 | +87.2 | +87.4 |

### A2. The oxygen term by hand (independent of group contribution)

Formation energies (NBS tables, Wagman et al., *J. Phys. Chem. Ref. Data* 11, Suppl. 2, 1982), 25 °C, 1 M standard state for solutes:
- H2O(l): −237.129 kJ/mol
- H2O2(aq, undissociated; pKa 11.6, so it is undissociated at pH 7): −134.03 kJ/mol
- O2(aq): +16.4 kJ/mol

ΔG° (H2O + ½ O2(aq) → H2O2(aq)) = −134.03 − (−237.129) − ½ (16.4) = **+94.90 kJ/mol**. No proton is involved, so the value does not depend on pH.

At 1 mM (water at unit activity): ΔG'm = ΔG° + RT ln(10⁻³ / (10⁻³)^½) = 94.90 + (2.4790 × −3.4539) = 94.90 − 8.56 = **+86.34 kJ/mol** (eQuilibrator: +87.4 ± 3.4).

Cross-check with redox potentials at pH 7 (Wood 1988, O2 at 1 atm): E'°(O2/H2O2) = +0.281 V and E'°(H2O2/2 H2O) = +1.349 V.
- O2(g) + 2 H2O → 2 H2O2: ΔG'° = −2F(0.281 − 1.349) = −2 × 96.485 × (−1.068) = +206.1 kJ/mol, that is +103.05 per H2O2 with ½ O2(g).
- Changing O2 from gas (1 atm) to solute (1 M) takes away ½ × 16.4 = 8.2, which gives **+94.85 kJ/mol**. This agrees with the 16.4 kJ/mol that also follows from the solubility of O2, about 1.3 mM/atm: −RT ln(1.3 × 10⁻³) = +16.4 kJ/mol.
- These potentials come from the same thermochemical data as the NBS values. This is a consistency check, not an independent measurement.

### A3. Decomposition

r_0953 (as written) = −[PM + O2 + H2O → PL + NH4+ + H2O2] + 3 × [H2O + ½ O2 → H2O2]

Check of the stoichiometry: PL + NH4+ + H2O2 + 3 H2O + 1.5 O2 → PM + O2 + H2O + 3 H2O2, which simplifies to PL + NH4+ + 2 H2O + 0.5 O2 → PM + 2 H2O2 ✓.

- With eQuilibrator values: +112.4 + 3 × 87.4 = **+374.6 kJ/mol**, identical to the direct estimate in A1.
- With the NBS oxygen term: +112.4 + 3 × 86.34 = **+371.4 kJ/mol**.
- Lower bound if the vitamer part were assumed to be 0 instead of +112: **+259 kJ/mol**.

How far concentrations can shift it: ΔG' = ΔG'm + RT ln(Q/Q₁ₘₘ), with RT ln 10 = 5.71 kJ/mol at 25 °C. Bringing +374.6 down to 0 needs Q/Q₁ₘₘ = 10^(−374.6/5.71) ≈ 10⁻⁶⁶.

### A4. yETFL table in the repository (`data/databases/model_rxnDeltaG.tsv`)

The values are read as kcal/mol and converted with 4.184 kJ/kcal. The table does not state its units. Reading them as kcal/mol is an inference supported by the agreement in the last column: read as kJ/mol, all four would disagree with eQuilibrator by a factor of about 4.

| Reaction | Table | × 4.184 (kJ/mol) | eQuilibrator ΔG'° | Difference |
|---|---|---|---|---|
| r_0953 | 84.02 | +351.5 | +383.2 | 32 |
| r_0954 | −25.29 | −105.8 | −97.0 | 9 |
| r_0955 | −28.31 | −118.4 | −113.4 | 5 |
| r_0956 | −28.15 | −117.8 | −106.3 | 12 |

### A5. Model tests (cobrapy 0.32.1, GLPK, yeast-GEM `develop` b54e27f)

- **Mass and charge balance** (`check_mass_balance`): balanced for r_0953 as written and as corrected.
- **Can r_0953 carry flux?** (`verification/fva_permissive.py`)
  - With every boundary reaction open both ways and the reactions that are switched off by default opened, r_0953 ranges over [0, 1000].
  - In the canonical medium (no B6 vitamer supplied) it carries no flux, because there is no source of pyridoxal.
  - Pyridoxamine[c] (s_1394) takes part only in r_2027 (made), r_2026 (used) and r_0953.
- **PLP and biomass:** PLP (s_1393) is not in the biomass or cofactor pseudoreactions.
- **Energy-generating cycles** (`verification/energy_cycles.py`): with no uptake at all, and with only O2 and H2O uptake, maximal ATP maintenance (r_4046) is 0 for develop, for r_0953 removed and for r_0953 corrected.
- **The hand-edited YAML** (`verification/semantic_equality.py`): equals develop with r_0953 removed by cobrapy. Reactions (4104), metabolites (2748), genes (1143), bounds, GPRs and subsystems are identical, and growth is the same (0.0811003).
- **Blocked reactions** (`verification/blocked_reactions.py`, all boundaries open): develop 991, removed 992 (+ r_0954), corrected 993 (+ r_0953, r_0954); none becomes unblocked.
- **Pyridoxine as a carbon source** (`verification/pyridoxine_growth.py`): growth 0.0811 → 0.2057 h⁻¹ with 1 mmol/gDW/h pyridoxine, identical in develop, with r_0953 removed, and with r_0953 and r_0956 off. The flux runs r_4330 → r_2025 → r_4212 backwards.
- **yeast-GEM benchmarks** (`verification/benchmarks.py`, `yeastgem.model_tests`, `processes=1`):

| Variant | Growth, canonical medium | Growth + pyridoxine (1 mmol/gDW/h) | Essential genes: accuracy / MCC | TP / TN / FP / FN | Growth R² (Tobias 2013) |
|---|---|---|---|---|---|
| develop | 0.081100 | 0.205670 | 0.9015 / 0.5323 | 933 / 65 / 94 / 15 | 0.9019 |
| r_0953 removed (this proposal) | 0.081100 | 0.205670 | 0.9015 / 0.5323 | 933 / 65 / 94 / 15 | 0.9019 |
| r_0953 corrected (alternative) | 0.081100 | 0.205670 | 0.9015 / 0.5323 | 933 / 65 / 94 / 15 | 0.9019 |

The per-gene calls (the TP, TN, FP and FN sets) are identical to develop in both variants.

### A6. Sources

- Wagman D.D. et al. (1982) The NBS tables of chemical thermodynamic properties. *J. Phys. Chem. Ref. Data* 11, Suppl. 2.
- Wood P.M. (1988) The potential diagram for oxygen at pH 7. *Biochem. J.* 253:287–289. https://doi.org/10.1042/bj2530287
- Notheis C., Drewke C., Leistner E. (1995) Purification and characterization of the pyridoxol-5'-phosphate:oxygen oxidoreductase (deaminating) from *Escherichia coli*. *Biochim. Biophys. Acta* 1247:265–271. https://doi.org/10.1016/0167-4838(94)00235-9
- Loubbardi A. et al. (1995) Sterol uptake induced by an impairment of pyridoxal phosphate synthesis in *Saccharomyces cerevisiae*: cloning and sequencing of the PDX3 gene encoding pyridoxine (pyridoxamine) phosphate oxidase. *J. Bacteriol.* 177:1817–1823. https://doi.org/10.1128/jb.177.7.1817-1823.1995 (abstract only)
- Beber M.E. et al. (2022) eQuilibrator 3.0. *Nucleic Acids Res.* 50:D603–D609. https://doi.org/10.1093/nar/gkab1106
- IUBMB EC 1.4.3.5: https://www.enzyme-database.org/query.php?ec=1.4.3.5 ; Rhea RHEA:15817, RHEA:15149; UniProt P38075.
- BiGG `PYDXO`: http://bigg.ucsd.edu/universal/reactions/PYDXO ; MetaNetX `MNXR103365`; ModelSEED `rxn09566`.
