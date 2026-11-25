# Human-GEM 2.1.1: half-reactions with no electron acceptor allow unlimited anaerobic NAD(P)H oxidation

Found while validating the exometabolomics fit (`js/core/fitexchange.js`). Reproduce with
`node tools/validation/humangem_charge_cycle.mjs`.

## What happens

Four reactions close a cycle whose net equation is

    NAD(P)H  ->  NAD(P)+ + H+

so the model can oxidise NAD(P)H with no electron acceptor at all. Combined with glycolysis and the
pentose phosphate pathway it performs

    glucose + 6 H2O  ->  6 CO2 + 24 H+          (no O2)

In our test case (HEK293 measured exchange rates, glucose uptake deliberately set 5x too high) the cycle
runs at 20.35 mmol/gDW/h, the cell secretes 41.1 mmol/gDW/h of H+ and the respiratory quotient reaches
CO2/O2 = 33.6 with the oxygen uptake pinned at 0.335.

## The cycle

| Reaction | Name | Equation |
|---|---|---|
| `MAR03167` | glutathione transferase (urate) | urate [C5H4N4O3] -> urate radical [C5H3N4O3] + H+ |
| `MAR06606` | ascorbate to monodehydroascorbate conversion | ascorbate + urate radical -> monodehydroascorbate + urate |
| `MAR03265` | prostaglandin-F synthase (monodehydroascorbate) | monodehydroascorbate [C6H6O6] -> dehydroascorbic acid [C6H5O6] + H+ |
| `MAR06396` | glutaredoxin (dehydroascorbic acid) | dehydroascorbic acid + H+ + NADPH -> ascorbate + NADP+ |

`MAR03868` (glutathione:NAD+ oxidoreductase) with `MAR08620` (glutathione:dehydroascorbate oxidoreductase)
does the same with NADH, at 3.50 mmol/gDW/h in the same solution.

## Why a stoichiometric check does not catch it

`MAR03167` and `MAR03265` balance in C, H, N and O. A hydrogen **atom** on the left becomes a **proton**
on the right, and the electron is simply dropped: each reaction creates +1 of charge out of nothing. Only a
charge balance sees it, and the formulas in `data/<slug>/model.json` carry no charge, so our own element
check reports these reactions as balanced.

Independent corroboration from this repository: `data/human-gem/thermo.json`, built with eQuilibrator 3
(component contribution), has **no ΔG for MAR03167, MAR06606, MAR03265, MAR03868 or MAR08620** — only
`MAR06396` gets one (6.2 ± 1.6 kJ/mol). eQuilibrator refuses a reaction that does not balance, so the
absence of a free energy is itself evidence.

## Chemistry

Urate and ascorbate radicals are real intermediates, but a one-electron oxidation needs an acceptor:

    urate + A  ->  urate radical + AH           (not: urate -> urate radical + H+)

The two reactions seem to have been written as the proton-releasing half of a redox pair whose other half
(the acceptor) was left out.

## Suggested fix upstream

Either give both reactions an explicit acceptor, or mark them as half-reactions and remove them from the
model. Since `MAR06606` already couples the urate radical to ascorbate and `MAR06396`/`MAR08620` regenerate
ascorbate at the cost of NADPH or GSH, deleting `MAR03167` and `MAR03265` would leave the ascorbate–urate
antioxidant chemistry intact while closing the cycle.

## How this affects MetfluxLive

Nothing is changed in the model: the defect only shows up when exchange fluxes are forced far from a
physiological solution, which is exactly what the consistency fit does. `exchangeBalance()` therefore
reports the respiratory quotient, the electron budget (counting a proton as carrying **no** electrons,
unlike its formula "H") and the net proton flux, and warns when the proton efflux exceeds the carbon
uptake. Those checks are what make the fit's verdict trustworthy.
