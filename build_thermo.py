"""Gibbs free energies of reaction for MetfluxLive (data/<slug>/thermo.json).

    tools/.venv-thermo/bin/python tools/build_thermo.py human-gem

Uses eQuilibrator 3 / component contribution (Beber et al. 2022; Noor et al. 2013), run locally.
The first run downloads the eQuilibrator compound database from Zenodo (doi 10.5281/zenodo.4128543)
into the user cache; nothing about the model is sent anywhere.

For every reaction whose metabolites are all found in eQuilibrator (by MetaNetX, ChEBI, HMDB, BiGG or KEGG
identifier, else by the InChIKey of the model's SMILES) it stores dG'm: the transformed Gibbs energy with every reactant at 1 mM (water at 1), at
37 °C, ionic strength 0.25 M, pMg 3 and the pH of its compartment (COND below). Reactions that span two
compartments are computed with the multicompartment method of Haraldsdottir et al. 2012 (as implemented
in eQuilibrator), which adds the work of moving charge across the membrane potential and protons across
the pH difference. The cytochrome c and ubiquinone couples use measured midpoint potentials (COUPLES). Exchange, biomass and pseudo-reactions are left out.

Output: { "cond": {...}, "rxns": { "<rxn id>": [dGm, sigma] } } in kJ/mol, for the forward direction
as written in the model.
"""
import json
import math
import sys
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parent.parent

# pH: Casey, Grinstein & Orlowski 2010 (Nat Rev Mol Cell Biol 11:50) except the intermembrane space
# (Porcelli et al. 2005, Biochem Biophys Res Commun 326:799). Golgi: mid-way between cis (6.7) and trans (6.0).
PH = {"c": 7.2, "e": 7.4, "m": 8.0, "i": 6.9, "n": 7.2, "r": 7.2, "g": 6.4, "l": 4.7, "x": 7.0}
# Electrical potential relative to the cytosol (V). Mitochondrial matrix -0.16 (typical 150-180 mV,
# Zorova et al. 2018); plasma membrane: cytosol at -0.06 V against the medium (Yang & Brackenbury 2013
# give -10 to -90 mV depending on the cell). Intermembrane space at cytosolic potential (porous outer
# membrane). Other organelles are set to 0: their potentials are small or poorly known.
PSI = {"c": 0.0, "e": 0.06, "m": -0.16, "i": 0.0, "n": 0.0, "r": 0.0, "g": 0.0, "l": 0.0, "x": 0.0}
IONIC, PMG, TEMP = 0.25, 3.0, 310.15
# Redox couples taken from measured midpoint potentials instead of eQuilibrator, when both forms appear in a
# reaction: (oxidised, reduced, electrons n, protons h, E'0 at pH 7 in V, dE/dpH in V).
#  - cytochrome c is a protein, absent from eQuilibrator: +0.25 V, pH independent (canonical mitochondrial
#    cytochrome c, Loveridge et al. 2026, J Mol Biol 438:169996).
#  - ubiquinone/ubiquinol: eQuilibrator's group-contribution estimate for the long isoprenoid quinone puts the
#    couple near +0.41 V (succinate dehydrogenase would be +74 kJ/mol); measured in submitochondrial
#    particles: +0.065 V at pH 7, -60 mV per pH unit (Urban & Klingenberg 1969, Eur J Biochem 9:519).
# The last field is where the couple really sits: Human-GEM files cytochrome c under the matrix, but it is
# in the intermembrane space, so the electrons it takes or gives cross the inner membrane (work -F*dpsi each).
COUPLES = [("ferricytochrome C", "ferrocytochrome C", 1, 0, 0.25, 0.0, "i"),
           ("ubiquinone", "ubiquinol", 2, 2, 0.065, -0.060, None)]
FARADAY = 96.485  # kJ/(mol V)
MAX_SIGMA = 1000.0  # eQuilibrator reports ~1e5 kJ/mol when a group is not covered: no estimate


def main():
    slug = sys.argv[1] if len(sys.argv) > 1 else "human-gem"
    from equilibrator_api import ComponentContribution, Q_
    from equilibrator_api.phased_reaction import PhasedReaction

    cc = ComponentContribution()
    cc.ionic_strength = Q_(f"{IONIC}M")
    cc.p_mg = Q_(PMG)
    cc.temperature = Q_(f"{TEMP}K")
    RT = 8.314462618e-3 * TEMP

    model = json.loads((ROOT / "data" / slug / "model.json").read_text())
    bases, mets, rxns = model["bases"], model["mets"], model["rxns"]

    # InChIKeys from the model's own SMILES (tools/cache/<slug>/smiles.json, written by build_model.py)
    smiles = {}
    sp = ROOT / "tools" / "cache" / slug / "smiles.json"
    if sp.exists():
        smiles = json.loads(sp.read_text())
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
    except ImportError:
        Chem = None

    def by_inchikey(base_id):
        smi = smiles.get(base_id)
        if not smi or Chem is None:
            return None
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        key = Chem.MolToInchiKey(mol)
        if not key:
            return None
        # same skeleton and stereo, any protonation state; then same skeleton only (component contribution
        # ignores stereochemistry, so stereoisomers get the same estimate)
        for prefix in (key[:25], key[:14]):
            ids = key_index().get(prefix)
            if ids:
                return cc.ccache.get_compound_by_internal_id(min(ids))
        return None

    _keys = {}

    def key_index():
        """InChIKey prefixes (25 and 14 characters) -> compound ids, read once from the cache database."""
        if not _keys:
            from sqlalchemy import text
            rows = cc.ccache.session.execute(text("SELECT id, inchi_key FROM compounds WHERE inchi_key IS NOT NULL"))
            for cid, k in rows:
                for pre in (k[:25], k[:14]):
                    _keys.setdefault(pre, []).append(cid)
            print(f"  InChIKey index: {len(_keys)} prefixes", flush=True)
        return _keys

    def lookup(b):
        keys = []
        if b.get("mnx"):
            keys.append(f"metanetx.chemical:{b['mnx']}")
        if b.get("chebi"):
            keys.append(f"chebi:{b['chebi'] if b['chebi'].startswith('CHEBI:') else 'CHEBI:' + b['chebi']}")
        if b.get("hmdb"):
            keys.append(f"hmdb:{b['hmdb']}")
        if b.get("bigg"):
            keys.append(f"bigg.metabolite:{b['bigg']}")
        if b.get("kegg"):
            keys.append(f"kegg:{b['kegg']}")
        for k in keys:
            try:
                c = cc.get_compound(k)
            except Exception:
                c = None
            if c is not None and c.inchi_key:
                return c
        return by_inchikey(b["id"])

    t0 = time.time()
    cpd = [lookup(b) for b in bases]
    found = sum(c is not None for c in cpd)
    print(f"compounds: {found}/{len(bases)} found ({time.time() - t0:.0f} s)", flush=True)
    water = next((c for c, b in zip(cpd, bases) if c is not None and b["n"] == "H2O"), None)
    proton = next((c for c, b in zip(cpd, bases) if c is not None and b["n"] == "H+"), None)

    byname = {b["n"]: i for i, b in enumerate(bases)}
    couples = [(byname[o], byname[r], n, h, e0, slope, at) for o, r, n, h, e0, slope, at in COUPLES if o in byname and r in byname]
    electron = None

    def conc_correction(sparse):
        """RT * sum(nu * ln 1 mM) over reactants other than water, protons and electrons."""
        s = sum(nu for c, nu in sparse.items() if c is not water and c is not proton and c is not electron)
        return RT * s * math.log(1e-3)

    out, why = {}, {"exchange": 0, "unmapped": 0, "unbalanced": 0, "uncovered": 0, "multi": 0, "error": 0}
    t0 = time.time()
    for j, r in enumerate(rxns):
        if len(r["m"]) < 2 or r.get("obj") or "biomass" in r["n"].lower():
            why["exchange"] += 1
            continue
        by = {}
        ok = True
        extra = 0.0  # kJ/mol added for the redox couples
        electrons = 0.0  # electrons left out of the reaction when a couple is replaced
        skip = set()
        rb = [mets[m]["b"] for m in r["m"]]
        for ox, red, n, h, e0, slope, at in couples:
            if ox in rb and red in rb:
                # R = R'' + nu * (ox + h H+ + n e- -> red): R'' has the couple removed and h*nu protons added in
                # the compartment of the reduced form; the electrons have formation energy 0 and are dropped
                q = rb.index(red)
                nu, comp = r["c"][q], mets[r["m"][q]]["c"]
                e_comp = e0 + slope * (PH.get(comp, 7.2) - 7.0)
                extra += nu * -n * FARADAY * e_comp
                if at and at != comp:  # electrons (charge -1) moving from the reaction's compartment to the couple's
                    extra += nu * n * -FARADAY * (PSI.get(at, 0.0) - PSI.get(comp, 0.0))
                electrons += n * nu
                if h:
                    side = by.setdefault(comp, {})
                    side[proton] = side.get(proton, 0.0) + h * nu
                skip.update((ox, red))
        for m, nu in zip(r["m"], r["c"]):
            if mets[m]["b"] in skip:
                continue
            c = cpd[mets[m]["b"]]
            if c is None:
                ok = False
                break
            side = by.setdefault(mets[m]["c"], {})
            side[c] = side.get(c, 0.0) + nu
        if not ok:
            why["unmapped"] += 1
            continue
        for side in by.values():
            for c in [c for c, nu in side.items() if abs(nu) < 1e-9]:
                del side[c]
        comps = [k for k, side in by.items() if side]
        try:
            if len(comps) == 1:
                k = comps[0]
                rx = PhasedReaction(by[k])
                if not rx.is_balanced(ignore_atoms=("H", "e-") if skip else ("H",)):
                    why["unbalanced"] += 1
                    continue
                cc.p_h = Q_(PH.get(k, 7.2))
                dg = cc.standard_dg_prime(rx)
                val, sd = dg.value.m_as("kJ/mol"), dg.error.m_as("kJ/mol")
                val += conc_correction(by[k]) + extra
            elif len(comps) == 2:
                a, b = comps
                # the inner compartment is the one with the more negative potential (then the higher pH)
                inner, outer = sorted((a, b), key=lambda k: (PSI.get(k, 0.0), -PH.get(k, 7.2)))
                total = {}
                for side in by.values():
                    for c, nu in side.items():
                        total[c] = total.get(c, 0.0) + nu
                if not PhasedReaction(total).is_balanced(ignore_atoms=("H", "e-") if skip else ("H",)):
                    why["unbalanced"] += 1
                    continue
                cc.p_h = Q_(PH.get(inner, 7.2))
                dg = cc.multicompartmental_standard_dg_prime(
                    PhasedReaction(by[inner]), PhasedReaction(by[outer]),
                    # eQuilibrator's sign: psi(inner) - psi(outer). Checked with H+ moving from the intermembrane
                    # space into the matrix: -22 kJ/mol with -0.16 V (15.4 electrical + 6.5 from the pH step)
                    e_potential_difference=Q_(f"{PSI.get(inner, 0.0) - PSI.get(outer, 0.0)}V"),
                    p_h_outer=Q_(PH.get(outer, 7.2)), ionic_strength_outer=Q_(f"{IONIC}M"),
                    p_mg_outer=Q_(PMG), tolerance=0.5 + abs(electrons))
                val, sd = dg.value.m_as("kJ/mol"), dg.error.m_as("kJ/mol")
                val += conc_correction(by[inner]) + conc_correction(by[outer]) + extra
            else:
                why["multi"] += 1
                continue
        except Exception as err:
            why["error"] += 1
            if why["error"] <= 5:
                print("  error", r["id"], err)
            continue
        if not (math.isfinite(val) and math.isfinite(sd)) or sd > MAX_SIGMA:
            why["uncovered"] += 1
            continue
        out[r["id"]] = [round(val, 1), round(sd, 1)]
        if j % 1000 == 0:
            print(f"  {j}/{len(rxns)}  {len(out)} done ({time.time() - t0:.0f} s)", flush=True)

    print(f"reactions with dG'm: {len(out)}/{len(rxns)}; skipped: {why}")
    doc = {
        "source": "eQuilibrator 3 (component contribution), compound database doi 10.5281/zenodo.4128543",
        "cond": {"pH": PH, "psi": PSI, "ionicStrength": IONIC, "pMg": PMG, "T": TEMP, "conc": 1e-3},
        "rxns": out,
    }
    (ROOT / "data" / slug / "thermo.json").write_text(json.dumps(doc, separators=(",", ":")))


if __name__ == "__main__":
    main()
