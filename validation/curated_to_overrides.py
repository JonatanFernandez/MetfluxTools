#!/usr/bin/env python3
"""Turn curated Recon 3D atom mappings (Preciat Gonzalez et al. 2017) into rules of
tools/atommap_overrides.json for the Human-GEM reactions where the validation found our map wrong.

Each curated molecule is matched to the Human-GEM metabolite on the same side with the same skeleton
(heavy atoms, element, H count on carbon). Among the graph isomorphisms the one that preserves the
most CIP labels (R/S) is used, so molecules that differ only in stereochemistry (myo-inositol phosphate)
are not mirrored. Atoms that share a map number correspond; equivalent atoms are paired in order.

Usage: curated_to_overrides.py <standardisedPredictions dir> RECON_ID:MAR_ID [...]
"""
import json, sys
from collections import defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler  # noqa: F401  (CIP labels)

src = open(Path(__file__).with_name("atommaps_vs_recon3d.py")).read()
src = src.split("# ---------------------------------------------------------------- Recon3D id")[0]
sys_argv = sys.argv
sys.argv = [sys_argv[0], sys_argv[1]]
exec(compile(src, "v", "exec"))
pairs_wanted = [x.split(":") for x in sys_argv[2:]]


def skel_mol(m):
    heavy = [a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum() > 1]
    pos = {i: k for k, i in enumerate(heavy)}
    g = Chem.RWMol()
    for i in heavy:
        a = m.GetAtomWithIdx(i)
        n = Chem.Atom(a.GetAtomicNum())
        if a.GetAtomicNum() == 6:
            n.SetIsotope(1 + a.GetTotalNumHs())
        n.SetNoImplicit(True)
        g.AddAtom(n)
    for b in m.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in pos and j in pos:
            g.AddBond(pos[i], pos[j], Chem.BondType.SINGLE)
    g.UpdatePropertyCache(strict=False)
    return g, heavy


def cip(m):
    m = Chem.Mol(m)
    try:
        Chem.AssignStereochemistryFrom3D(m) if m.GetNumConformers() and m.GetConformer().Is3D() else None
        if m.GetNumConformers():
            Chem.AssignChiralTypesFromBondDirs(m)
        Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    except Exception:
        pass
    return {a.GetIdx(): a.GetProp("_CIPCode") for a in m.GetAtoms() if a.HasProp("_CIPCode")}


def iso(cur, ours):
    """Map curated atom index -> our atom index (heavy atoms), best CIP agreement."""
    gc, hc = skel_mol(cur)
    go, ho = skel_mol(ours)
    if gc.GetNumAtoms() != go.GetNumAtoms():
        return None
    params = Chem.SubstructMatchParameters()
    params.uniquify = False
    params.maxMatches = 5000
    for a in gc.GetAtoms():
        a.SetProp("dummy", "1")
    matches = go.GetSubstructMatches(gc, params)  # tuple index: curated skel atom -> our skel atom
    if not matches:
        return None
    cc, co = cip(cur), cip(ours)
    best = max(matches, key=lambda mt: sum(1 for k, j in enumerate(mt)
                                           if hc[k] in cc and ho[j] in co and cc[hc[k]] == co[ho[j]]))
    return {hc[k]: ho[j] for k, j in enumerate(best)}


rules = []
for recon, mar in pairs_wanted:
    subs, prods = read_rxn(STD / "Curated" / (recon + ".rxn"))
    r = rx_by_id[mar]
    loc = []
    for q, i in enumerate(r["m"]):
        b = bases[mets[i]["b"]]
        m, sk = our_mol(mets[i]["b"])
        loc.append((q, b["n"], m, sk[0], r["c"][q]))
    # curated reactants may correspond to our products if the reaction is written the other way
    best = None
    for orient in (1, -1):
        used, assign, ok = set(), [], True
        for side, mols in ((-1, subs), (1, prods)):
            for cm in mols:
                if all(a.GetAtomicNum() <= 1 for a in cm.GetAtoms()):
                    assign.append(None)
                    continue
                key = skeleton(cm)[0]
                cand = [x for x in loc if x[3] == key and x[4] * side * orient > 0]
                cand.sort(key=lambda x: x[0] in used)  # a metabolite with coefficient > 1 appears several times
                if not cand:
                    ok = False
                    break
                used.add(cand[0][0])
                assign.append(cand[0])
            if not ok:
                break
        if ok:
            best = (orient, assign)
            break
    if best is None:
        print("skip (molecules do not correspond):", recon, mar)
        continue
    orient, assign = best
    allmols = subs + prods
    nsub = len(subs)
    by_num = defaultdict(lambda: ([], []))
    for k, cm in enumerate(allmols):
        x = assign[k]
        if x is None:
            continue
        mp = iso(cm, x[2])
        if mp is None:
            print("skip (no isomorphism):", recon, x[1])
            break
        for a in cm.GetAtoms():
            if a.GetAtomicNum() > 1 and a.GetAtomMapNum() > 0:
                by_num[a.GetAtomMapNum()][0 if k < nsub else 1].append((x[1], mp[a.GetIdx()]))
    else:
        pairs = []
        for num, (rs, ps) in sorted(by_num.items()):
            if not rs or not ps:
                continue
            if orient == -1:
                rs, ps = ps, rs
            for k, p in enumerate(ps):
                s = rs[min(k, len(rs) - 1)] if len(rs) == len(ps) else rs[0]
                if [s[0], s[1], p[0], p[1]] not in pairs:
                    pairs.append([s[0], s[1], p[0], p[1]])
        names_s = sorted({p[0] for p in pairs})
        names_p = sorted({p[2] for p in pairs})
        rules.append({
            "name": f"{r['n']} ({mar}): curated map of Recon 3D {recon} (Preciat Gonzalez et al. 2017, J Cheminform 9:39)",
            "sub": names_s, "prod": names_p, "pairs": pairs,
            "smiles": {n: canon[next(b["id"] for b in bases if b["n"] == n)] for n in names_s + names_p},
        })
        print("rule:", recon, mar, len(pairs), "atom pairs")

json.dump(rules, open(Path(__file__).with_name("curated_rules.json"), "w"), indent=1)
