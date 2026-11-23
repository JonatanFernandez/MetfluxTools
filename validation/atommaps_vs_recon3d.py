#!/usr/bin/env python3
"""Validation 1: MetfluxLive atom maps (Human-GEM) against the manually curated atom mappings of
Preciat Gonzalez et al. 2017 (J Cheminform 9:39, Additional file 3, standardisedPredictions/Curated,
512 Recon 3D reactions). The same comparison is run on the standardised predictions of the published
algorithms in that file (DREAM, RDT, CLCA, ...), so every method is scored with one criterion.

Criterion (per reaction, carbon and nitrogen atoms only, i.e. the atoms a 13C/15N tracer follows):
  every C/N atom of a reactant or product is reduced to (molecule, symmetry class). Molecules are
  compared as heavy-atom graphs labelled by element and, for carbon, the number of attached H
  (invariant to protonation state), with bond orders ignored. A reaction map becomes the SET of
  transitions {(molecule, class) <-> (molecule, class)}. Two maps agree when the sets are equal.
  Reactions are comparable when both sides contain the same C/N-bearing molecules.
  Equivalent atoms are handled by the symmetry classes; the curated files already give equivalent
  atoms a common map number.

Usage: tools/.venv/bin/python tools/validation/atommaps_vs_recon3d.py <dir with standardisedPredictions>
"""
import csv, io, json, sys, zipfile
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[2]
STD = Path(sys.argv[1])
EFFECTIVE = "--stored" not in sys.argv  # default: the transitions the tracer actually uses


def skeleton(m):
    """(key, class per atom index) of a sanitised RDKit molecule; H atoms dropped."""
    heavy = [a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum() > 1]
    pos = {i: k for k, i in enumerate(heavy)}
    g = Chem.RWMol()
    for i in heavy:
        a = m.GetAtomWithIdx(i)
        n = Chem.Atom(a.GetAtomicNum() if a.GetAtomicNum() > 0 else 0)
        if a.GetAtomicNum() == 6:
            n.SetIsotope(1 + a.GetTotalNumHs())
        n.SetNoImplicit(True)
        g.AddAtom(n)
    for b in m.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in pos and j in pos:
            g.AddBond(pos[i], pos[j], Chem.BondType.SINGLE)
    g.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(g)
    ranks = list(Chem.CanonicalRankAtoms(g, breakTies=False))
    key = Chem.MolToSmiles(g)
    return key, {i: ranks[pos[i]] for i in heavy}


def sanitise(m):
    try:
        Chem.SanitizeMol(m)
    except Exception:
        m.UpdatePropertyCache(strict=False)
    return m


def read_rxn(path):
    """Reactant and product molecules of an MDL RXN file, each with atom-map numbers."""
    txt = path.read_text(errors="replace")
    blocks = txt.split("$MOL")
    head = blocks[0].splitlines()
    nr, npr = int(head[4][:3]), int(head[4][3:6])
    mols = []
    for blk in blocks[1:]:
        m = Chem.MolFromMolBlock(blk.split("\n", 1)[1], sanitize=False, removeHs=False)
        if m is None:
            return None
        mols.append(sanitise(m))
    return mols[:nr], mols[nr:nr + npr]


def transitions_rxn(subs, prods):
    """Transition set and molecule keys of a mapped RXN (atoms sharing a map number correspond)."""
    by_num = defaultdict(lambda: ([], []))
    keys = [[], []]
    for side, mols in ((0, subs), (1, prods)):
        for m in mols:
            key, cls = skeleton(m)
            keys[side].append((key, m))
            for a in m.GetAtoms():
                if a.GetAtomicNum() in (6, 7) and a.GetAtomMapNum() > 0:
                    by_num[a.GetAtomMapNum()][side].append((key, cls[a.GetIdx()]))
    t = set()
    for s_list, p_list in by_num.values():
        for s in s_list:
            for p in p_list:
                t.add(frozenset([("s",) + s, ("p",) + p]))
    return t, keys


def cn_keys(keys):
    return Counter(k for k, m in keys if any(a.GetAtomicNum() in (6, 7) for a in m.GetAtoms()))


def comparable(ka, kb):
    """Same C/N-bearing molecules on each side (either direction). Returns +1, -1 or 0."""
    if cn_keys(ka[0]).keys() == cn_keys(kb[0]).keys() and cn_keys(ka[1]).keys() == cn_keys(kb[1]).keys():
        return 1
    if cn_keys(ka[0]).keys() == cn_keys(kb[1]).keys() and cn_keys(ka[1]).keys() == cn_keys(kb[0]).keys():
        return -1
    return 0


def flip(t):
    out = set()
    for pair in t:
        out.add(frozenset((("p" if x[0] == "s" else "s"),) + x[1:] for x in pair))
    return out


# ---------------------------------------------------------------- MetfluxLive side
data = ROOT / "data" / "human-gem"
model = json.loads((data / "model.json").read_text())
amap = json.loads((data / "atommaps.json").read_text())
canon = json.loads((ROOT / "tools" / "cache" / "human-gem" / "canon.json").read_text())
bases, mets = model["bases"], model["mets"]
rx_by_id = {r["id"]: r for r in model["rxns"]}
mol_cache = {}


def our_mol(bi):
    bid = bases[bi]["id"]
    if bid not in mol_cache:
        smi = canon.get(bid)
        m = Chem.MolFromSmiles(smi) if smi else None
        mol_cache[bid] = (m, skeleton(m) if m is not None else None)
    return mol_cache[bid]


def transitions_ours(rid):
    r = rx_by_id[rid]
    loc = [our_mol(mets[i]["b"]) for i in r["m"]]
    keys = [[], []]
    for q, (m, sk) in enumerate(loc):
        if m is None:
            return None, None, "no structure"
        if abs(r["c"][q]) > 0:
            keys[0 if r["c"][q] < 0 else 1].append((sk[0], m))
    e = amap.get(rid)
    if e is None:
        return None, keys, "no map"
    t = set()
    used_s, used_p = defaultdict(set), defaultdict(set)
    for seg in e["s"]:
        sl, pl, flat = seg[0], seg[1], seg[2]
        (ms, (ks, cs)), (mp, (kp, cp)) = loc[sl], loc[pl]
        pairs = [(i, i) for i in range(ms.GetNumAtoms())] if flat == -1 else list(zip(flat[0::2], flat[1::2]))
        for i, j in pairs:
            used_s[sl].add(i); used_p[pl].add(j)
            if ms.GetAtomWithIdx(i).GetAtomicNum() in (6, 7):
                t.add(frozenset([("s", ks, cs[i]), ("p", kp, cp[j])]))
    if EFFECTIVE:
        # what js/core/tracer.js does with atoms left without a transition: each such product atom
        # receives the mean label of the unmapped substrate atoms of the same element (unmapped CoA
        # atoms of a product stay unlabelled instead)
        left = defaultdict(list)
        for q, (m, (k, c)) in enumerate(loc):
            if r["c"][q] < 0:
                for a in m.GetAtoms():
                    if a.GetAtomicNum() in (6, 7) and a.GetIdx() not in used_s[q]:
                        left[a.GetAtomicNum()].append(("s", k, c[a.GetIdx()]))
        for q, (m, (k, c)) in enumerate(loc):
            if r["c"][q] > 0 and bases[mets[r["m"][q]]["b"]]["n"] != "CoA":
                for a in m.GetAtoms():
                    if a.GetAtomicNum() in (6, 7) and a.GetIdx() not in used_p[q]:
                        for src in left[a.GetAtomicNum()]:
                            t.add(frozenset([src, ("p", k, c[a.GetIdx()])]))
    return t, keys, e["k"]


# ---------------------------------------------------------------- Recon3D id -> Human-GEM reactions
zf = zipfile.ZipFile(ROOT / "Human-GEM-2.1.1.zip")
name = next(n for n in zf.namelist() if n.endswith("model/reactions.tsv"))
rows = list(csv.DictReader(io.TextIOWrapper(zf.open(name), encoding="utf-8"), delimiter="\t"))
to_hg = defaultdict(list)
for row in rows:
    for x in (row["rxnRecon3DID"] or "").split(";"):
        x = x.strip()
        if x and row["rxns"] in rx_by_id:
            to_hg[x].append(row["rxns"])
    hmr = (row.get("rxnHMR2ID") or "").strip()
    if hmr and row["rxns"] in rx_by_id:
        to_hg.setdefault(hmr, []).append(row["rxns"]) if row["rxns"] not in to_hg.get(hmr, []) else None


def hg_for(recon_id):
    for cand in (recon_id, recon_id + "f", recon_id + "r", recon_id + "i"):
        if cand in to_hg:
            return to_hg[cand]
    return []


# ---------------------------------------------------------------- compare
curated = {p.stem: read_rxn(p) for p in sorted((STD / "Curated").glob("*.rxn"))}
methods = [d.name for d in sorted(STD.iterdir()) if d.is_dir() and d.name != "Curated"]
cur_t = {}
for rid, rx in curated.items():
    if rx is not None:
        cur_t[rid] = transitions_rxn(*rx)

stats = Counter()
detail = []
ours_ok = {}
for rid, (ct, ckeys) in cur_t.items():
    if not any(t for t in ct):
        stats["curated without C/N transitions"] += 1
        continue
    hg = hg_for(rid)
    if not hg:
        stats["not in Human-GEM"] += 1
        detail.append((rid, "", "not in Human-GEM"))
        continue
    verdict, used = None, ""
    for h in hg:
        ot, okeys, kind = transitions_ours(h)
        if okeys is None:
            v = "no structure"
        elif ot is None:
            v = "no map"
        else:
            d = comparable(ckeys, okeys)
            if d == 0:
                v = "not comparable"
            else:
                v = "agree" if (ot if d == 1 else flip(ot)) == ct else "differ"
        rank = ["agree", "differ", "not comparable", "no map", "no structure"]
        if verdict is None or rank.index(v) < rank.index(verdict):
            verdict, used = v, h + ":" + (kind if okeys is not None and ot is not None else "")
    stats[verdict] += 1
    detail.append((rid, used, verdict))
    if verdict in ("agree", "differ"):
        ours_ok[rid] = verdict == "agree"

print("MetfluxLive vs curated (Recon 3D reactions):", dict(stats))
n = len(ours_ok)
print(f"  comparable: {n}; agree {sum(ours_ok.values())} ({100*sum(ours_ok.values())/max(n,1):.1f} %)")
by_kind = defaultdict(Counter)
for rid, used, v in detail:
    if v in ("agree", "differ"):
        by_kind[used.split(":")[1]][v] += 1
print("  by MetfluxLive method:", {k: dict(v) for k, v in by_kind.items()})

# same criterion, same reactions, for the published algorithms
res = {"MetfluxLive": (sum(ours_ok.values()), n)}
for meth in methods:
    ok = tot = 0
    for rid in ours_ok:
        p = STD / meth / (rid + ".rxn")
        if not p.exists():
            continue
        rx = read_rxn(p)
        if rx is None:
            continue
        mt, mkeys = transitions_rxn(*rx)
        d = comparable(cur_t[rid][1], mkeys)
        if d == 0:
            continue
        tot += 1
        ok += (mt if d == 1 else flip(mt)) == cur_t[rid][0]
    res[meth] = (ok, tot)
    # on the subset this method covers, how does MetfluxLive do?
for meth, (ok, tot) in res.items():
    print(f"  {meth:12s} {ok}/{tot} = {100*ok/max(tot,1):.1f} %")

out = ROOT / "tools" / "validation" / ("atommaps_vs_recon3d.tsv" if EFFECTIVE else "atommaps_vs_recon3d_stored.tsv")
with out.open("w") as f:
    f.write("recon3d_rxn\thumangem_rxn:method\tverdict\n")
    for row in detail:
        f.write("\t".join(row) + "\n")
print("details:", out)

# all methods on the reactions every one of them maps (and MetfluxLive, on the same set)
per = {}
for meth in methods:
    per[meth] = {}
    for rid in ours_ok:
        p = STD / meth / (rid + ".rxn")
        rx = read_rxn(p) if p.exists() else None
        if rx is None:
            continue
        mt, mkeys = transitions_rxn(*rx)
        d = comparable(cur_t[rid][1], mkeys)
        if d:
            per[meth][rid] = (mt if d == 1 else flip(mt)) == cur_t[rid][0]
common = set(ours_ok).intersection(*[set(v) for v in per.values()])
print(f"common subset: {len(common)} reactions")
print(f"  {'MetfluxLive':12s} {sum(ours_ok[r] for r in common)}/{len(common)}")
for meth, v in per.items():
    print(f"  {meth:12s} {sum(v[r] for r in common)}/{len(common)}")
