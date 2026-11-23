"""Print the C/N transitions that differ between MetfluxLive and the curated Recon 3D map of one reaction."""
import sys
sys.argv = [sys.argv[0], sys.argv[1]] + sys.argv[2:]
recon, hg = sys.argv[2], sys.argv[3]
src = open(__file__.replace("diff_one.py", "atommaps_vs_recon3d.py")).read().split("# ---------------------------------------------------------------- Recon3D id")[0]
exec(compile(src, "v", "exec"))
ct, ck = transitions_rxn(*read_rxn(STD / "Curated" / (recon + ".rxn")))
ot, ok, kind = transitions_ours(hg)
d = comparable(ck, ok)
ot = ot if d == 1 else flip(ot)
r = rx_by_id[hg]
names = {}
for i in r["m"]:
    b = bases[mets[i]["b"]]; m, sk = our_mol(mets[i]["b"])
    names[sk[0]] = b["n"]
def show(t):
    out = []
    for pair in sorted(t, key=str):
        out.append("  " + "  <->  ".join(f"{x[0]}:{names.get(x[1], x[1][:30])}#{x[2]}" for x in sorted(pair)))
    return "\n".join(out)
print(hg, r["n"], "direction", d)
print("only MetfluxLive:\n" + show(ot - ct)); print("only curated:\n" + show(ct - ot))
# molecules with their symmetry class as atom-map number (class+1), to read the transitions above
shown = set()
for i in r["m"]:
    b = bases[mets[i]["b"]]; m, sk = our_mol(mets[i]["b"])
    if b["n"] in shown or not any(a.GetAtomicNum() in (6, 7) for a in m.GetAtoms()):
        continue
    shown.add(b["n"])
    mm = Chem.Mol(m)
    for a in mm.GetAtoms():
        if a.GetAtomicNum() in (6, 7):
            a.SetAtomMapNum(sk[1][a.GetIdx()] + 1)
    print(f"  {b['n']}: {Chem.MolToSmiles(mm)}")
