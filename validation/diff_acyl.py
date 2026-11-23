"""Like diff_one.py but hides atoms of the CoA moiety and prints classes as annotated in the SMILES (class+1)."""
import sys
recon, hg = sys.argv[2], sys.argv[3]
src = open(__file__.replace("diff_acyl.py", "atommaps_vs_recon3d.py")).read().split("# ---------------------------------------------------------------- Recon3D id")[0]
exec(compile(src, "v", "exec"))
coa = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(~O)OP(=O)(~O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(~O)~O")
ct, ck = transitions_rxn(*read_rxn(STD / "Curated" / (recon + ".rxn")))
ot, ok, kind = transitions_ours(hg)
d = comparable(ck, ok); ot = ot if d == 1 else flip(ot)
r = rx_by_id[hg]; names = {}; coacls = defaultdict(set)
for i in r["m"]:
    b = bases[mets[i]["b"]]; m, sk = our_mol(mets[i]["b"]); names[sk[0]] = b["n"]
    for match in m.GetSubstructMatches(coa):
        coacls[sk[0]] |= {sk[1][x] for x in match if m.GetAtomWithIdx(x).GetAtomicNum() in (6, 7)}
    if b["n"] != "CoA":
        mm = Chem.Mol(m)
        for a in mm.GetAtoms():
            if a.GetAtomicNum() in (6, 7) and sk[1][a.GetIdx()] not in coacls[sk[0]]: a.SetAtomMapNum(sk[1][a.GetIdx()] + 1)
        smi = Chem.MolToSmiles(mm); print(f"  {b['n']}: {smi[smi.find('S'):] if 'CoA' in b['n'] else smi}")
def fmt(t):
    out = []
    for pair in t:
        p = next(x for x in pair if x[0] == "p") if any(x[0] == "p" for x in pair) else None
        s = next(x for x in pair if x[0] == "s") if any(x[0] == "s" for x in pair) else None
        if p is None or s is None: continue
        if p[2] in coacls[p[1]] and s[2] in coacls[s[1]]: continue
        out.append(f"{names.get(p[1],'?')}:{p[2]+1} <- {names.get(s[1],'?')}:{s[2]+1}")
    return sorted(out)
print("ours   :", fmt(ot - ct)); print("curated:", fmt(ct - ot))
coa_o = sorted({(names.get(p[1]), names.get(s[1])) for pair in ot for p in pair for s in pair if p[0]=='p' and s[0]=='s' and p[2] in coacls[p[1]] and s[2] in coacls[s[1]]})
coa_c = sorted({(names.get(p[1]), names.get(s[1])) for pair in ct for p in pair for s in pair if p[0]=='p' and s[0]=='s' and p[2] in coacls[p[1]] and s[2] in coacls[s[1]]})
print("CoA moiety ours   :", coa_o); print("CoA moiety curated:", coa_c)
