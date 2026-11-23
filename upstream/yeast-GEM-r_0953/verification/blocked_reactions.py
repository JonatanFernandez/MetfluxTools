import sys, warnings, logging; warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
sys.path.insert(0, sys.argv[1] + "/ygdeps")
import cobra
from cobra.flux_analysis import find_blocked_reactions
def load(v):
    m = cobra.io.load_yaml_model(sys.argv[1] + "/yg/yeast_develop.yml")
    if v == "removed": m.remove_reactions(["r_0953"])
    if v == "corrected":
        r = m.reactions.r_0953; g = m.metabolites.get_by_id; r.subtract_metabolites(r.metabolites)
        r.add_metabolites({g("s_1394"): -1, g("s_1275"): -1, g("s_0803"): -1, g("s_1392"): 1, g("s_0419"): 1, g("s_0837"): 1})
    return m
res = {}
for v in ("develop", "removed", "corrected"):
    m = load(v)
    res[v] = set(find_blocked_reactions(m, open_exchanges=True, processes=1))
    print(v, "blocked (all exchanges open):", len(res[v]), flush=True)
for v in ("removed", "corrected"):
    nb = sorted(res[v] - res["develop"]); ub = sorted(res["develop"] - res[v])
    m = load("develop")
    print(f"\n{v}: newly blocked {len(nb)}:", [(i, m.reactions.get_by_id(i).name) for i in nb])
    print(f"{v}: newly unblocked {len(ub)}:", ub)
