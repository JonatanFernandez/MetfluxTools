# yeast-GEM's own benchmarks (essential_genes, growth) and B6 checks for: develop, r_0953 removed, r_0953 corrected
import sys, os, warnings, logging, time
SP = sys.argv[1]; sys.path[:0] = [SP + "/ygdeps", SP + "/yg/repo/code/python"]
warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
os.environ["YEAST_GEM_PATH"] = SP + "/yg/repo"
import cobra
from cobra.flux_analysis import flux_variability_analysis as fva
from yeastgem.model_tests import essential_genes, growth
print("cobrapy", cobra.__version__, "| solver", cobra.Configuration().solver.__name__, flush=True)
def variant(name):
    m = cobra.io.load_yaml_model(SP + "/yg/yeast_develop.yml")   # develop as published (same content as 9.1.1)
    if name == "removed":
        m.remove_reactions(["r_0953"], remove_orphans=True)
    elif name == "corrected":
        r = m.reactions.r_0953; g = m.metabolites.get_by_id
        r.subtract_metabolites(r.metabolites)
        r.add_metabolites({g("s_1394"): -1, g("s_1275"): -1, g("s_0803"): -1, g("s_1392"): 1, g("s_0419"): 1, g("s_0837"): 1})
    return m
res = {}
B6 = ["r_0953", "r_0954", "r_0955", "r_0956", "r_2025", "r_2026", "r_2027", "r_4330"]
for name in ("develop", "removed", "corrected"):
    t0 = time.time(); m = variant(name)
    mu = m.optimize().objective_value
    with m:
        m.reactions.r_2028.lower_bound = -1          # pyridoxine uptake allowed
        mu_pn = m.optimize().objective_value
        f = fva(m, [x for x in B6 if x in m.reactions], fraction_of_optimum=0.0, processes=1)
    eg = essential_genes(m, processes=1)
    r2 = growth(m)
    res[name] = (mu, mu_pn, eg, r2)
    print(f"\n== {name}  ({time.time()-t0:.0f} s)\n growth canonical medium {mu:.6f} | + pyridoxine {mu_pn:.6f}")
    print(" FVA with pyridoxine (min, max):", {i: (round(f.loc[i,'minimum'],3), round(f.loc[i,'maximum'],3)) for i in f.index})
    print(f" essential genes: accuracy {eg.accuracy:.4f}  MCC {eg.mcc:.4f}  TP {len(eg.tp)} TN {len(eg.tn)} FP {len(eg.fp)} FN {len(eg.fn)}")
    print(f" growth benchmark (Tobias 2013 chemostats) R2 {r2:.4f}", flush=True)
d = res["develop"][2]
for name in ("removed", "corrected"):
    e = res[name][2]
    same = all(set(getattr(d, k)) == set(getattr(e, k)) for k in ("tp", "tn", "fp", "fn"))
    print(f"{name}: per-gene essentiality calls identical to develop: {same}")
