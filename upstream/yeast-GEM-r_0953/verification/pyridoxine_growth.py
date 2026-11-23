import sys, warnings, logging; warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
sys.path.insert(0, sys.argv[1] + "/ygdeps")
import cobra
m = cobra.io.load_yaml_model(sys.argv[1] + "/yg/yeast_develop.yml")
ex = [r for r in m.reactions if len(r.metabolites) == 1]
print("open uptakes (canonical):", sorted((r.id, r.name, r.lower_bound) for r in ex if r.lower_bound < 0)[:40])
def run(tag, mm):
    with mm:
        mm.reactions.r_2028.lower_bound = -1
        s = cobra.flux_analysis.pfba(mm)
        mu = s.fluxes["r_2111"]
        print(f"\n[{tag}] growth with pyridoxine {mu:.4f}")
        top = s.fluxes[[r.id for r in ex]]; top = top[top.abs() > 1e-6].sort_values()
        print(" exchanges:", {k: round(v,3) for k, v in top.items()})
        b6 = ["r_0953","r_0954","r_0955","r_0956","r_2025","r_2026","r_2027","r_4330","r_2029","r_4211","r_4212"]
        print(" B6:", {k: round(s.fluxes[k],3) for k in b6 if k in mm.reactions})
        return s
s0 = run("develop", m)
with m:
    m.remove_reactions(["r_0953"]); run("r_0953 removed", m)
with m:
    for r in ["r_0953","r_0956"]: pass
    m.reactions.r_0953.bounds=(0,0); m.reactions.r_0956.bounds=(0,0); run("r_0953 and r_0956 off", m)
