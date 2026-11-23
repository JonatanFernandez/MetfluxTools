import sys, warnings, logging; warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
sys.path.insert(0, sys.argv[1])
import cobra
m = cobra.io.load_yaml_model(sys.argv[2])
print("PLP (s_1393) in biomass/pseudoreactions:", [r.id for r in m.metabolites.s_1393.reactions if r.id in ("r_4041","r_4598") or "pseudoreaction" in r.name])
print("biomass objective:", [r.id for r in m.reactions if r.objective_coefficient])
def egc(tag, mm):
    with mm:
        for r in mm.reactions:
            if len(r.metabolites) == 1: r.bounds = (0, 1000)      # nothing can enter
        mm.objective = "r_4046"; mm.reactions.r_4046.bounds = (0, 1000)
        v = mm.optimize().objective_value
        # with oxygen and water allowed in (r_0953 consumes them), still no carbon/nitrogen
        for x in ("r_1992", "r_2100"):  # O2 and H2O exchange
            mm.reactions.get_by_id(x).lower_bound = -1000
        v2 = mm.optimize().objective_value
    print(f"[{tag}] max ATP maintenance with no uptake: {v:.6g}; with only O2 and H2O uptake: {v2:.6g}")
print("O2 exchange:", m.reactions.r_1992.reaction, "| H2O exchange:", m.reactions.r_2100.reaction)
egc("current", m)
with m:
    m.remove_reactions(["r_0953"]); egc("r_0953 removed", m)
r = m.reactions.r_0953
r.subtract_metabolites(r.metabolites)
g = m.metabolites.get_by_id
r.add_metabolites({g("s_1394"): -1, g("s_1275"): -1, g("s_0803"): -1, g("s_1392"): 1, g("s_0419"): 1, g("s_0837"): 1})
egc("r_0953 corrected", m)
