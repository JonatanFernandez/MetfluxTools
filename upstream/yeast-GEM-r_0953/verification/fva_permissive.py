import sys, warnings, logging; warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
sys.path.insert(0, sys.argv[1])
import cobra
from cobra.flux_analysis import flux_variability_analysis as fva
m = cobra.io.load_yaml_model(sys.argv[2])
ex = [r for r in m.reactions if len(r.metabolites) == 1]
print(len(ex), "boundary reactions")
ids = ["r_0953","r_0954","r_0955","r_0956","r_2025","r_2026","r_2027","r_4330","r_4211","r_4212"]
with m:
    for r in ex: r.bounds = (-1000, 1000)   # every boundary open both ways: the most permissive medium
    # also release every reaction that is switched off by default (lb=ub=0), e.g. glycine cleavage
    off = [r for r in m.reactions if r.lower_bound == 0 and r.upper_bound == 0]
    for r in off: r.upper_bound = 1000
    print("reactions switched off by default and opened:", [r.id for r in off])
    f = fva(m, ids, fraction_of_optimum=0.0, loopless=False)
    print(f.round(6).to_string())
    print("s_1394 (pyridoxamine[c]) reactions:", [(r.id, r.metabolites[m.metabolites.s_1394]) for r in m.metabolites.s_1394.reactions])
    print("s_1392 (pyridoxal[c]) reactions:", [(r.id, r.metabolites[m.metabolites.s_1392]) for r in m.metabolites.s_1392.reactions])
