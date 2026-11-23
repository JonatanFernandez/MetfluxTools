# The hand-edited yml must equal "release model minus r_0953" exactly (cobrapy view).
import sys, warnings, logging; warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
sys.path.insert(0, sys.argv[1])
import cobra
a = cobra.io.load_yaml_model(sys.argv[2]); b = cobra.io.load_yaml_model(sys.argv[3])
a.remove_reactions(["r_0953"], remove_orphans=True)
sig = lambda m: ({r.id: (r.name, tuple(sorted((x.id, c) for x, c in r.metabolites.items())), r.bounds, r.gene_reaction_rule, tuple(r.subsystem) if isinstance(r.subsystem, list) else r.subsystem) for r in m.reactions},
                 {x.id: (x.name, x.formula, x.charge, x.compartment) for x in m.metabolites}, sorted(g.id for g in m.genes))
ra, ma, ga = sig(a); rb, mb, gb = sig(b)
print("reactions:", len(ra), len(rb), "identical" if ra == rb else "DIFFER")
print("metabolites:", len(ma), len(mb), "identical" if ma == mb else "DIFFER")
print("genes:", len(ga), len(gb), "identical" if ga == gb else "DIFFER")
print("growth (canonical medium):", a.optimize().objective_value, b.optimize().objective_value)
