import warnings; warnings.filterwarnings("ignore")
from equilibrator_api import ComponentContribution, Q_
import importlib.metadata as md
print("equilibrator-api", md.version("equilibrator-api"), "| equilibrator-cache", md.version("equilibrator-cache"), "| component-contribution", md.version("component-contribution"))
cc = ComponentContribution()
R = {
 "r_0953 current (PL + NH4 + 2 H2O + 0.5 O2 -> PM + 2 H2O2)": "kegg:C00250 + kegg:C00014 + 2 kegg:C00001 + 0.5 kegg:C00007 = kegg:C00534 + 2 kegg:C00027",
 "r_0953 proposed (PM + O2 + H2O -> PL + NH4 + H2O2)": "kegg:C00534 + kegg:C00007 + kegg:C00001 = kegg:C00250 + kegg:C00014 + kegg:C00027",
 "r_0954 PMP oxidase (PMP + O2 + H2O -> PLP + NH4 + H2O2)": "kegg:C00647 + kegg:C00007 + kegg:C00001 = kegg:C00018 + kegg:C00014 + kegg:C00027",
 "r_0955 PNP oxidase (PNP + O2 -> PLP + H2O2)": "kegg:C00627 + kegg:C00007 = kegg:C00018 + kegg:C00027",
 "r_0956 PN oxidase (PN + O2 -> PL + H2O2)": "kegg:C00314 + kegg:C00007 = kegg:C00250 + kegg:C00027",
 "H2O + 0.5 O2 -> H2O2": "kegg:C00001 + 0.5 kegg:C00007 = kegg:C00027",
}
for T in (298.15, 303.15):
  for pH in (7.2, 6.7):
    cc.p_h = Q_(pH); cc.ionic_strength = Q_("0.25M"); cc.p_mg = Q_(3.0); cc.temperature = Q_(f"{T}K")
    print(f"\n-- T={T} K, pH={pH}, I=0.25 M, pMg=3")
    for k, s in R.items():
        r = cc.parse_reaction_formula(s)
        assert r.is_balanced(ignore_atoms=("H",)) , k
        bal = r.is_balanced()
        g0 = cc.standard_dg_prime(r); gm = cc.physiological_dg_prime(r)
        print(f"{k:62s} dG'0 = {g0.value.m_as('kJ/mol'):7.1f} ± {g0.error.m_as('kJ/mol'):4.1f}   dG'm = {gm.value.m_as('kJ/mol'):7.1f} ± {gm.error.m_as('kJ/mol'):4.1f} kJ/mol  balanced(with H)={bal}")
