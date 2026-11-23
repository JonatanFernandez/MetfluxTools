"""Atom-map a reaction with RXNMapper, for the model editor of MetfluxLive.

    tools/.venv/bin/python tools/map_reaction.py "CC(=O)C(=O)[O-].CC(=O)C(=O)[O-].[H+]>>CC(=O)C(C)O.O=C=O.O=C=O"

Prints the atom-mapped reaction SMILES and RXNMapper's confidence. Paste the SMILES into
"Atom map -> atom-mapped reaction SMILES" of the reaction dialog. Each stoichiometric copy must be written
out (2 pyruvate = two pyruvate molecules). RXNMapper runs locally; nothing is sent anywhere.
"""
import sys


def main():
    if len(sys.argv) < 2 or ">>" not in sys.argv[1]:
        print(__doc__)
        sys.exit(1)
    from rxnmapper import RXNMapper
    res = RXNMapper().get_attention_guided_atom_maps([sys.argv[1]], canonicalize_rxns=False)[0]
    print(res["mapped_rxn"])
    print(f"confidence {res['confidence']:.3f}", file=sys.stderr)


if __name__ == "__main__":
    main()
