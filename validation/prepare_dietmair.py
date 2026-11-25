#!/usr/bin/env python3
"""Intracellular metabolome of HEK293, producer against non-producer -> dietmair2012_intracellular.json.

Source: Table S3 of Dietmair S, Hodson MP, Quek LE, Timmins NE, Gray P, Nielsen LK, "A multi-omics analysis
of recombinant protein production in Hek293 cells", PLoS ONE 7:e43394 (2012), doi 10.1371/journal.pone.0043394,
file pone.0043394.s003.xlsx (CC BY 4.0). Download it from the article page or from
https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3427347/supplementaryFiles (a zip of every item).

The same paper's Table 2 is already in the application as the "hek293" measured exchange set
(js/core/measured.js), so the exometabolome and the intracellular metabolome come from one experiment.

Layout of the sheet: a name column, then 16 numbers — 4 time points x (mean, standard error) for the
producer line, then the same 4 for the non-producer. Values are RELATIVE concentrations, not absolute, so
only the ratio between the two lines means anything. Rows naming a technique ("HPLC", "LC-MS/MS") separate
the blocks and carry no data.

    python3 tools/validation/prepare_dietmair.py <path to pone.0043394.s003.xlsx>
"""
import json, re, sys, zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "tools/validation/dietmair2012_intracellular.json"

# the abbreviations of the sheet -> the names Human-GEM uses. "|" marks a measurement that cannot tell two
# isomers apart, which the importer keeps as one ambiguous node instead of guessing (docs, exometabolomics).
NAMES = {
    "NAD": "NAD+", "NADP": "NADP+", "AMP": "AMP", "ADP": "ADP", "ATP": "ATP", "GDP": "GDP", "GTP": "GTP",
    "UTP": "UTP", "CTP": "CTP", "UDP-GlcA": "UDP-glucuronate", "UDPG": "UDP-glucose",
    "ASP": "aspartate", "GLU": "glutamate", "ASN": "asparagine", "SER": "serine", "GLN": "glutamine",
    "HIS": "histidine", "GLY": "glycine", "THR": "threonine", "ARG": "arginine", "ALA": "alanine",
    "TYR": "tyrosine", "VAL": "valine", "MET": "methionine", "TRP": "tryptophan", "PHE": "phenylalanine",
    "ILE": "isoleucine", "ORN": "ornithine", "LEU": "leucine", "LYS": "lysine", "PRO": "proline",
    "LAC": "L-lactate", "PYR": "pyruvate", "PEP": "PEP",
    "ACO": "cis-aconitate", "ISO": "isocitrate", "CIT": "citrate", "FUM": "fumarate", "SUC": "succinate",
    "OAA": "OAA", "MAL": "malate", "KGA": "AKG",
    "GA3P": "glyceraldehyde-3-phosphate", "DHAP": "DHAP",
    "3PG and 2PG": "3-phospho-D-glycerate|2-phospho-D-glycerate",
    "G1P": "glucose-1-phosphate", "G6P": "glucose-6-phosphate", "F6P": "fructose-6-phosphate",
    "F16DP": "fructose-1,6-bisphosphate", "R5P": "ribose-5-phosphate", "RL5P": "ribulose-5-phosphate",
    "X5P": "D-xylulose-5-phosphate",
}
TECHNIQUES = {"HPLC", "LC-MS/MS", "LC-MS", "GC-MS"}


def col(ref):
    """'AB12' -> 27 (zero-based column index)."""
    letters = re.match(r"([A-Z]+)", ref or "").group(1)
    n = 0
    for ch in letters:
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def sheet_rows(path):
    z = zipfile.ZipFile(path)
    shared = ["".join(re.findall(r"<t[^>]*>(.*?)</t>", si, re.S))
              for si in re.findall(r"<si>(.*?)</si>", z.read("xl/sharedStrings.xml").decode("utf8", "replace"), re.S)]
    xml = z.read("xl/worksheets/sheet1.xml").decode("utf8", "replace")
    for row in re.findall(r"<row[^>]*>(.*?)</row>", xml, re.S):
        cells = {}
        for cell in re.findall(r"<c\b[^>]*>.*?</c>|<c\b[^>]*/>", row, re.S):
            ref = re.search(r'r="([A-Z]+\d+)"', cell)
            val = re.search(r"<v>(.*?)</v>", cell, re.S)
            if not ref or not val:
                continue
            raw = val.group(1)
            is_str = 't="s"' in cell
            cells[col(ref.group(1))] = shared[int(raw)] if is_str and raw.isdigit() and int(raw) < len(shared) else raw
        yield cells


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    src = Path(sys.argv[1])
    rows = list(sheet_rows(src))
    mets, technique, unknown = [], None, []
    for cells in rows:
        label = (cells.get(0) or "").strip()
        if not label:
            continue
        nums = {c: v for c, v in cells.items() if c > 0}
        if label in TECHNIQUES or not nums:
            if label in TECHNIQUES:
                technique = label
            continue
        # 16 values: producer TP1..4 (mean, SE), then non-producer TP1..4
        vals = [float(nums[c]) if c in nums and nums[c] not in ("", None) else None for c in range(1, 17)]
        if sum(v is not None for v in vals) < 8:
            unknown.append(label)
            continue
        pair = lambda i: {"mean": vals[2 * i], "se": vals[2 * i + 1]}
        mets.append({
            "label": label,
            "metabolite": NAMES.get(label),
            "technique": technique,
            "producer": [pair(i) for i in range(4)],
            "nonProducer": [pair(i + 4) for i in range(4)],
        })
    missing = [m["label"] for m in mets if not m["metabolite"]]
    doc = {
        "source": ("Table S3 of Dietmair et al., PLoS ONE 7:e43394 (2012), doi 10.1371/journal.pone.0043394 "
                   "(file pone.0043394.s003.xlsx, CC BY 4.0)"),
        "cells": "HEK293, a stable recombinant-protein producer line and its non-producing parental line",
        "quantity": "relative intracellular concentration (not absolute): only the ratio between the two lines is meaningful",
        "timePoints": 4,
        "spread": "standard error",
        "exometabolome": "the same experiment's exchange rates are the \"hek293\" measured set in js/core/measured.js (their Table 2)",
        "nameMap": "the sheet's abbreviations are mapped to Human-GEM names; \"|\" marks a pair of isomers the measurement cannot separate",
        "metabolites": mets,
    }
    OUT.write_text(json.dumps(doc, indent=1, ensure_ascii=False))
    print(f"{len(mets)} metabolites -> {OUT.relative_to(ROOT)}")
    if missing:
        print("no Human-GEM name mapped for:", ", ".join(missing))
    if unknown:
        print("rows skipped (too few numbers):", ", ".join(unknown))


if __name__ == "__main__":
    main()
