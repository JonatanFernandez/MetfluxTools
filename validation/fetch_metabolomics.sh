#!/bin/sh
# Download the metabolomics studies used for validation and for the worked examples, from the
# Metabolomics Workbench REST API (https://www.metabolomicsworkbench.org/tools/mw_rest.php).
#
# Each study is fetched as its mwTab text file, which carries the metabolite table, the sample factors
# (the experimental conditions) and the units in one file:
#     /rest/study/analysis_id/<AN>/mwtab/txt
# The /datatable endpoint of the documentation does not exist on the server; mwtab is the one that works.
#
# The files are NOT kept in the repository: run this to get them, then the prepare_*.py scripts turn each
# one into the JSON the validation scripts read. Check each study's own terms before redistributing it.
#
#   sh tools/validation/fetch_metabolomics.sh
set -e
cd "$(dirname "$0")"
mkdir -p metabolomics
cd metabolomics

# study | analysis | what it is
for pair in \
  "ST002076:AN003383" \
  "ST004848:AN008205" \
  "ST004339:AN007246" \
  "ST001609:AN002643" \
  "ST003960:AN006508" \
; do
  ST=$(echo "$pair" | cut -d: -f1)
  AN=$(echo "$pair" | cut -d: -f2)
  if [ -s "$ST.mwtab" ]; then
    echo "$ST already here"
    continue
  fi
  echo "fetching $ST ($AN)…"
  curl -sS -m 300 "https://www.metabolomicsworkbench.org/rest/study/analysis_id/$AN/mwtab/txt" -o "$ST.mwtab"
  # the server answers 200 with "This function does not exist" when a path is wrong
  if [ "$(wc -c < "$ST.mwtab")" -lt 2000 ]; then
    echo "  WARNING: $ST.mwtab is suspiciously small, check it:"
    head -c 200 "$ST.mwtab"
    echo
  fi
done
ls -la *.mwtab
