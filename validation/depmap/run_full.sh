#!/bin/sh
# All DepMap 26Q3 lines with a cell type in the catalogue, objective "growth without the cofactor pool" and medium
# "standard + serum lipids" (CORR=nocof,serum), 4 shards in parallel, block then complex.
# Then: TAG=nocof+serum node tools/validation/depmap_evaluate.mjs
cd "$(dirname "$0")/../../.."
for v in block complex; do
  for k in 0 1 2 3; do
    DEPMAP=depmap26q3_metabolic.json CORR=nocof,serum STEP=1 nice -n 10 node tools/validation/depmap_essentiality.mjs $v $k 4 > tools/validation/depmap/logfull_${v}_$k.txt 2>&1 &
  done
  wait
done
echo ALL DONE
