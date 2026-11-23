#!/bin/sh
# Every 6th DepMap line, 4 shards in parallel, block then complex. Logs in tools/validation/depmap/log_*.txt
cd "$(dirname "$0")/../../.."
for v in block complex; do
  for k in 0 1 2 3; do
    STEP=6 nice -n 10 node tools/validation/depmap_essentiality.mjs $v $k 4 > tools/validation/depmap/log_${v}_$k.txt 2>&1 &
  done
  wait
done
echo ALL DONE
