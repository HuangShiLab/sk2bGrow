#!/usr/bin/env bash
# How fragmented is too fragmented? The same genome cut into 2, 5, 10, 20, 50 and
# 100 contigs, at one depth, so the collapse has a dose-response curve rather
# than a single point -- which is what any contig-count QC threshold needs.
set -uo pipefail
cd "$(dirname "$0")"
ROOT=${ROOT:-$HOME/Downloads/sk2bGrow}
BIN=$ROOT/target/release/sk2bgrow
WORK=${WORK:?set WORK}
SUB=${SUB:?set SUB}
GENOME=${GENOME:?set GENOME}
DEPTH=${DEPTH:-10}
COUNTS=${COUNTS:-"2 5 10 20 50 100"}
export PYTHONPATH=$ROOT/python
mkdir -p "$WORK/sweep"

for n in $COUNTS; do
  [ -f "$WORK/sweep/n$n.fna" ] || \
      python3 fragment.py "$GENOME" -n "$n" --seed 0 -o "$WORK/sweep/n$n.fna" | sed 's/^/  /'
  [ -f "$WORK/sweep/db_n$n/manifest.json" ] || \
      $BIN index "$WORK/sweep/n$n.fna" -o "$WORK/sweep/db_n$n" --quiet
done

for f in "$SUB"/*_${DEPTH}x.fq; do
  s=$(basename "$f" .fq)
  for n in $COUNTS; do
    o="$WORK/sweep/out/n${n}_${s}"
    [ -f "$o/output.tsv" ] && continue
    $BIN profile "$f" -d "$WORK/sweep/db_n$n" -o "$o" --quiet --python python3 \
        >/dev/null 2>&1 || echo "FAIL n=$n $s"
  done
done
echo "sweep: $(ls -d "$WORK"/sweep/out/*/ 2>/dev/null | wc -l) runs"
