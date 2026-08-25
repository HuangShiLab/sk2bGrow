#!/usr/bin/env bash
# Arm E: Pilea's FracMinHash sketch through sk2bGrow's estimator, unmodified.
# Completes the 2x2 of sketch x estimator:
#     A = anchors + V-fit      E = FracMinHash + V-fit
#     B = anchors + rank       C = FracMinHash + rank  (Pilea itself)
set -uo pipefail
cd "$(dirname "$0")"
ROOT=/Users/shihuang/Downloads/sk2bGrow
PILEA_PY=${PILEA_PY:-../pilea_env/bin/python}

$PILEA_PY armE_counts.py -r ../ecoli.fna -o counts sub/*.fq || exit 1

mkdir -p out
for f in counts/*.counts.tsv; do
  s=$(basename "$f" .counts.tsv); o="out/E_$s"
  [ -f "$o/output.tsv" ] && continue
  PYTHONPATH=$ROOT/python python3 -m sk2bgrow.cli profile "$f" \
      --db db --output "$o" --min-coverage 0 >/dev/null 2>&1 || echo "FAIL E $s"
done
echo "arm E: $(ls -d out/E_*/ 2>/dev/null | wc -l) cells"
