#!/usr/bin/env bash
# Zheng et al. 2020 E. coli benchmark.
#   arm A: sk2bGrow default        (16-enzyme anchors, adaptive windows, V-shape fit)
#   arm B: sk2bGrow Pilea-parity   (25 kb windows, sorted+RANSAC)  -> isolates the estimator
#   arm C: Pilea                   (its own FracMinHash sketch + estimator)
set -uo pipefail
cd "$(dirname "$0")"
ROOT=/Users/shihuang/Downloads/sk2bGrow
BIN=$ROOT/target/release/sk2bgrow
GLEN=4641652
COVS="${COVS:-0.5 1 2 5 10}"
mkdir -p sub out

# --- index once -------------------------------------------------------------
if [ ! -f db/manifest.json ]; then
  $BIN index ../ecoli.fna -o db --enzymes all --quiet
fi

for f in fq/*.fq; do
  s=$(basename "$f" .fq); medium="${s%%.*}"
  have=$(( $(wc -l < "$f") / 4 ))
  for cov in $COVS; do
    n=$(python3 -c "print(int($cov*$GLEN/150))")
    [ "$n" -gt "$have" ] && continue
    sfq="sub/${medium}_${cov}x.fq"
    [ -s "$sfq" ] || head -n $((n*4)) "$f" > "$sfq"

    # arm A
    o="out/A_${medium}_${cov}x"
    if [ ! -f "$o/output.tsv" ]; then
      PYTHONPATH=$ROOT/python $BIN profile "$sfq" -d db -o "$o" --quiet --python python3 \
        >/dev/null 2>"$o.log" || echo "FAIL A $medium $cov"
    fi
    # arm B: Pilea-parity windowing + estimator
    o="out/B_${medium}_${cov}x"
    if [ ! -f "$o/output.tsv" ]; then
      $BIN profile "$sfq" -d db -o "$o" --quiet --no-stats --windowing bp --window-bp 25000 >/dev/null 2>&1
      PYTHONPATH=$ROOT/python python3 -m sk2bgrow.cli profile "$o"/*.counts.tsv \
        --db db --output "$o" --use-rust-windows --method sorted --min-coverage 0 \
        >/dev/null 2>"$o.log" || echo "FAIL B $medium $cov"
    fi
  done
done
echo "sk2bGrow arms done: $(ls -d out/A_* 2>/dev/null | wc -l) A runs, $(ls -d out/B_* 2>/dev/null | wc -l) B runs"
