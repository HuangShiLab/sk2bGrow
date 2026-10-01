#!/usr/bin/env bash
# R3: profile the simulated grid against every reference condition.
# References: complete, correctly-ordered (ord10, ord100), scrambled (2..100).
# Samples: log2PTR in {0, 0.5, 1.0, 1.5} x depth in {1, 5, 10}.
set -uo pipefail
cd "$(dirname "$0")"
ROOT=${ROOT:-$HOME/Downloads/sk2bGrow}
BIN=$ROOT/target/release/sk2bgrow
WORK=${WORK:-$ROOT/benches/work/r3}
PY=${PYTHON:-$ROOT/.venv/bin/python}
export PYTHONPATH=$ROOT/python

mkdir -p "$WORK/out"
for ref in "$WORK"/refs/*.fna "$WORK"/refsR/*.fna; do
  name=$(basename "$ref" .fna)
  [ -f "$WORK/db_$name/manifest.json" ] || \
      $BIN index "$ref" -o "$WORK/db_$name" --quiet || echo "FAIL index $name"
done

for f in "$WORK"/reads/*.fq "$WORK"/readsR/*.fq; do
  s=$(basename "$f" .fq)
  case "$f" in *readsR*) refs="$WORK"/refsR/*.fna ;; *) refs="$WORK"/refs/*.fna ;; esac
  for ref in $refs; do
    name=$(basename "$ref" .fna)
    o="$WORK/out/${name}_${s}"
    [ -f "$o/output.tsv" ] && continue
    # R3 deliberately evaluates the coordinate fit on broken coordinates. Under
    # the current fail-closed auto rule this choice must be explicit.
    $BIN profile "$f" -d "$WORK/db_$name" -o "$o" --quiet --python "$PY" \
        --method v_shape --min-coverage 0 \
        >/dev/null 2>&1 || echo "FAIL $name $s"
  done
done
echo "R3 runs: $(ls -d "$WORK"/out/*/ 2>/dev/null | wc -l)"
