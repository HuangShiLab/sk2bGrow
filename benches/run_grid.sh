#!/usr/bin/env bash
# Multi-strain grid, following Pilea's Fig 3 design at laptop scale.
# Reads are generated, consumed by both tools, then deleted -- never more than
# one sample on disk at a time.
cd "$(dirname "$0")"
ROOT=/Users/shihuang/Downloads/sk2bGrow
BIN=$ROOT/target/release/sk2bgrow
P=../pilea_env/bin/pilea
STRAINS="${STRAINS:-4 8 16}"; COVS="${COVS:-1 2 4 8}"; REPS="${REPS:-2}"
mkdir -p res
for s in $STRAINS; do for c in $COVS; do for r in $(seq 1 $REPS); do
  tag="s${s}_c${c}_r${r}"
  [ -f "res/${tag}.done" ] && continue
  fq="/tmp/sim_${tag}.fq"
  python3 simulate.py --n-strains $s --coverage $c --seed $((s*1000+c*10+r)) \
      --out "$fq" --truth "res/${tag}.truth" >/dev/null

  /usr/bin/time -l $BIN profile "$fq" -d db -o "res/A_${tag}" --quiet --no-stats \
      2>"res/A_${tag}.time" >/dev/null
  PYTHONPATH=$ROOT/python python3 -m sk2bgrow.cli profile "res/A_${tag}"/*.counts.tsv \
      --db db --output "res/A_${tag}" --min-coverage 0 >/dev/null 2>&1

  /usr/bin/time -l $P profile "$fq" -d pileadb -o "res/C_${tag}" --single -t 8 \
      2>"res/C_${tag}.time" >/dev/null
  # Pilea's shipped --min-cove 5 returns nothing below 5x, so also run it with
  # the gates off; otherwise the low-coverage cells are empty by construction.
  /usr/bin/time -l $P profile "$fq" -d pileadb -o "res/D_${tag}" --single -t 8 \
      -x 0 -z 0 -c 0 2>"res/D_${tag}.time" >/dev/null

  rm -f "$fq" "res/A_${tag}"/*.counts.tsv "res/C_${tag}"/*.kmc "res/D_${tag}"/*.kmc
  touch "res/${tag}.done"
  echo "  done $tag"
done; done; done
echo "grid complete: $(ls res/*.done 2>/dev/null | wc -l) cells"
