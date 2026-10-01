#!/usr/bin/env bash
# Current-commit mixed-strain validation: 4/8/16 strains x 1/2/4x x 2 reps.
# This is deliberately smaller than the historical C3 grid so it can run
# repeatably on a workstation while exercising the full production stack.
set -euo pipefail
cd "$(dirname "$0")"
ROOT="$(cd ../.. && pwd)"
PY="${PY:-$ROOT/.venv/bin/python}"
BIN="${BIN:-$ROOT/target/release/sk2bgrow}"
GENOMES="${GENOMES:-$ROOT/benches/genomes}"
STRAINS="${STRAINS:-4 8 16}"
COVS="${COVS:-1 2 4}"
REPS="${REPS:-2}"
mkdir -p db sim res logs

if [[ ! -f db/manifest.json ]]; then
  "$BIN" index "$GENOMES"/*.fna -o db --enzymes all --quiet
fi

for s in $STRAINS; do
  for c in $COVS; do
    for r in $(seq 1 "$REPS"); do
      tag="s${s}_c${c}_r${r}"
      [[ -f "res/${tag}.done" ]] && continue
      fq="sim/${tag}.fq"
      "$PY" "$ROOT/benches/simulate.py" \
        --genomes "$GENOMES" --n-strains "$s" --coverage "$c" \
        --seed "$((s*1000 + c*10 + r))" \
        --out "$fq" --truth "res/${tag}.truth" > "logs/${tag}.simulate.log"
      out="res/A_${tag}"
      mkdir -p "$out"
      "$BIN" profile "$fq" -d db -o "$out" --threads 4 --no-stats --quiet \
        > "logs/${tag}.count.log" 2> "logs/${tag}.count.time"
      PYTHONPATH="$ROOT/python" "$PY" -m sk2bgrow.cli profile \
        "$out"/*.counts.tsv --db db --output "$out" --min-coverage 0 \
        --method v_shape \
        > "logs/${tag}.stats.log" 2> "logs/${tag}.stats.time"
      rm -f "$out"/*.counts.tsv "$out"/*.em.tsv "$fq"
      touch "res/${tag}.done"
      echo "done $tag"
    done
  done
done

# score.py includes absent Pilea arms from the historical grid; this validation
# is current-sk2bGrow only.
PYTHONPATH="$ROOT/python" "$PY" "$PWD/summarize_mixed.py"
