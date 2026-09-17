#!/usr/bin/env bash
# Experiment C7 / open issue A4: does the ZTP/ZTNB window-rate layer bias
# window rates at near-zero counts?
#
# Real E. coli K-12 index -> simulated tent-gradient reads with known truth ->
# real counting + window-rate pipeline -> analysis against exact per-window
# truth (analyze.py).
#
# Arms:   pois = Poisson counts;  nb = overdispersed (sigma_eff = 0.5 per kb).
# Depths: 0.5 1 2 5  (5x is the high-depth control)
# Truth:  log2PTR b in 0 0.5 1 1.5 2, three seeds per cell.
# GC control: pois arm at 0.5x/1x, rep 0, re-run through the stats layer with
# GC correction ON (reads are identical; no GC bias is simulated).
set -uo pipefail
cd "$(dirname "$0")"
ROOT="$(cd ../.. && pwd)"
BIN="$ROOT/target/release/sk2bgrow"
GENOME="$ROOT/benches/genomes/Escherichia_coli_K12.fna"
DB=work/db
PY="$ROOT/.venv/bin/python"
[ -x "$PY" ] || PY=python3
export PYTHONPATH="$ROOT/python"

[ -x "$BIN" ] || { echo "build the release binary first: cargo build --release"; exit 1; }

if [ ! -f "$DB/manifest.json" ]; then
  echo "indexing $GENOME"
  "$BIN" index "$GENOME" -o "$DB" --enzymes all --write-tgt --quiet || exit 1
fi

DEPTHS="0.5 1 2 5"
BS="0 0.5 1 1.5 2"
REPS="0 1 2"

bidx=-1
for b in $BS; do
  bidx=$((bidx + 1))
  didx=-1
  for depth in $DEPTHS; do
    didx=$((didx + 1))
    for arm in pois nb; do
      sigma=0; armidx=0
      if [ "$arm" = nb ]; then sigma=0.5; armidx=1; fi
      for rep in $REPS; do
        s="${arm}_d${depth}_b${b}_r${rep}"
        seed=$((10000000 + armidx * 1000000 + didx * 100000 + bidx * 10000 + rep))
        fq="work/fq/${s}.fq"
        out="work/out/${s}"
        if [ ! -s "$fq" ]; then
          if [ "$arm" = nb ]; then
            "$PY" simulate_tent.py --genome "$GENOME" --depth "$depth" --log2ptr "$b" \
              --seed "$seed" --out "$fq" --sigma-eff "$sigma" --eff-out "work/eff/${s}.npy" || { echo "FAIL sim $s"; continue; }
          else
            "$PY" simulate_tent.py --genome "$GENOME" --depth "$depth" --log2ptr "$b" \
              --seed "$seed" --out "$fq" || { echo "FAIL sim $s"; continue; }
          fi
        fi
        if [ ! -f "$out/${s}.counts.tsv" ]; then
          mkdir -p "$out"
          "$BIN" profile "$fq" -d "$DB" -o "$out" --quiet --no-stats >/dev/null 2>"$out/count.log" \
            || { echo "FAIL count $s"; continue; }
        fi
        if [ ! -f "$out/windows.rates.tsv" ]; then
          "$PY" -m sk2bgrow.cli profile "$out/${s}.counts.tsv" --db "$DB" --output "$out" \
            --no-gc-correct >/dev/null 2>"$out/stats.log" || echo "FAIL stats $s"
        fi
      done
    done
  done
done

# GC-correction control: same counts, stats layer with GC correction on.
for depth in 0.5 1; do
  for b in $BS; do
    s="pois_d${depth}_b${b}_r0"
    src="work/out/${s}/${s}.counts.tsv"
    out="work/out_gc/${s}"
    [ -f "$src" ] || continue
    if [ ! -f "$out/windows.rates.tsv" ]; then
      mkdir -p "$out"
      cp "$src" "$out/"
      "$PY" -m sk2bgrow.cli profile "$out/${s}.counts.tsv" --db "$DB" --output "$out" \
        >/dev/null 2>"$out/stats.log" || echo "FAIL gc $s"
    fi
  done
done

echo "runs: $(ls -d work/out/*/ 2>/dev/null | wc -l) primary, $(ls -d work/out_gc/*/ 2>/dev/null | wc -l) gc-control"
"$PY" analyze.py
