#!/usr/bin/env bash
#SBATCH --job-name=A4refresh
#SBATCH --partition=intel
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930/logs/a4.%A_%a.out
set -euo pipefail
source /lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930/env.sh
A4="$WORK/a4"
line=$(awk -v i=$((SLURM_ARRAY_TASK_ID+1)) 'NR==i{print}' "$A4/work/cells.tsv")
read SAMPLE ARM DEPTH B REP SEED OUTROOT GCFLAG < <(printf '%s\n' "$line")
dtag=$(printf '%g' "$DEPTH"); btag=$(printf '%g' "$B")
echo "$line"
if [ "$GCFLAG" = 0 ]; then
  out="$A4/$OUTROOT/$SAMPLE"; mkdir -p "$out"
  fq="$A4/work/fq/$SAMPLE.fq"; eff="$A4/work/eff/$SAMPLE.npy"
  sigma=0; sigma_args=()
  if [ "$ARM" = nb ]; then sigma=0.5; sigma_args=(--sigma-eff 0.5 --eff-out "$eff"); fi
  if [ ! -s "$out/$SAMPLE.counts.tsv" ]; then
    [ -s "$fq" ] || "$PY" "$PAPER/benches/windowrate/simulate_tent.py" --genome "$ECOLI" \
      --depth "$DEPTH" --log2ptr "$B" --seed "$SEED" "${sigma_args[@]}" --out "$fq"
    "$BIN" profile "$fq" -d "$ECOLI_DB" -o "$out" --threads 2 --no-stats --quiet \
      > "$out/count.log" 2> "$out/count.time.txt"
  fi
  if [ ! -s "$out/windows.rates.tsv" ]; then
    "$PY" -m sk2bgrow.cli profile "$out/$SAMPLE.counts.tsv" --db "$ECOLI_DB" \
      --output "$out" --no-gc-correct > "$out/stats.log" 2> "$out/stats.time.txt"
  fi
else
  src_sample="$SAMPLE"; out="$A4/$OUTROOT/$SAMPLE"; mkdir -p "$out"
  src="$A4/out/$src_sample/$src_sample.counts.tsv"
  [ -s "$src" ] || { echo "missing GC source $src"; exit 4; }
  cp "$src" "$out/$SAMPLE.counts.tsv"
  [ -s "$out/windows.rates.tsv" ] || "$PY" -m sk2bgrow.cli profile "$out/$SAMPLE.counts.tsv" \
    --db "$ECOLI_DB" --output "$out" > "$out/stats.log" 2> "$out/stats.time.txt"
fi
[ -s "$out/windows.rates.tsv" ] || { echo "missing rates: $SAMPLE"; exit 5; }
echo DONE
