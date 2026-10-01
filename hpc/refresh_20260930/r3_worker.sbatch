#!/usr/bin/env bash
#SBATCH --job-name=R3refresh
#SBATCH --partition=intel
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930/logs/r3.%A_%a.out
set -euo pipefail
source /lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930/env.sh
R3="$WORK/r3"
ref=$(awk -v i=$((SLURM_ARRAY_TASK_ID+1)) 'NR==i{print}' "$R3/refs.txt")
name=$(basename "$ref" .fna)
db="$R3/db_$name"; outroot="$R3/out"
[[ "$name" == *R ]] && outroot="$R3/outX"
echo "ref=$ref db=$db"
if [ ! -s "$db/manifest.json" ]; then
  "$BIN" index "$ref" -o "$db" --enzymes all --threads 2 --quiet
fi
for fq in "$R3"/reads/*.fq "$R3"/readsR/*.fq "$R3"/readsRX/*.fq; do
  sample=$(basename "$fq" .fq)
  out="$outroot/${name}_${sample}"; mkdir -p "$out"
  if [ ! -s "$out/output.tsv" ]; then
    "$BIN" profile "$fq" -d "$db" -o "$out" --threads 2 --quiet --python "$PY" \
      --method v_shape --min-coverage 0 > "$out/run.log" 2> "$out/time.txt"
  fi
  rm -f "$out"/*.counts.tsv "$out"/*.em.tsv
  [ -s "$out/windows.rates.tsv" ] || { echo "missing rates $out"; exit 6; }
done
echo DONE
