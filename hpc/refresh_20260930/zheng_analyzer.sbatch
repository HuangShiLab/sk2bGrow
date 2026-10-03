#!/usr/bin/env bash
#SBATCH --job-name=Zanalysis
#SBATCH --partition=intel
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930/logs/zheng_analysis.%A.out
set -euo pipefail
source /lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930/env.sh
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
for seed_dir in "$WORK"/zheng/s*; do
  seed=$(basename "$seed_dir")
  mkdir -p "$seed_dir"
  cp "$PAPER/benches/zheng2020/analyze.py" "$seed_dir/"
  cp "$PAPER/benches/zheng2020/growth_rates.tsv" "$seed_dir/"
  cp "$WORK/selected_runs.tsv" "$seed_dir/"
  cd "$seed_dir"
  "$PY" analyze.py > "analyze_$seed.log" 2>&1
done
cd "$WORK/zheng"
"$PY" "$HPCDIR/summarize_refresh.py" \
  --zheng-root "$WORK/zheng" --a4 "$WORK/a4/results" --r3 "$WORK/r3" \
  --out "$WORK/results" > "$WORK/results/refresh_summary.log" 2>&1
echo DONE
