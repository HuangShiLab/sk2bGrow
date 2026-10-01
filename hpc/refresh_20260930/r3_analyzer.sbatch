#!/usr/bin/env bash
#SBATCH --job-name=R3analysis
#SBATCH --partition=intel
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930/logs/r3_analysis.%A.out
set -euo pipefail
source /lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930/env.sh
cd "$WORK/r3"
WORK="$WORK/r3" "$PY" r3_score.py > r3_analysis.log 2>&1
echo DONE
