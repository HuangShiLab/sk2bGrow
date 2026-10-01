#!/usr/bin/env bash
#SBATCH --job-name=A4analysis
#SBATCH --partition=intel
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930/logs/a4_analysis.%A.out
set -euo pipefail
source /lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930/env.sh
cd "$WORK/a4"
"$PY" analyze.py > analyze.log 2>&1
echo DONE
