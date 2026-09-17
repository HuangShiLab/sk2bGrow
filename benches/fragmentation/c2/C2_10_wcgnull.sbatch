#!/usr/bin/env bash
#SBATCH --job-name=C2_wcgnull
#SBATCH --partition=intel
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/logs/C2/wcgnull.%A.out
# C2 Task 1: score the WCG statistic on every C1 + C1b cell (all nulls)
# offline from existing windows.rates.tsv / output.tsv. No estimator re-runs.
set -euo pipefail
BASE=/lustre1/g/aos_shihuang/sk2bgrow-hpc
PY=$BASE/micromamba/envs/sk2bgrow/bin/python
cd $BASE/bench/C2/scripts
"$PY" c2_wcg_null.py
echo "DONE wcgnull"
