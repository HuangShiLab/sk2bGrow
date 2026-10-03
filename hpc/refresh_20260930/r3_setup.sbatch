#!/usr/bin/env bash
#SBATCH --job-name=R3setup
#SBATCH --partition=intel
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930/logs/r3_setup.%A.out
set -euo pipefail
source /lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930/env.sh
R3="$WORK/r3"; mkdir -p "$R3/refs" "$R3/refsR" "$R3/refsRX" "$R3/reads" "$R3/readsR" "$R3/readsRX"
cp "$PAPER/benches/fragmentation/"*.py "$R3/"
cp "$PAPER/benches/simulate.py" "$R3/"
cp "$ECOLI" "$R3/refs/complete.fna"
cd "$R3"
"$PY" r3_refs.py "$ECOLI" -o refs --counts "2 5 10 20 50 100" --ordered "10 100" --seed 0
"$PY" r3_refs.py "$ECOLI" -o refsR --counts "2 5 10 20 50 100" --ordered "10 100" --seed 0 --rotate 0.05 --suffix R
"$PY" r3_refs.py "$ECOLI" -o refsRX --counts "2 5 10 20 50 100" --ordered "10 100" --seed 0 --rotate 0.15 --suffix RX
"$PY" r3_simulate.py "$ECOLI" -o reads --seed 1000
"$PY" r3_simulate.py refsR/completeR.fna -o readsR --seed 2000
"$PY" r3_simulate.py refsRX/completeRX.fna -o readsRX --seed 3000
: > refs.txt
for d in refs refsR refsRX; do find "$R3/$d" -maxdepth 1 -name '*.fna' -printf "$R3/$d/%f\n"; done | sort >> "$R3/refs.txt"
echo "refs=$(wc -l < refs.txt)"
echo DONE
