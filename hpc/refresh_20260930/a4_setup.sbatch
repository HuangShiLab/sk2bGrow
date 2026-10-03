#!/usr/bin/env bash
#SBATCH --job-name=A4setup
#SBATCH --partition=intel
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930/logs/a4_setup.%A.out
set -euo pipefail
source /lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930/env.sh
A4="$WORK/a4"; mkdir -p "$A4/work/out" "$A4/work/out_gc" "$A4/work/fq" "$A4/work/eff"
cp "$PAPER/benches/windowrate/analyze.py" "$A4/"
cp "$PAPER/benches/windowrate/simulate_tent.py" "$A4/"
ln -snf "$ECOLI_DB" "$A4/work/db"
python=python3
"$python" - "$A4" <<'PY'
import sys
from pathlib import Path
root=Path(sys.argv[1]); arms=['pois','nb']; depths=[0.5,1,2,5]; bs=[0,0.5,1,1.5,2]; reps=[0,1,2]
with (root/'work/cells.tsv').open('w') as fh:
    for arm_i,arm in enumerate(arms):
      for d_i,d in enumerate(depths):
       for b_i,b in enumerate(bs):
        for rep in reps:
         seed=10000000+arm_i*1000000+d_i*100000+b_i*10000+rep
         sample=f'{arm}_d{d:g}_b{b:g}_r{rep}'
         fh.write(f'{sample}\t{arm}\t{d}\t{b}\t{rep}\t{seed}\tout\t0\n')
    # GC controls, same underlying primary count files.
    for d in [0.5,1]:
      for b in bs:
        sample=f'pois_d{d:g}_b{b:g}_r0'
        fh.write(f'{sample}\tpois\t{d}\t{b}\t0\t0\tout_gc\t1\n')
PY
echo DONE
