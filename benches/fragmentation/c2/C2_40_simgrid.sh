#!/usr/bin/env bash
#SBATCH --job-name=C2_sim
#SBATCH --partition=intel
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/logs/C2/sim.%A_%a.out
# C2 Task 2.3: the laptop 16-genome multi-strain simulation grid (Pilea
# Methods design: ori at 0, ter at midpoint, V-profile, log2PTR ~ U[0,2],
# 150 bp single-end) with FRAGMENTED references. Reads are simulated on the
# COMPLETE genomes (truth on the real coordinate) and profiled against the
# frag16 DB (100-contig drafts, fragment.py protocol seed 0); the complete16
# arm on identical reads is the control.
# Array 0-23: s in {4,8,16} x c in {1,2,4,8} x r in {1,2}, laptop seeds.
set -euo pipefail
BASE=/lustre1/g/aos_shihuang/sk2bgrow-hpc
BIN=$BASE/src/target/release/sk2bgrow
PY=$BASE/micromamba/envs/sk2bgrow/bin/python
C2=$BASE/bench/C2
G16=$BASE/refs/genomes16
STRAINS=(4 8 16); COVS=(1 2 4 8)

t=$SLURM_ARRAY_TASK_ID
s=${STRAINS[$((t / 8))]}
r=$(( (t % 8) / 4 + 1 ))
c=${COVS[$((t % 4))]}
tag="s${s}_c${c}_r${r}"
pick_seed=$((s * 1000 + c * 10 + r))
seed=$((pick_seed * 100 + 1))
mkdir -p $C2/sim/res $C2/sim/tmp
fq=$C2/sim/tmp/${tag}.fq
[ -f $C2/sim/res/${tag}.done ] && { echo "skip $tag"; exit 0; }

"$PY" $C2/scripts/simulate_c3.py --genomes "$G16" --n-strains $s \
  --coverage $c --pick-seed $pick_seed --seed $seed \
  --out "$fq" --truth $C2/sim/res/${tag}.truth \
  --picks $C2/sim/res/picks_${tag}.tsv

for arm in frag16 complete16; do
  out=$C2/sim/res/${arm}_${tag}
  mkdir -p "$out"
  [ -s "$out/output.tsv" ] && continue
  /usr/bin/time -v "$BIN" profile "$fq" -d $C2/db/$arm -o "$out" \
    --no-stats --threads 8 --quiet 2> "$out/time.txt"
  /usr/bin/time -v "$PY" -m sk2bgrow.cli profile "$out"/*.counts.tsv \
    --db $C2/db/$arm --output "$out" --windows "$out/windows.tsv" \
    > "$out/stats.log" 2> "$out/stats.time.txt"
  rm -f "$out"/*.counts.tsv
done
rm -f "$fq"
touch $C2/sim/res/${tag}.done
echo "DONE $tag"
