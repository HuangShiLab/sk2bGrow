#!/usr/bin/env bash
#SBATCH --job-name=C2_prof
#SBATCH --partition=intel
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/logs/C2/prof.%A_%a.out
# C2 real-read grid (Tasks 2.1, 2.2, 2.4): 17 runs x 5 depths = 85 array
# tasks, each looping every reference condition in refs.tsv (17 rows):
# frag100, 6 incomplete/contaminated MAG variants, 4 scaffolded-ladder refs,
# frag/complete at k = 2/4/8 enzymes. Reads: C1 subsampled Zheng fq (reuse,
# no re-subsampling). Count stage = sk2bgrow --no-stats; stats = FIXED
# estimator (same invocation as C1/C1b). Skip-if-exists makes reruns
# incremental.
set -euo pipefail
BASE=/lustre1/g/aos_shihuang/sk2bgrow-hpc
BIN=$BASE/src/target/release/sk2bgrow
PY=$BASE/micromamba/envs/sk2bgrow/bin/python
C2=$BASE/bench/C2
FQ=$BASE/bench/C1/fq
DEPTHS=(0.5 1 2 5 10)

t=$SLURM_ARRAY_TASK_ID
run=$(awk -v i=$((t / 5 + 1)) 'NR>1 && NR-1==i {print $1}' $C2/runs.tsv)
depth=${DEPTHS[$((t % 5))]}
cell="${run}.${depth}x"
[ -n "$run" ] || { echo "no run for task $t"; exit 1; }
[ -s "$FQ/${cell}_1.fq.gz" ] || { echo "no fq for $cell"; exit 1; }

while IFS=$'\t' read -r cond fasta db; do
  [ "$cond" = "condition" ] && continue
  out=$C2/counts/$cond/$cell
  mkdir -p "$out"
  if [ ! -s "$out/output.tsv" ]; then
    if ! ls "$out"/*.counts.tsv >/dev/null 2>&1; then
      /usr/bin/time -v "$BIN" profile "$FQ/${cell}_1.fq.gz" "$FQ/${cell}_2.fq.gz" \
        -d "$db" -o "$out" --no-stats --threads 8 --quiet \
        2> "$out/time.txt"
    fi
    /usr/bin/time -v "$PY" -m sk2bgrow.cli profile "$out"/*.counts.tsv \
      --db "$db" --output "$out" --windows "$out/windows.tsv" \
      > "$out/stats.log" 2> "$out/stats.time.txt"
  fi
  echo "cell $cell $cond done"
done < $C2/refs.tsv
echo "DONE task=$t cell=$cell"
