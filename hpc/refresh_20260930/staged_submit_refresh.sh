#!/usr/bin/env bash
set -euo pipefail
source /lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930/env.sh
ntasks() { squeue -u "$USER" -h -r | wc -l; }
ids_z=(); ids_a=()
chunks=(
  "zheng 1-10" "zheng 11-20" "zheng 21-30" "zheng 31-40" "zheng 41-50"
  "a4 1-10" "a4 11-20" "a4 21-30" "a4 31-40" "a4 41-50"
  "a4 51-60" "a4 61-70" "a4 71-80" "a4 81-90" "a4 91-100"
  "a4 101-110" "a4 111-120" "a4 121-129"
)
for chunk in "${chunks[@]}"; do
  kind=${chunk%% *}; range=${chunk#* }
  script="$HPCDIR/zheng_worker.sbatch"; [[ "$kind" = a4 ]] && script="$HPCDIR/a4_worker.sbatch"
  while [ "$(ntasks)" -gt 20 ]; do sleep 180; done
  if [ "$kind" = zheng ]; then
    id=$(sbatch --array="$range%10" "$script" | awk '{print $4}')
    ids_z+=("$id")
  else
    id=$(sbatch --array="$range%10" "$script" | awk '{print $4}')
    ids_a+=("$id")
  fi
  echo "$(date -Is) submitted $kind $range id=$id queue=$(ntasks)"
  sleep 30
done
dep_z=$(IFS=:; echo "${ids_z[*]}")
dep_a=$(IFS=:; echo "${ids_a[*]}")
sbatch --dependency=afterany:4193591:"$dep_z" "$HPCDIR/zheng_analyzer.sbatch"
sbatch --dependency=afterany:4193542:"$dep_a" "$HPCDIR/a4_analyzer.sbatch"
printf 'ZHENG_CHUNKS=%s\nA4_CHUNKS=%s\n' "$dep_z" "$dep_a" > "$REFRESH_BASE/staged_ids.env"
echo ALL_SUBMITTED
