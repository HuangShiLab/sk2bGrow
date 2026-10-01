#!/usr/bin/env bash
set -euo pipefail
source /lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930/env.sh
ntasks() { squeue -u "$USER" -h -r | wc -l; }
chunks=(
  "zheng 1-16" "zheng 18-30" "zheng 31-45" "zheng 46-50"
  "a4 1-10" "a4 11-20" "a4 21-30" "a4 31-40" "a4 41-50"
  "a4 51-60" "a4 61-70" "a4 71-80" "a4 81-90" "a4 91-100"
  "a4 101-110" "a4 111-120" "a4 121-129"
)
zids=(); aids=()
for chunk in "${chunks[@]}"; do
  kind=${chunk%% *}; range=${chunk#* }
  while [ "$(ntasks)" -gt 20 ]; do sleep 180; done
  script="$HPCDIR/zheng_worker.sbatch"; [[ "$kind" = a4 ]] && script="$HPCDIR/a4_worker.sbatch"
  id=$(sbatch --array="$range%10" "$script" | awk '{print $4}')
  if [ "$kind" = zheng ]; then zids+=("$id"); else aids+=("$id"); fi
  echo "$(date -Is) $kind $range id=$id queue=$(ntasks)"
  sleep 60
done
dz=$(IFS=:; echo "${zids[*]}"); da=$(IFS=:; echo "${aids[*]}")
sbatch --dependency=afterany:"$dz" "$HPCDIR/zheng_analyzer.sbatch"
sbatch --dependency=afterany:"$da" "$HPCDIR/a4_analyzer.sbatch"
printf 'RERUN_ZHENG=%s\nRERUN_A4=%s\n' "$dz" "$da" > "$REFRESH_BASE/rerun_staged_ids.env"
echo ALL_RERUN_SUBMITTED
