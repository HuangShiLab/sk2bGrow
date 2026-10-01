#!/usr/bin/env bash
set -euo pipefail
source /lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930/env.sh
ids=(); for i in $(seq 1 16); do ids+=("$i"); done; for i in $(seq 18 50); do ids+=("$i"); done
for id in "${ids[@]}"; do
  while [ "$(squeue -u "$USER" -h -r | wc -l)" -gt 20 ]; do sleep 120; done
  j=$(sbatch --array="$id" "$HPCDIR/zheng_worker.sbatch" | awk '{print $4}')
  echo "$(date -Is) task=$id job=$j queue=$(squeue -u "$USER" -h -r|wc -l)"
  sleep 20
done
echo ALL_ZHENG_INDIVIDUAL
