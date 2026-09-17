#!/usr/bin/env bash
#SBATCH --job-name=C2_submit
#SBATCH --partition=intel
#SBATCH --time=48:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/logs/C2/submit.%A.out
# C2 grid submitter: polls the per-user squeue TASK count (array tasks count
# individually against MaxSubmitJobsPerUser=50) and submits the C2_30_profile
# chunks in documented priority order as slots free, then the Task-2.3 sim
# grid. Retries on QOSMaxSubmitJobPerUserLimit. Pattern copied from
# C5 scripts/02_submit_grid.sh (itself from C3 02_submit_chunk2.sh).
# Already queued by hand, NOT resubmitted here: chunk 1-8 (job 3945712),
# chunk 80-84 (job 3945655), task 0 (pilot job 3945569).
set -u
cd /lustre1/g/aos_shihuang/sk2bgrow-hpc/bench/C2
mkdir -p .submit
ntasks() { squeue -u $USER -h -r | wc -l; }

CHUNKS=("9-16%8" "17-24%8" "25-32%8" "33-40%8" "41-48%8" "49-56%8" "56-63%8" "64-71%8" "72-79%8")
for i in $(seq 1 520); do
  # one submission attempt per round, strict priority order
  pending=""
  for ci in "${!CHUNKS[@]}"; do
    [ -f ".submit/chunk_$ci" ] || { pending="chunk $ci ${CHUNKS[$ci]}"; break; }
  done
  if [ -z "$pending" ] && [ ! -f .submit/sim ]; then
    pending="sim 0-23%12"
  fi
  if [ -n "$pending" ]; then
    kind=${pending%% *}; rest=${pending#* }
    ci=${rest%% *}; range=${rest#* }
    [ "$kind" = sim ] && range=$ci
    limit=40; [ "$kind" = sim ] && limit=35
    if [ "$(ntasks)" -le "$limit" ]; then
      script=scripts/C2_30_profile.sh
      [ "$kind" = sim ] && script=scripts/C2_40_simgrid.sh
      if sbatch --array="$range" "$script"; then
        touch ".submit/$kind$([ "$kind" = chunk ] && echo "_$ci")"
        echo "submitted $kind $range (round $i)"
      else
        echo "round $i: $kind $range QOS-blocked, retry in 600 s"
        sleep 600   # QOS race: someone else took the slots; keep polling
      fi
    fi
  fi
  # pilot task 0 safety net: if the pilot job vanished without producing
  # output, submit task 0 on its own
  if [ ! -f .submit/task0 ] && ! squeue -u $USER -h | grep -q 3945569; then
    first_cell=$(awk 'NR>1{print $1; exit}' runs.tsv)
    if [ ! -s "counts/frag100/${first_cell}.0.5x/output.tsv" ]; then
      sbatch --array=0 scripts/C2_30_profile.sh && touch .submit/task0
    else
      touch .submit/task0
    fi
  fi
  nleft=$(ls .submit 2>/dev/null | wc -l)
  if [ "$nleft" -ge 11 ]; then echo ALL_SUBMITTED; exit 0; fi
  sleep 300
done
echo GAVE_UP; exit 1
