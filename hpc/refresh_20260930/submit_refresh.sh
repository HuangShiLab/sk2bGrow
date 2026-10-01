#!/usr/bin/env bash
set -euo pipefail
source /lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930/env.sh
cd "$SRC"
mkdir -p "$WORK/results" "$WORK/zheng" "$WORK/fq" "$WORK/a4" "$WORK/r3"
# 17 selected media plus RUN_OUT; source picks.tsv columns are medium/run/url.
awk 'BEGIN{FS=OFS="\t"}{print $1,$2}' "$SRC/benches/zheng2020/picks.tsv" | sort -k2,2 > "$WORK/selected_runs.tsv"
[ "$(wc -l < "$WORK/selected_runs.tsv")" = 17 ] || { echo 'expected 17 selected runs'; exit 1; }

jz=$(sbatch --array=0-50%25 "$HPCDIR/zheng_worker.sbatch" | awk '{print $4}')
sbatch --dependency=afterany:"$jz" "$HPCDIR/zheng_analyzer.sbatch"

sbatch "$HPCDIR/a4_setup.sbatch"
ja=$(sbatch --array=0-129%30 "$HPCDIR/a4_worker.sbatch" | awk '{print $4}')
sbatch --dependency=afterany:"$ja" "$HPCDIR/a4_analyzer.sbatch"

jr=$(sbatch "$HPCDIR/r3_setup.sbatch" | awk '{print $4}')
# setup lines count is known: 9 refs per family = 27.
jrw=$(sbatch --dependency=afterok:"$jr" --array=0-26%20 "$HPCDIR/r3_worker.sbatch" | awk '{print $4}')
sbatch --dependency=afterany:"$jrw" "$HPCDIR/r3_analyzer.sbatch"

cat > "$REFRESH_BASE/job_ids.env" <<EID
ZHENG_WORKER=$jz
A4_WORKER=$ja
R3_SETUP=$jr
R3_WORKER=$jrw
EID
echo "submitted: Zheng=$jz A4=$ja R3setup=$jr R3worker=$jrw"
