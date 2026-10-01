#!/usr/bin/env bash
#SBATCH --job-name=Zstage
#SBATCH --partition=intel
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=24:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930/logs/staged_submit.%A.out
bash /lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930/staged_submit_refresh.sh
