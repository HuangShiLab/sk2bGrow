#!/usr/bin/env bash
set -euo pipefail
export REFRESH_BASE=/lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930
export SRC="$REFRESH_BASE/src"
export HPCDIR="$REFRESH_BASE/hpc/refresh_20260930"
export WORK="$REFRESH_BASE/work"
export LOGS="$REFRESH_BASE/logs"
export BIN="$SRC/target/release/sk2bgrow"
export PY=/lustre1/g/aos_shihuang/sk2bgrow-hpc/micromamba/envs/sk2bgrow/bin/python
export PPILEA=/lustre1/g/aos_shihuang/sk2bgrow-hpc/bench/C1/vendor/env-pilea138/bin/python
export PILEA=/lustre1/g/aos_shihuang/sk2bgrow-hpc/bench/C1/vendor/env-pilea138/bin/pilea
export SEQKIT=/lustre1/g/aos_shihuang/sk2bgrow-hpc/tools/seqkit
export ECOLI=/lustre1/g/aos_shihuang/sk2bgrow-hpc/refs/Escherichia_coli_K12.fna
export ECOLI_DB=/lustre1/g/aos_shihuang/sk2bgrow-hpc/bench/C1/db
export PILEA_DB=/lustre1/g/aos_shihuang/sk2bgrow-hpc/bench/C1/pileadb
export C1_FQ=/lustre1/g/aos_shihuang/sk2bgrow-hpc/bench/C1/fq
export PYTHONPATH="$SRC/python"
mkdir -p "$WORK" "$LOGS"
