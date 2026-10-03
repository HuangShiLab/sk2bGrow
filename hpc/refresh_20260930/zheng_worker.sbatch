#!/usr/bin/env bash
#SBATCH --job-name=Zrefresh
#SBATCH --partition=intel
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930/logs/zheng.%A_%a.out
set -euo pipefail
source /lustre1/g/aos_shihuang/sk2bgrow-hpc/refresh_20260930/env.sh
MANI="$WORK/selected_runs.tsv"
SEED_IDX=$(( SLURM_ARRAY_TASK_ID / 17 ))
RIDX=$(( SLURM_ARRAY_TASK_ID % 17 + 1 ))
read RUN MEDIUM < <(awk -v i="$RIDX" 'NR==i{print $2"\t"$1}' "$MANI")
SEED="s$SEED_IDX"
GLEN=4641652
DEPTHS=(0.5 1 2 5 10)
CODES=(50 100 200 500 1000)
echo "run=$RUN medium=$MEDIUM seed=$SEED"

for di in "${!DEPTHS[@]}"; do
  depth=${DEPTHS[$di]}; code=${CODES[$di]}
  cell="${MEDIUM}_${SEED}_${depth}x"
  if [ "$SEED_IDX" = 0 ]; then
    f1="$C1_FQ/${RUN}.${depth}x_1.fq.gz"
    f2="$C1_FQ/${RUN}.${depth}x_2.fq.gz"
  else
    dir="$WORK/fq/$SEED"; mkdir -p "$dir"
    f1="$dir/${RUN}.${depth}x_1.fq.gz"; f2="$dir/${RUN}.${depth}x_2.fq.gz"
    p1="$C1_FQ/${RUN}.20x_1.fq.gz"; p2="$C1_FQ/${RUN}.20x_2.fq.gz"
    if [ ! -s "$f1" ] || [ ! -s "$f2" ]; then
      # seqkit stats -T columns: file, format, type, num_seqs, sum_len, ...
      pool_bases=$("$SEQKIT" stats -T "$p1" "$p2" 2>/dev/null | awk 'NR>1{sum+=$5} END{print sum}')
      frac=$(awk -v d="$depth" -v g="$GLEN" -v b="$pool_bases" 'BEGIN{printf "%.12f",d*g/b}')
      rs=$((900000 + RIDX*100 + code + SEED_IDX*10))
      for mate in 1 2; do
        src=$([ "$mate" = 1 ] && echo "$p1" || echo "$p2")
        tmp="$dir/.${RUN}.${depth}x_${mate}.$SLURM_JOB_ID.fq.gz"
        "$SEQKIT" sample -j 2 -p "$frac" -s "$rs" "$src" -o "$tmp"
        mv "$tmp" "${f1%_1.fq.gz}_${mate}.fq.gz"
      done
    fi
    h1=$("$SEQKIT" seq -n -i "$f1" | md5sum | cut -d' ' -f1)
    h2=$("$SEQKIT" seq -n -i "$f2" | md5sum | cut -d' ' -f1)
    [ "$h1" = "$h2" ] || { echo "mate sync failed: $f1 $f2"; exit 2; }
  fi

  # Arm A: anchors + current auto/V-fit. Default QC remains active; failed rows stay in output.
  a="$WORK/zheng/$SEED/out/A_$cell"
  if [ ! -s "$a/output.tsv" ]; then
    mkdir -p "$a"
    /usr/bin/time -v "$BIN" profile "$f1" "$f2" -d "$ECOLI_DB" -o "$a" \
      --threads 4 --python "$PY" --quiet > "$a/run.log" 2> "$a/time.txt"
  fi
  rm -f "$a"/*.counts.tsv "$a"/*.em.tsv

  # Arm B: anchors + 25 kb windows + sorted/RANSAC parity.
  b="$WORK/zheng/$SEED/out/B_$cell"
  if [ ! -s "$b/output.tsv" ]; then
    mkdir -p "$b"
    "$BIN" profile "$f1" "$f2" -d "$ECOLI_DB" -o "$b" --threads 4 --no-stats \
      --windowing bp --window-bp 25000 --quiet
    "$PY" -m sk2bgrow.cli profile "$b"/*.counts.tsv --db "$ECOLI_DB" \
      --output "$b" --windows "$b/windows.tsv" --use-rust-windows --method sorted \
      --min-coverage 0 > "$b/stats.log" 2> "$b/stats.time.txt"
  fi
  rm -f "$b"/*.counts.tsv "$b"/*.em.tsv

  # Arm E: Pilea FracMinHash sketch rewritten into sk2bGrow count format.
  e="$WORK/zheng/$SEED/out/E_$cell"
  if [ ! -s "$e/output.tsv" ]; then
    ecdir="$WORK/zheng/$SEED/ecounts/$cell"
    mkdir -p "$e" "$ecdir"
    rm -f "$ecdir"/*.counts.tsv
    "$PPILEA" "$PAPER/benches/zheng2020/armE_counts.py" -r "$ECOLI" \
      --sample-name "$cell" -o "$ecdir" "$f1" "$f2" > "$e/armE.log" 2>&1
    raw=$(find "$ecdir" -maxdepth 1 -name '*.counts.tsv' | head -1)
    [ -n "$raw" ] || { echo "armE produced no counts: $cell"; exit 3; }
    fixed="$raw"
    "$PY" -m sk2bgrow.cli profile "$fixed" --db "$ECOLI_DB" --output "$e" \
      --min-coverage 0 > "$e/stats.log" 2> "$e/stats.time.txt"
    rm -f "$ecdir"/*.counts.tsv "$ecdir"/sketch_loci.pkl
  fi

  # Pilea arms.
  for mode in default relaxed; do
    p="$WORK/zheng/$SEED/out/C_${mode}_$cell"
    if [ ! -s "$p/output.tsv" ]; then
      mkdir -p "$p"
      if [ "$mode" = relaxed ]; then extra=(-x 0 -z 0 -c 0); else extra=(); fi
      /usr/bin/time -v "$PILEA" profile "$f1" "$f2" -d "$PILEA_DB" -o "$p" \
        -t 4 "${extra[@]}" > "$p/pilea.log" 2> "$p/time.txt"
    fi
  done
  echo "done $cell"
done
echo DONE
