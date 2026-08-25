#!/usr/bin/env bash
# A2 / experiment C2: does a fragmented reference break the coordinate V-fit,
# and does `sk2bgrow scaffold` put it back?
#
# Four reference conditions over the same reads (Zheng E. coli, 16 media + the
# stationary control, subsampled to 0.5-10x):
#
#   complete   the finished chromosome                       (upper bound)
#   frag       100 lognormal contigs, shuffled and flipped    (the MAG case)
#   scafSelf   frag scaffolded against the same genome        (upper bound on scaffolding)
#   scafRel    frag scaffolded against a different strain     (the realistic MAG case)
#
# plus Pilea gates-off on `frag`, since sorted-rank regression needs no
# coordinates and should be indifferent -- that is the claim to check, not assume.
set -uo pipefail
cd "$(dirname "$0")"
ROOT=${ROOT:-$HOME/Downloads/sk2bGrow}
BIN=$ROOT/target/release/sk2bgrow
WORK=${WORK:?set WORK to a scratch directory}
SUB=${SUB:?set SUB to the directory of subsampled .fq files}
GENOME=${GENOME:?set GENOME to the complete reference FASTA}
RELDB=${RELDB:?set RELDB to a database holding the relative}
SELF_NAME=${SELF_NAME:-Escherichia_coli_K12}
REL_NAME=${REL_NAME:-Escherichia_coli_O157H7}
export PYTHONPATH=$ROOT/python

mkdir -p "$WORK/out"

# --- references -------------------------------------------------------------
[ -f "$WORK/frag.fna" ] || python3 fragment.py "$GENOME" -n 100 --seed 0 -o "$WORK/frag.fna"
for pair in "scafSelf:$SELF_NAME" "scafRel:$REL_NAME"; do
  cond=${pair%%:*}; ref=${pair#*:}
  [ -f "$WORK/$cond.fna" ] && continue
  $BIN scaffold "$WORK/frag.fna" -d "$RELDB" -r "$ref" -o "$WORK/$cond.tgt" --quiet
  python3 rescaffold.py "$WORK/frag.fna" "$WORK/$cond.scaffold.json" --score \
      -o "$WORK/$cond.fna" --label ecoli | sed "s/^/  [$cond] /"
done
for cond in frag scafSelf scafRel; do
  [ -f "$WORK/db_$cond/manifest.json" ] || $BIN index "$WORK/$cond.fna" -o "$WORK/db_$cond" --quiet
done
[ -f "$WORK/db_complete/manifest.json" ] || $BIN index "$GENOME" -o "$WORK/db_complete" --quiet

# --- sk2bGrow over every condition -------------------------------------------
for f in "$SUB"/*.fq; do
  s=$(basename "$f" .fq)
  for cond in complete frag scafSelf scafRel; do
    o="$WORK/out/${cond}_${s}"
    [ -f "$o/output.tsv" ] && continue
    $BIN profile "$f" -d "$WORK/db_$cond" -o "$o" --quiet --python python3 \
        >/dev/null 2>&1 || echo "FAIL $cond $s"
  done
done
echo "sk2bgrow: $(ls -d "$WORK"/out/*/ 2>/dev/null | wc -l) runs"

# --- Pilea on the fragmented reference ---------------------------------------
PILEA=${PILEA:-$WORK/../pilea_env/bin/pilea}
if [ -x "$PILEA" ]; then
  [ -f "$WORK/pileadb_frag/sketches.pdb" ] || \
      $PILEA index "$WORK/frag.fna" -o "$WORK/pileadb_frag" -t 8 >/dev/null 2>&1
  for f in "$SUB"/*.fq; do
    s=$(basename "$f" .fq); o="$WORK/out/pileaFrag_${s}"
    [ -f "$o/profile.tsv" ] && continue
    $PILEA profile "$f" -d "$WORK/pileadb_frag" -o "$o" --single -t 8 -x 0 -z 0 -c 0 \
        >/dev/null 2>&1 || true
  done
  echo "pilea: $(ls -d "$WORK"/out/pileaFrag_*/ 2>/dev/null | wc -l) runs"
fi
echo A2 DONE
