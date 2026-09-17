#!/usr/bin/env bash
#SBATCH --job-name=C2_refs
#SBATCH --partition=intel
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/logs/C2/refs.%A.out
# C2 reference construction: frag100 (fragment.py protocol, seed 0, identical
# to the laptop draft), incomplete/contaminated MAG variants, the scaffolding
# distance ladder (O157:H7 -> Shigella -> Salmonella -> Vibrio), the frag16 /
# complete16 simulation references, and every anchor DB. Writes refs.tsv and
# runs.tsv consumed by C2_30_profile.sh.
set -euo pipefail
BASE=/lustre1/g/aos_shihuang/sk2bgrow-hpc
BIN=$BASE/src/target/release/sk2bgrow
PY=$BASE/micromamba/envs/sk2bgrow/bin/python
C2=$BASE/bench/C2
G16=$BASE/refs/genomes16
K12=$G16/Escherichia_coli_K12.fna
cd $C2/scripts

# --- 1. the 100-contig draft (must reproduce the laptop's 43,707 anchors) ---
[ -f $C2/refs/frag100.fna ] || \
  "$PY" fragment.py "$K12" -n 100 --seed 0 -o $C2/refs/frag100.fna

# --- 2. incomplete / contaminated MAG variants (Task 2.1) -------------------
"$PY" c2_make_refs.py --frag $C2/refs/frag100.fna --genomes-dir "$G16" \
  --out $C2/refs --seed 0

# --- 3. frag16 references for the simulation grid (Task 2.3) -----------------
"$PY" c2_make_refs.py --frag16 --genomes-dir "$G16" --out $C2/refs/frag16 \
  --seed 0 -n 100

# --- 4. scaffolding distance ladder (Task 2.2) -------------------------------
# one DB holding the four scaffolding references, nearest to most distant
LAD=$C2/db/ladder
if [ ! -f $LAD/manifest.json ]; then
  "$BIN" index $G16/Escherichia_coli_O157H7.fna $G16/Shigella_dysenteriae.fna \
    $G16/Salmonella_enterica_LT2.fna $G16/Vibrio_cholerae.fna \
    -o $LAD --threads 16 --quiet
fi
# genome names inside the DB, in index order
mapfile -t NAMES < <("$PY" -c "
import json
m = json.load(open('$LAD/manifest.json'))
for g in m['genomes']: print(g['name'])")
echo "ladder DB genomes: ${NAMES[*]}"
for short in O157H7 Shigella Salmonella Vibrio; do
  [ -f $C2/refs/scaf_${short}.fna ] && continue
  ref=$(printf '%s\n' "${NAMES[@]}" | grep -i "$short" | head -1)
  [ -n "$ref" ] || { echo "ERROR: no ladder genome matches $short"; exit 1; }
  "$BIN" scaffold $C2/refs/frag100.fna -d $LAD -r "$ref" \
    -o $C2/refs/scaf_${short}.tgt --threads 16 --quiet
  JSON=$C2/refs/scaf_${short}.scaffold.json
  [ -f "$JSON" ] || JSON=$C2/refs/scaf_${short}.tgt.scaffold.json
  [ -f "$JSON" ] || { echo "ERROR: scaffold json not found for $short"; ls $C2/refs/; exit 1; }
  "$PY" rescaffold.py $C2/refs/frag100.fna "$JSON" --score \
    -o $C2/refs/scaf_${short}.fna --label ecoli \
    | tee $C2/refs/scaf_${short}.score.txt
done

# --- 5. anchor DBs -----------------------------------------------------------
index1() { # name fasta enzymes
  [ -f $C2/db/$1/manifest.json ] && return 0
  "$BIN" index "$2" -o $C2/db/$1 --enzymes "$3" --threads 16 --quiet
}
index1 frag100        $C2/refs/frag100.fna all
for cond in comp90 comp75 comp50 comp75_c10_O157H7 comp50_c10_O157H7 \
            comp75_c10_Salmonella; do
  index1 $cond $C2/refs/$cond.fna all
done
for short in O157H7 Shigella Salmonella Vibrio; do
  index1 scaf_${short} $C2/refs/scaf_${short}.fna all
done
# panel-size subsets x fragmentation (Task 2.4). Nested, density-ranked by
# usable anchors on K-12 (laptop subset definition was not recoverable --
# deviation recorded in METHODS.md): k2 = the pair whose index size (17,055
# anchors) matched the laptop k=2 exactly.
K2=CjeI,CjePI
K4=$K2,HaeIV,Hin4I
K8=$K4,BcgI,BslFI,AlfI,Bsp24I
index1 frag_k2 $C2/refs/frag100.fna $K2
index1 frag_k4 $C2/refs/frag100.fna $K4
index1 frag_k8 $C2/refs/frag100.fna $K8
index1 complete_k2 "$K12" $K2
index1 complete_k4 "$K12" $K4
index1 complete_k8 "$K12" $K8
# simulation DBs (Task 2.3)
[ -f $C2/db/frag16/manifest.json ] || \
  "$BIN" index $C2/refs/frag16 -o $C2/db/frag16 --threads 16 --quiet
[ -f $C2/db/complete16/manifest.json ] || \
  "$BIN" index "$G16" -o $C2/db/complete16 --threads 16 --quiet

# --- 6. manifests for the profile/score stages -------------------------------
{
  echo -e "condition\tfasta\tdb"
  echo -e "frag100\t$C2/refs/frag100.fna\t$C2/db/frag100"
  for cond in comp90 comp75 comp50 comp75_c10_O157H7 comp50_c10_O157H7 \
              comp75_c10_Salmonella; do
    echo -e "$cond\t$C2/refs/$cond.fna\t$C2/db/$cond"
  done
  for short in O157H7 Shigella Salmonella Vibrio; do
    echo -e "scaf_${short}\t$C2/refs/scaf_${short}.fna\t$C2/db/scaf_${short}"
  done
  for k in 2 4 8; do
    echo -e "frag_k$k\t$C2/refs/frag100.fna\t$C2/db/frag_k$k"
    echo -e "complete_k$k\t$K12\t$C2/db/complete_k$k"
  done
} > $C2/refs.tsv

# 17 runs: the laptop picks (16 media + one RUN_OUT replicate), depths 0.5-10x
{
  echo -e "run\tmedium"
  awk '{print $2"\t"$1}' $BASE/bench/C1/picks.tsv
} > $C2/runs.tsv

# --- 7. anchor-count check against the laptop numbers -------------------------
"$PY" - <<'EOF'
import json
base = "/lustre1/g/aos_shihuang/sk2bgrow-hpc/bench/C2/db"
expect = {"frag100": 43707}   # laptop: 43,707 of 43,735 survive the cuts
for name in sorted(__import__("os").listdir(base)):
    try:
        m = json.load(open(f"{base}/{name}/manifest.json"))
    except Exception:
        continue
    n = m["n_anchors"]
    flag = ""
    if name in expect:
        flag = "OK" if n == expect[name] else f"MISMATCH (laptop {expect[name]})"
    print(f"{name:24s} {n:8d} anchors {flag}")
EOF
echo "DONE refs"
