#!/usr/bin/env bash
# For each still-short medium, try every replicate run until one yields the
# target read count. ENA throttles some paths unpredictably; a sibling run is
# usually fine.
cd "$(dirname "$0")"
N=600000; LINES=$((N*4))
for m in M6 M13 M19 M24 M25 M27; do
  best=0; bestf=""
  for f in fq/${m}.*.fq; do [ -s "$f" ] && { n=$(( $(wc -l < "$f")/4 )); [ $n -gt $best ] && { best=$n; bestf=$f; }; }; done
  [ "$best" -ge "$N" ] && { echo "$m ok ($best)"; continue; }
  for run in $(awk -F'\t' -v m=$m 'NR>1{d=$3; sub(/.*gDNA_/,"",d); if(d==m) print $1}' runs.tsv); do
    url=$(awk -F'\t' -v r=$run '$1==r{split($5,a,";"); print a[1]}' runs.tsv)
    [ -z "$url" ] && continue
    echo "$m: trying $run"
    curl -sS --retry 5 --retry-delay 5 --retry-all-errors --speed-time 60 --speed-limit 20000 \
         --max-time 3600 "https://${url}" 2>/dev/null | gunzip -c 2>/dev/null | head -n $LINES > "fq/${m}.${run}.fq.part"
    n=$(( $(wc -l < "fq/${m}.${run}.fq.part")/4 ))
    echo "   -> $n reads"
    if [ "$n" -ge "$N" ]; then
      rm -f fq/${m}.*.fq; mv "fq/${m}.${run}.fq.part" "fq/${m}.${run}.fq"; break
    elif [ "$n" -gt "$best" ]; then
      best=$n; rm -f fq/${m}.*.fq; mv "fq/${m}.${run}.fq.part" "fq/${m}.${run}.fq"
    else rm -f "fq/${m}.${run}.fq.part"; fi
  done
done
echo "=== final ==="; for f in fq/*.fq; do printf "%-24s %8d\n" "$(basename $f .fq)" $(( $(wc -l < $f)/4 )); done | sort -k2 -n
