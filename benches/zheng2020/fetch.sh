#!/usr/bin/env bash
# Download the first N R1 reads per sample. Subsampling to <=10x anyway, so
# there is no point pulling 490x. SRA preserves original (flowcell) read order,
# which is random with respect to genome position -- fine for coverage profiles.
cd "$(dirname "$0")"
N_READS=600000            # ~90 Mbp ~= 19x of the 4.64 Mb genome
LINES=$((N_READS*4))
mkdir -p fq
# one run per medium, first replicate
awk -F'\t' 'NR>1{d=$3; sub(/.*gDNA_/,"",d); if(!(d in seen)){seen[d]=1; split($5,f,";"); print d"\t"$1"\t"f[1]}}' runs.tsv > picks.tsv
while IFS=$'\t' read -r medium run url; do
  out="fq/${medium}.${run}.fq"
  [ -s "$out" ] && { echo "have $medium"; continue; }
  curl -sSL --max-time 900 "https://${url}" 2>/dev/null | gunzip -c 2>/dev/null | head -n $LINES > "$out.part" || true
  n=$(( $(wc -l < "$out.part") / 4 ))
  mv "$out.part" "$out"
  echo "$medium $run  $n reads  $(du -h "$out" | cut -f1)"
done < picks.tsv
echo "DONE: $(ls fq/*.fq 2>/dev/null | wc -l) samples, $(du -sh fq | cut -f1)"
