#!/usr/bin/env bash
# Resolve each accession to its RefSeq FTP path and pull the genomic FASTA.
# Total is ~60 MB — deliberately light enough to iterate on a laptop.
cd "$(dirname "$0")"
while IFS=$'\t' read -r acc name; do
  out="genomes/${name}.fna"
  [ -s "$out" ] && { echo "have $name"; continue; }
  asm=$(curl -sSL --max-time 45 "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/${acc}/dataset_report" 2>/dev/null \
        | python3 -c "import json,sys; print(json.load(sys.stdin)['reports'][0]['assembly_info']['assembly_name'])" 2>/dev/null)
  [ -z "$asm" ] && { echo "FAIL resolve $acc"; continue; }
  p=$(echo "$acc" | sed 's/GCF_//; s/\..*//' | sed 's/\(...\)\(...\)\(...\)/\1\/\2\/\3/')
  url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${p}/${acc}_${asm}/${acc}_${asm}_genomic.fna.gz"
  curl -sSL --max-time 180 "$url" 2>/dev/null | gunzip -c > "$out.part" 2>/dev/null
  n=$(grep -c '^>' "$out.part" 2>/dev/null || echo 0)
  if [ "$n" -ge 1 ]; then mv "$out.part" "$out"; echo "  $name  $(du -h $out|cut -f1)  ${n} seq"
  else rm -f "$out.part"; echo "FAIL download $name"; fi
done < accs.txt
echo "=== $(ls genomes/*.fna 2>/dev/null | wc -l) genomes, $(du -sh genomes|cut -f1) ==="
