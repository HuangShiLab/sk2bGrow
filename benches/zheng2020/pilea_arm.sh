#!/usr/bin/env bash
# Pilea arm. Run twice per sample: at its shipped defaults (min-cove 5), and
# with the quality gates relaxed so it still reports below 5x -- otherwise the
# low-coverage comparison is "estimate vs no estimate" rather than a comparison.
cd "$(dirname "$0")"
P=../pilea_env/bin/pilea
mkdir -p out
for sfq in sub/*.fq; do
  s=$(basename "$sfq" .fq)
  for mode in default relaxed; do
    o="out/C_${mode}_${s}"
    [ -f "$o/profile.tsv" ] && continue
    if [ "$mode" = default ]; then extra=""; else extra="-x 0 -z 0 -c 0"; fi
    $P profile "$sfq" -d pileadb -o "$o" --single -t 8 $extra >"$o.log" 2>&1 || true
  done
done
echo "pilea runs: $(ls -d out/C_* 2>/dev/null | grep -v '\.log' | wc -l)"
ls out/C_default_*/ 2>/dev/null | head -5
