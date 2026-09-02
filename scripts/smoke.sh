#!/usr/bin/env bash
# Full-stack smoke test: synthetic genome -> index -> reads -> profile -> output.tsv.
#
# This is the one check that exercises the Rust/Python handoff for real. The test
# suites cover each half; only this covers the seam.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Prefer the release binary when it exists: counting is the slow step and the
# debug build is ~20x slower.
if [[ -n "${SK2BGROW_BIN:-}" ]]; then
    BIN="$SK2BGROW_BIN"
elif [[ -x "$ROOT/target/release/sk2bgrow" ]]; then
    BIN="$ROOT/target/release/sk2bgrow"
else
    BIN="$ROOT/target/debug/sk2bgrow"
fi
PY="${PY:-python3}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/sk2bgrow-smoke.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

if [[ ! -x "$BIN" ]]; then
    echo "error: $BIN not found; run 'make build' first" >&2
    exit 1
fi

echo "==> work dir: $WORK"

# --- a synthetic 1.5 Mb genome with a planted replication gradient -----------
"$PY" - "$WORK" <<'PYEOF'
import sys, random
from pathlib import Path

work = Path(sys.argv[1])
rng = random.Random(20240823)
LEN = 1_500_000                   # big enough that every enzyme clears the
ORI = 150_000                     # per-enzyme window minimum, as on a real genome
LOG2_PTR = 1.0

genome = "".join(rng.choice("ACGT") for _ in range(LEN))
(work / "genome.fna").write_text(
    ">synth_chr\n" + "\n".join(genome[i:i + 70] for i in range(0, LEN, 70)) + "\n"
)

def copy_number(pos):
    d = min((pos - ORI) % LEN, (ORI - pos) % LEN)
    return 2.0 ** (-LOG2_PTR * d / (LEN / 2))

# Tiling 150 bp reads, sampled in proportion to local copy number.
reads = []
step, base_depth = 40, 9
for start in range(0, LEN - 150, step):
    n = base_depth * copy_number(start + 75)
    k = int(n) + (1 if rng.random() < n - int(n) else 0)
    for _ in range(k):
        r = genome[start:start + 150]
        if rng.random() < 0.5:
            r = r.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        reads.append(r)
rng.shuffle(reads)

with open(work / "S1.fq", "w") as fh:
    for i, r in enumerate(reads):
        fh.write(f"@r{i}\n{r}\n+\n{'I' * len(r)}\n")
print(f"    genome {LEN:,} bp, ori {ORI:,}, planted log2(PTR) {LOG2_PTR}, {len(reads):,} reads")
PYEOF

echo "==> index"
"$BIN" index "$WORK/genome.fna" -o "$WORK/db" --enzymes all --write-tgt

echo "==> audit"
"$BIN" audit "$WORK/db" -o "$WORK/audit.tsv"

echo "==> profile (Rust counting + Python statistics)"
PYTHONPATH="$ROOT/python" "$BIN" profile "$WORK/S1.fq" -d "$WORK/db" -o "$WORK/out" --python "$PY"

echo "==> verify"
PYTHONPATH="$ROOT/python" "$PY" - "$WORK" <<'PYEOF'
import sys
from pathlib import Path
import pandas as pd

work = Path(sys.argv[1])
df = pd.read_csv(work / "out" / "output.tsv", sep="\t", na_values=["NA"])
assert len(df) == 1, f"expected one genome-sample row, got {len(df)}"
row = df.iloc[0]
print(row.to_string())

problems = []
if not (0.55 <= row["log2(PTR)"] <= 1.45):
    problems.append(f"log2(PTR) {row['log2(PTR)']:.3f} is far from the planted 1.0")
if row["n_enzymes"] < 14:
    problems.append(f"only {row['n_enzymes']} enzymes contributed")
if abs(((row["ori"] - 150_000) % 1_500_000 + 750_000) % 1_500_000 - 750_000) > 150_000:
    problems.append(f"ori {row['ori']:.0f} is far from the planted 150,000")
if problems:
    print("\nFAIL:")
    for p in problems:
        print("  -", p)
    sys.exit(1)
print("\nOK: recovered the planted gradient end to end.")
PYEOF
