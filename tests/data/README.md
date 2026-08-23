# Test data

Deliberately empty of real genomes.

Every test in this repository builds its own fixture, for two reasons: a
reference genome is tens of megabytes and does not belong in git, and a
generated fixture states its own ground truth — a planted PTR, a planted origin,
a planted GC gradient — so a test can assert recovery rather than
self-consistency.

Where the fixtures come from:

| suite | fixture |
|---|---|
| `crates/sk2bgrow-core/src/**` unit tests | hand-written sequences with planted recognition sites |
| `crates/sk2bgrow-core/tests/pipeline.rs` | xorshift-generated genomes plus simulated reads |
| `crates/sk2bgrow-cli/tests/cli.rs` | same, driven through the real binary |
| `tests/python/conftest.py` | `make_counts()` — synthetic count tables with a planted gradient |
| `scripts/smoke.sh` | a 1.5 Mb genome and reads sampled in proportion to copy number |

The pseudo-random generators are seeded, so every fixture is reproducible.

## Real data

For benchmarking against real references and reads, see
[`../../benches/README.md`](../../benches/README.md). Those datasets are
downloaded on demand into `benches/work/`, which is git-ignored.
