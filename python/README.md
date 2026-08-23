# sk2bgrow (Python layer)

The statistics half of [sk2bGrow](../README.md): everything from a per-anchor
count table to a PTR estimate.

```bash
pip install -e .
python -m sk2bgrow.cli profile out/*.counts.tsv --db db --output out/
```

Normally invoked by the Rust binary (`sk2bgrow profile` shells out to it), but it
runs standalone against any count table — which is the point of putting the
layer boundary at a file. Iterating on the statistical model costs a re-run of
this layer, not a re-count of the reads.

| module | role |
|---|---|
| `io` | read the interface files written by the Rust layer |
| `ztp` | zero-truncated Poisson / negative-binomial window rates |
| `gc_bias` | per-enzyme loess GC correction, with shrinkage |
| `fit` | V-shape MLE on real coordinates; sorted+RANSAC for parity |
| `fusion` | inverse-variance fusion across enzymes + Cochran's Q |
| `dynamics` | ΔPTR, anchor × sample matrix, trend tests |
| `report` | `output.tsv` and QC figures |
| `simulate` | Monte-Carlo harness reproducing the design report's §5 |

See [`../docs/design/02-algorithm.md`](../docs/design/02-algorithm.md) for the
statistical reasoning behind each step.
