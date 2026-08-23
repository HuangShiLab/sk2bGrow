"""sk2bGrow — statistics layer.

The Rust crates turn reference genomes into a deterministic 16-enzyme anchor
library and reads into a per-anchor count table. Everything from that table to a
PTR estimate lives here:

    counts.tsv ──ztp──> window rates ──gc_bias──> corrected rates
                                                       │
                              per-enzyme ──fit──> log2(PTR) ± se
                                                       │
                                            fusion ──> output.tsv

Module map, against the design report:

===============  ==========================================================
:mod:`io`        read/write the Rust interface files
:mod:`ztp`       zero-truncated Poisson / negative-binomial window rates (§7.1 step 3)
:mod:`gc_bias`   per-enzyme loess GC correction (fixes report defect D6)
:mod:`fit`       V-shape MLE on real coordinates, plus Pilea-parity sorted
                 regression (fixes D3)
:mod:`fusion`    inverse-variance fusion across enzymes + chi-square QC (fixes D4)
:mod:`dynamics`  multi-sample delta-PTR
:mod:`report`    output.tsv and QC figures
:mod:`simulate`  Monte-Carlo harness reproducing report §5
===============  ==========================================================
"""

__version__ = "0.1.0"

__all__ = [
    "io",
    "ztp",
    "gc_bias",
    "fit",
    "fusion",
    "dynamics",
    "report",
    "simulate",
]
