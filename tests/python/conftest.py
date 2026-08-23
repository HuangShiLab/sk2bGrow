"""Shared fixtures. Adds ``python/`` to the path so tests run from a checkout
without installing the package."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from sk2bgrow import io as sk_io  # noqa: E402

GENOME_LEN = 4_641_652
ORI = 3_923_883
ENZYMES = ["BcgI", "AlfI", "CjeI", "Hin4I"]


def make_counts(
    log2_ptr: float = 1.0,
    per_anchor_depth: float = 8.0,
    n_per_enzyme: int = 1200,
    enzymes=tuple(ENZYMES),
    sample: str = "S1",
    genome: str = "ecoli",
    sigma_eff: float = 0.1,
    gc_slope: float = 0.0,
    seed: int = 0,
) -> pd.DataFrame:
    """A synthetic count table with the columns the Rust layer writes.

    Anchors carry a real replication gradient around ``ORI``, so a correct
    pipeline recovers ``log2_ptr``.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for enzyme in enzymes:
        pos = np.sort(rng.integers(0, GENOME_LEN, size=n_per_enzyme))
        wrapped = np.mod(pos - ORI, GENOME_LEN)
        d = np.minimum(wrapped, GENOME_LEN - wrapped)
        log2_copy = -log2_ptr * d / (GENOME_LEN / 2)
        gc = rng.uniform(0.40, 0.60, size=n_per_enzyme)
        eff = rng.lognormal(-0.5 * sigma_eff**2, sigma_eff, size=n_per_enzyme) if sigma_eff > 0 else 1.0
        lam = per_anchor_depth * np.exp2(log2_copy + gc_slope * (gc - 0.5)) * eff
        rows.append(
            pd.DataFrame(
                {
                    "sample": sample,
                    "genome_id": 0,
                    "genome": genome,
                    "contig_id": 0,
                    "position": pos,
                    "global_position": pos,
                    "enzyme": enzyme,
                    "strand": "+",
                    "flags": 3,  # unique in genome and across the database
                    "local_gc": np.round(gc * 200).astype(int),
                    "window_id": pos // 25_000,
                    "count": rng.poisson(lam),
                }
            )
        )
    df = pd.concat(rows, ignore_index=True)
    df["gc"] = np.where(df["local_gc"] == sk_io.GC_UNDEFINED, np.nan, df["local_gc"] / 200.0)
    df["usable"] = (df["flags"] & sk_io.FLAG_USABLE_MASK) == 0
    return df


@pytest.fixture
def manifest() -> sk_io.DbManifest:
    return sk_io.DbManifest(
        enzymes=list(ENZYMES),
        gc_flank=250,
        genomes={
            0: sk_io.GenomeInfo(
                id=0, name="ecoli", genome_len=GENOME_LEN, n_contigs=1,
                taxonomy="d__Bacteria;s__Escherichia coli", ori=None, ori_confidence=0.0,
            )
        },
        version="0.1.0",
    )


@pytest.fixture
def manifest_with_ori(manifest) -> sk_io.DbManifest:
    g = manifest.genomes[0]
    return sk_io.DbManifest(
        enzymes=manifest.enzymes,
        gc_flank=manifest.gc_flank,
        genomes={0: sk_io.GenomeInfo(g.id, g.name, g.genome_len, g.n_contigs, g.taxonomy, ORI, 1.0)},
        version=manifest.version,
    )


@pytest.fixture
def fragmented_manifest(manifest) -> sk_io.DbManifest:
    g = manifest.genomes[0]
    return sk_io.DbManifest(
        enzymes=manifest.enzymes,
        gc_flank=manifest.gc_flank,
        genomes={0: sk_io.GenomeInfo(g.id, g.name, g.genome_len, 350, g.taxonomy, None, 0.0)},
        version=manifest.version,
    )
