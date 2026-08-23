"""Reading and writing the files that join the Rust and Python layers.

The two layers share no state: the Rust counter writes a count table, this layer
reads it. Keeping the boundary at a file (rather than an FFI call) is what lets
the statistics iterate without rebuilding the counter, and it makes every
intermediate inspectable with ``head``.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

#: Columns written by ``sk2bgrow_core::count::write_count_table``.
COUNT_COLUMNS = [
    "sample",
    "genome_id",
    "genome",
    "contig_id",
    "position",
    "global_position",
    "enzyme",
    "strand",
    "flags",
    "local_gc",
    "window_id",
    "count",
]

WINDOW_COLUMNS = [
    "window_id",
    "genome_id",
    "genome",
    "contig_id",
    "start",
    "end",
    "global_mid",
    "span",
    "n_anchors",
]

#: Anchor flag bits, mirroring ``sk2bgrow_core::anchor_db::flags``.
FLAG_UNIQUE_IN_GENOME = 1 << 0
FLAG_UNIQUE_ACROSS_DB = 1 << 1
FLAG_MASKED_MULTICOPY = 1 << 2
FLAG_MASKED_SHARED = 1 << 3
FLAG_NON_CHROMOSOMAL = 1 << 4
FLAG_GC_UNDEFINED = 1 << 5
FLAG_USABLE_MASK = FLAG_MASKED_MULTICOPY | FLAG_MASKED_SHARED | FLAG_NON_CHROMOSOMAL

#: Sentinel written for anchors that were not assigned to a window.
NO_WINDOW = 0xFFFFFFFF

#: Sentinel stored in ``local_gc`` when the +/-250 bp window was all N.
GC_UNDEFINED = 255


class Sk2bIoError(RuntimeError):
    """Raised when an interface file is missing columns or self-inconsistent."""


def read_counts(path: str | Path) -> pd.DataFrame:
    """Read one sample's count table.

    ``local_gc`` is de-quantised back to a fraction and the undefined sentinel
    becomes ``NaN``, so downstream code never mistakes 255 for a 127.5 % GC
    window.
    """
    path = Path(path)
    df = pd.read_csv(path, sep="\t", dtype={"enzyme": str, "genome": str, "sample": str})
    missing = [c for c in COUNT_COLUMNS if c not in df.columns]
    if missing:
        raise Sk2bIoError(f"{path}: count table is missing columns {missing}")
    df["gc"] = np.where(df["local_gc"] == GC_UNDEFINED, np.nan, df["local_gc"] / 200.0)
    df["usable"] = (df["flags"] & FLAG_USABLE_MASK) == 0
    return df


def read_windows(path: str | Path) -> pd.DataFrame:
    """Read the window table emitted by ``sk2bgrow profile``."""
    path = Path(path)
    df = pd.read_csv(path, sep="\t")
    missing = [c for c in WINDOW_COLUMNS if c not in df.columns]
    if missing:
        raise Sk2bIoError(f"{path}: window table is missing columns {missing}")
    return df


@dataclass(frozen=True)
class GenomeInfo:
    """One genome's entry from the database manifest."""

    id: int
    name: str
    genome_len: int
    n_contigs: int
    taxonomy: str | None
    ori: int | None
    ori_confidence: float

    @property
    def is_contiguous(self) -> bool:
        """Whether real-coordinate fitting is defensible for this reference.

        Pilea's own assembly-quality analysis puts the boundary near 100 contigs
        (report §1.2); beyond it, contig order is guesswork and the V-shape
        x-axis is not trustworthy.
        """
        return self.n_contigs <= 100


@dataclass(frozen=True)
class DbManifest:
    """The parts of ``manifest.json`` the statistics layer needs."""

    enzymes: list[str]
    gc_flank: int
    genomes: dict[int, GenomeInfo]
    version: str

    def genome(self, genome_id: int) -> GenomeInfo:
        try:
            return self.genomes[int(genome_id)]
        except KeyError as exc:
            raise Sk2bIoError(f"genome_id {genome_id} is not in the database manifest") from exc


#: Panel order, mirroring ``sk2bgrow_core::enzyme::PANEL``. Index == enzyme_idx.
ENZYME_PANEL = [
    "BcgI", "AlfI", "AloI", "BaeI", "BplI", "BsaXI", "BslFI", "Bsp24I",
    "CjeI", "CjePI", "CspCI", "FalI", "HaeIV", "Hin4I", "PpiI", "PsrI",
]


def read_manifest(db_dir: str | Path) -> DbManifest:
    """Read ``<db>/manifest.json``."""
    path = Path(db_dir) / "manifest.json"
    try:
        raw = json.loads(path.read_text())
    except FileNotFoundError as exc:
        raise Sk2bIoError(f"{path} not found; build one with `sk2bgrow index`") from exc

    params = raw.get("params", {})
    bits = int(params.get("enzymes", 0))
    enzymes = [name for i, name in enumerate(ENZYME_PANEL) if bits & (1 << i)]

    genomes: dict[int, GenomeInfo] = {}
    for g in raw.get("genomes", []):
        genomes[int(g["id"])] = GenomeInfo(
            id=int(g["id"]),
            name=g["name"],
            genome_len=int(g["genome_len"]),
            n_contigs=len(g.get("contigs", [])),
            taxonomy=g.get("taxonomy"),
            ori=None if g.get("ori") is None else int(g["ori"]),
            ori_confidence=float(g.get("ori_confidence", 0.0)),
        )
    return DbManifest(
        enzymes=enzymes,
        gc_flank=int(params.get("gc_flank", 250)),
        genomes=genomes,
        version=str(params.get("sk2bgrow_version", "unknown")),
    )


def read_stats(path: str | Path) -> dict:
    """Read a ``<sample>.stats.json`` sidecar."""
    return json.loads(Path(path).read_text())


def stats_path_for(counts_path: str | Path) -> Path:
    """Locate the stats sidecar next to a count table."""
    p = Path(counts_path)
    return p.with_name(p.name.replace(".counts.tsv", ".stats.json"))


def concat_counts(paths: Iterable[str | Path]) -> pd.DataFrame:
    """Read several count tables into one long frame.

    Every sample must have been counted against the same database — anchor
    identity is what makes the anchor x sample matrix in :mod:`dynamics` valid,
    so a mismatch is an error rather than an outer join.
    """
    frames = []
    key = None
    for p in paths:
        df = read_counts(p)
        this_key = (
            int(df["genome_id"].nunique()),
            int(len(df)),
            int(df["position"].sum()),
        )
        if key is None:
            key = this_key
        elif this_key != key:
            raise Sk2bIoError(
                f"{p} was counted against a different anchor database "
                f"(anchor signature {this_key} vs {key}); re-run `sk2bgrow profile` "
                "for every sample against one database"
            )
        frames.append(df)
    if not frames:
        raise Sk2bIoError("no count tables given")
    return pd.concat(frames, ignore_index=True)
