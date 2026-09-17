#!/usr/bin/env python3
"""C7/A4 decomposition: where the genome-level slope compression comes from.

Runs after ``analyze.py`` (needs ``results/per_window.tsv``). All fits go
through the real pipeline machinery (``fit_windows`` + ``fuse_table``); one
factor is changed at a time:

  * ``recovered``          — the pipeline as shipped (IV weights, Tukey, ori search)
  * ``true_all_windows``   — exact rates everywhere         (fit machinery check)
  * ``true_same_windows``  — exact rates, same dropout      (window-dropout check)
  * ``const_se``           — recovered rates, flat weights  (self-weighting term)
  * ``ori_fixed``          — recovered rates, true origin   (ori-search term)
  * ``no_tukey``           — recovered rates, no trimming   (trimming term;
                             collapsed bug-artifact windows hand-removed, since
                             without trimming they destroy the fits outright)

Outputs ``results/decomposition.tsv`` and prints the ledger.
"""
from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))

from sk2bgrow import fit as sk_fit  # noqa: E402
from sk2bgrow import fusion as sk_fusion  # noqa: E402
from sk2bgrow import io as sk_io  # noqa: E402

ORI = 3_923_883.0

FIT_COLS = ["sample", "genome_id", "genome", "enzyme", "window", "global_mid",
            "rate", "se", "log2_rate", "log2_se", "detected_fraction", "dispersion",
            "n_anchors", "true_rate", "ratio", "log2_rec"]


def genome_slopes(windows: pd.DataFrame, manifest, meta: pd.DataFrame, **kw) -> pd.DataFrame:
    """OLS of fused log2PTR on true b, per (arm, depth), through the real fit+fusion."""
    pe = sk_fit.fit_windows(windows, manifest, method="v_shape", **kw)
    fu = sk_fusion.fuse_table(pe).merge(meta, on="sample")
    rows = []
    for (arm, depth), grp in fu.groupby(["arm", "depth"]):
        x = grp["true_log2_ptr"].to_numpy()
        y = grp["log2_ptr"].to_numpy()
        ok = np.isfinite(x) & np.isfinite(y)
        slope = float(np.polyfit(x[ok], y[ok], 1)[0]) if ok.sum() >= 5 else np.nan
        floor = float(grp.loc[grp["true_log2_ptr"] == 0, "log2_ptr"].mean())
        rows.append({"arm": arm, "depth": depth, "slope": slope, "b0_floor": floor, "n": int(ok.sum())})
    return pd.DataFrame(rows)


def main() -> None:
    root = Path(__file__).resolve().parent
    pw = pd.read_csv(root / "results" / "per_window.tsv", sep="\t", low_memory=False)
    pw = pw[pw["arm"].isin(["pois", "nb"])].copy()
    pw["log2_rec"] = np.where(pw["rate"] > 0, np.log2(pw["rate"]), np.nan)
    pw["ratio"] = pw["log2_rec"] - pw["log2_true"]

    manifest = sk_io.read_manifest(root / "work" / "db")
    gid, g = next(iter(manifest.genomes.items()))
    manifest_ori = dataclasses.replace(
        manifest, genomes={gid: dataclasses.replace(g, ori=int(ORI), ori_confidence=1.0)}
    )
    meta = pw[["sample", "arm", "depth", "true_log2_ptr"]].drop_duplicates("sample")

    variants = {}
    base = pw[FIT_COLS].copy()

    truth = base.copy()
    truth["log2_rate"] = np.where(truth["true_rate"] > 0, np.log2(truth["true_rate"]), np.nan)
    truth["rate"] = truth["true_rate"]
    med = truth.groupby(["sample", "enzyme"])["log2_se"].transform("median")
    truth["log2_se"] = truth["log2_se"].fillna(med).fillna(1.0)
    variants["true_all_windows"] = (truth, manifest, {})

    same = truth.copy()
    same.loc[~np.isfinite(base["log2_rec"]), "log2_rate"] = np.nan
    variants["true_same_windows"] = (same, manifest, {})

    variants["recovered"] = (base, manifest, {})

    const = base.copy()
    const["log2_se"] = 1.0
    variants["const_se"] = (const, manifest, {})
    variants["ori_fixed"] = (base, manifest_ori, {})
    variants["const_se+ori_fixed"] = (const, manifest_ori, {})

    # Trimming term: the collapsed all-ones/ZTNB windows (ratio < -3) are a
    # separate bug artifact; remove them by hand so the no-trimming arm is
    # interpretable, then toggle Tukey.
    noco = base[np.isfinite(base["ratio"]) & (base["ratio"] > -3)].copy()
    noco["log2_se"] = 1.0
    variants["const_se+tukey_off"] = (noco, manifest, {"tukey_k": 1e9})
    variants["const_se+ori_fixed+tukey_off"] = (noco, manifest_ori, {"tukey_k": 1e9})

    frames = []
    for name, (w, man, kw) in variants.items():
        df = genome_slopes(w, man, meta, **kw)
        df["variant"] = name
        frames.append(df)
        print(f"  {name} done", file=sys.stderr)
    out = pd.concat(frames, ignore_index=True)
    out.to_csv(root / "results" / "decomposition.tsv", sep="\t", index=False, float_format="%.10g")

    print("\ngenome-level slope, estimated log2PTR ~ true b (b0 floor = mean estimate at b=0):")
    piv = out.pivot_table(index=["variant"], columns=["arm", "depth"], values="slope")
    print(piv.round(3).to_string())
    print("\nb0 floors:")
    print(out.pivot_table(index=["variant"], columns=["arm", "depth"], values="b0_floor").round(3).to_string())


if __name__ == "__main__":
    main()
