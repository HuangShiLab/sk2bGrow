#!/usr/bin/env python3
"""Zheng low-depth -> 10x V-fit calibration prototype.

This is a small, transparent ridge-regression prototype, not a production model.
It uses leave-one-medium-out evaluation and a nested inner loop for alpha choice.
A C1b cross-species transfer test is included as an external warning benchmark.
"""
from __future__ import annotations

import subprocess
from collections import Counter
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, t

ROOT = Path(__file__).resolve().parents[2]
ZHENG = ROOT / "benches/refresh_20260930/zheng"
C1B = ROOT / "benches/realcommunity_20261001/remote_summaries/c1b_results_raw.tsv"
OUT = Path(__file__).resolve().parent
OUT.mkdir(parents=True, exist_ok=True)

SEEDS = ["s0", "s1", "s2"]
LOW_DEPTHS = [0.5, 1.0, 2.0]
ARM_A = "A"
ALPHAS = [0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0]
A_FEATURES = [
    f"A_{depth:g}_{stat}"
    for depth in LOW_DEPTHS
    for stat in ["mean", "std", "estcov", "pass"]
] + [
    "A_slope",
    "A_range",
    "A_last_drop",
    "A_linear_extrap10",
    "A_quad_extrap10",
]
TRANSFER_FEATURES = []
for depth in LOW_DEPTHS:
    TRANSFER_FEATURES.extend([f"A_{depth:g}_mean", f"A_{depth:g}_pass"])
TRANSFER_FEATURES.extend(
    [
        "A_slope",
        "A_range",
        "A_last_drop",
        "A_linear_extrap10",
        "A_quad_extrap10",
    ]
)

mpl.use("Agg")
mpl.rcParams.update(
    {
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "legend.fontsize": 7,
        "figure.dpi": 180,
        "savefig.dpi": 300,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def add_depth_features(row: dict) -> dict:
    x = np.log(LOW_DEPTHS)
    y = np.array([row[f"A_{depth:g}_mean"] for depth in LOW_DEPTHS], dtype=float)
    linear = np.polyfit(x, y, 1)
    quadratic = np.polyfit(x, y, 2)
    row["A_slope"] = linear[0]
    row["A_range"] = float(y.max() - y.min())
    row["A_last_drop"] = float(y[0] - y[-1])
    row["A_linear_extrap10"] = float(np.polyval(linear, np.log(10.0)))
    row["A_quad_extrap10"] = float(np.polyval(quadratic, np.log(10.0)))
    return row


def load_zheng() -> pd.DataFrame:
    frames = []
    for seed in SEEDS:
        path = ZHENG / seed / "results_raw.tsv"
        frame = pd.read_csv(path, sep="\t")
        frame = frame[frame["medium"] != "RUN_OUT"].copy()
        frames.append(frame)
    raw = pd.concat(frames, ignore_index=True)

    rows = []
    for medium, group in raw.groupby("medium", sort=True):
        row = {
            "medium": medium,
            "growth_rate": float(group["growth_rate"].iloc[0]),
            "y10": float(
                group.loc[
                    (group["arm"] == ARM_A) & (group["cov"] == 10.0), "log2ptr"
                ].mean()
            ),
        }
        for arm in ["A", "B", "E", "C_relaxed"]:
            for depth in LOW_DEPTHS:
                selected = group[(group["arm"] == arm) & (group["cov"] == depth)]
                prefix = f"{arm}_{depth:g}"
                row[f"{prefix}_mean"] = float(selected["log2ptr"].mean())
                row[f"{prefix}_std"] = float(selected["log2ptr"].std(ddof=1))
                row[f"{prefix}_estcov"] = float(selected["est_cov"].mean())
                row[f"{prefix}_pass"] = float(selected["passed"].mean())
        rows.append(add_depth_features(row))
    return pd.DataFrame(rows).sort_values("medium").reset_index(drop=True)


def load_c1b() -> pd.DataFrame:
    raw = pd.read_csv(C1B, sep="\t")
    raw = raw[raw["arm"] == "A_sk2bgrow"].copy()
    rows = []
    for run, group in raw.groupby("run", sort=True):
        row = {
            "run": run,
            "sample": group["sample"].iloc[0],
            "species": group["species"].iloc[0],
            "growth_rate": float(group["mu_per_h"].iloc[0]),
            "excluded_negative": bool(group["excluded_negative"].iloc[0]),
            "y10": float(group.loc[group["depth"] == 10.0, "log2ptr"].iloc[0]),
        }
        for depth in LOW_DEPTHS:
            selected = group[group["depth"] == depth]
            row[f"A_{depth:g}_mean"] = float(selected["log2ptr"].iloc[0])
            row[f"A_{depth:g}_pass"] = float(selected["pass_qc"].mean())
        rows.append(add_depth_features(row))
    return pd.DataFrame(rows).sort_values(["species", "run"]).reset_index(drop=True)


def ridge_fit(x: np.ndarray, y: np.ndarray, alpha: float) -> tuple:
    mean = x.mean(axis=0)
    scale = x.std(axis=0)
    scale[scale < 1e-8] = 1.0
    standardized = (x - mean) / scale
    target_mean = y.mean()
    matrix = standardized.T @ standardized + alpha * np.eye(x.shape[1])
    vector = standardized.T @ (y - target_mean)
    beta = np.linalg.solve(matrix, vector)
    return mean, scale, target_mean, beta


def ridge_predict(model: tuple, x: np.ndarray) -> np.ndarray:
    mean, scale, target_mean, beta = model
    return target_mean + ((x - mean) / scale) @ beta


def point_metrics(y: np.ndarray, prediction: np.ndarray) -> dict[str, float]:
    return {
        "pearson_r": float(pearsonr(y, prediction).statistic),
        "spearman_rho": float(spearmanr(y, prediction).statistic),
        "rmse": float(np.sqrt(np.mean((y - prediction) ** 2))),
        "mae": float(np.mean(np.abs(y - prediction))),
        "r2": float(
            1
            - np.sum((y - prediction) ** 2)
            / np.sum((y - y.mean()) ** 2)
        ),
    }


def bootstrap_metrics(
    y: np.ndarray,
    prediction: np.ndarray,
    n_boot: int = 10000,
    seed: int = 20261003,
) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    n = len(y)
    pearsons, rmses = [], []
    for _ in range(n_boot):
        index = rng.integers(0, n, n)
        yi, pi = y[index], prediction[index]
        if np.std(yi) == 0 or np.std(pi) == 0:
            continue
        pearsons.append(pearsonr(yi, pi).statistic)
        rmses.append(np.sqrt(np.mean((yi - pi) ** 2)))
    return {
        "pearson_ci_low": float(np.quantile(pearsons, 0.025)),
        "pearson_ci_high": float(np.quantile(pearsons, 0.975)),
        "rmse_ci_low": float(np.quantile(rmses, 0.025)),
        "rmse_ci_high": float(np.quantile(rmses, 0.975)),
    }


def metrics_row(
    evaluation_target: str,
    model: str,
    y: np.ndarray,
    prediction: np.ndarray,
    n: int | None = None,
) -> dict:
    values = point_metrics(y, prediction)
    values.update(bootstrap_metrics(y, prediction))
    values.update({"evaluation_target": evaluation_target, "model": model})
    values["n"] = int(n if n is not None else len(y))
    return values


def fixed_ridge_loomo(
    frame: pd.DataFrame, features: list[str], target: str, alpha: float
) -> np.ndarray:
    x = frame[features].to_numpy(float)
    y = frame[target].to_numpy(float)
    predictions = []
    for index in range(len(frame)):
        train = np.arange(len(frame)) != index
        model = ridge_fit(x[train], y[train], alpha)
        predictions.append(ridge_predict(model, x[index]))
    return np.asarray(predictions)


def nested_ridge_loomo(
    frame: pd.DataFrame, features: list[str], target: str
) -> pd.DataFrame:
    x = frame[features].to_numpy(float)
    y = frame[target].to_numpy(float)
    records = []
    for index in range(len(frame)):
        outer_train = np.arange(len(frame)) != index
        scores = {}
        for alpha in ALPHAS:
            squared_errors = []
            for inner_index in np.where(outer_train)[0]:
                inner_train = outer_train & (np.arange(len(frame)) != inner_index)
                model = ridge_fit(x[inner_train], y[inner_train], alpha)
                prediction = ridge_predict(model, x[inner_index])
                squared_errors.append((y[inner_index] - prediction) ** 2)
            scores[alpha] = float(np.mean(squared_errors))
        selected_alpha = min(ALPHAS, key=lambda alpha: scores[alpha])
        model = ridge_fit(x[outer_train], y[outer_train], selected_alpha)
        prediction = float(ridge_predict(model, x[index]))

        inner_residuals = []
        for inner_index in np.where(outer_train)[0]:
            inner_train = outer_train & (np.arange(len(frame)) != inner_index)
            inner_model = ridge_fit(
                x[inner_train], y[inner_train], selected_alpha
            )
            inner_prediction = ridge_predict(inner_model, x[inner_index])
            inner_residuals.append(y[inner_index] - inner_prediction)
        residual_sd = float(np.sqrt(np.mean(np.square(inner_residuals))))
        interval_halfwidth = float(
            t.ppf(0.975, len(inner_residuals) - 1) * residual_sd
        )
        records.append(
            {
                "prediction": prediction,
                "alpha": selected_alpha,
                "residual_sd": residual_sd,
                "ci_low": prediction - interval_halfwidth,
                "ci_high": prediction + interval_halfwidth,
            }
        )
    return pd.DataFrame(records)


def nested_alpha_mode(frame: pd.DataFrame, features: list[str], target: str) -> float:
    result = nested_ridge_loomo(frame, features, target)
    return float(Counter(result["alpha"]).most_common(1)[0][0])


def make_figure(
    zheng: pd.DataFrame, predictions: pd.DataFrame, c1b: pd.DataFrame, c1b_predictions: pd.DataFrame
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.15), constrained_layout=True)
    ax = axes[0]
    ax.plot(
        zheng["y10"],
        predictions["baseline_2x"],
        "o",
        ms=4,
        color="#8c8c8c",
        label="2× baseline",
    )
    ax.plot(
        zheng["y10"],
        predictions["ridge_fixed_alpha1"],
        "o",
        ms=4,
        color="#1f77b4",
        label="ridge α=1",
    )
    limits = [min(zheng.y10.min(), predictions[["baseline_2x", "ridge_fixed_alpha1"]].min().min()) - 0.05,
              max(zheng.y10.max(), predictions[["baseline_2x", "ridge_fixed_alpha1"]].max().max()) + 0.05]
    ax.plot(limits, limits, "--", color="#444444", lw=0.8)
    ax.set(xlabel="Observed 10× V-fit", ylabel="Predicted 10× V-fit",
           title="A  Zheng leave-one-medium-out", xlim=limits, ylim=limits)
    ax.legend(frameon=False, loc="upper left")

    ax = axes[1]
    x = np.arange(len(zheng))
    ax.bar(x - 0.2, abs(predictions["baseline_2x"] - zheng["y10"]), 0.38,
           label="2× baseline", color="#8c8c8c")
    ax.bar(x + 0.2, abs(predictions["ridge_fixed_alpha1"] - zheng["y10"]), 0.38,
           label="ridge α=1", color="#1f77b4")
    ax.set(xticks=x[::2], xlabel="Medium index", ylabel="Absolute error (log2PTR)",
           title="B  Zheng absolute error")
    ax.legend(frameon=False)

    ax = axes[2]
    species_colors = dict(zip(sorted(c1b.species.unique()), mpl.cm.tab10.colors))
    for species, group in c1b_predictions.groupby("species"):
        color = species_colors[species]
        ax.plot(group.observed_10x_vfit, group.baseline_2x, "o", ms=4, color=color, alpha=0.55)
        ax.plot(group.observed_10x_vfit, group.ridge_transfer, "s", ms=4, color=color)
    limits = [
        min(c1b.y10.min(), c1b_predictions[["baseline_2x", "ridge_transfer"]].min().min()) - 0.1,
        max(c1b.y10.max(), c1b_predictions[["baseline_2x", "ridge_transfer"]].max().max()) + 0.1,
    ]
    ax.plot(limits, limits, "--", color="#444444", lw=0.8)
    ax.set(xlabel="Observed C1b 10× V-fit", ylabel="Zheng-trained prediction",
           title="C  C1b external transfer", xlim=limits, ylim=limits)
    ax.text(0.03, 0.97, "circle: 2× baseline\nsquare: ridge", transform=ax.transAxes,
            va="top", fontsize=7)
    save(fig, "vfit_lowdepth_prototype")


def save(fig, stem: str) -> None:
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"{stem}.{suffix}", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    zheng = load_zheng()
    c1b = load_c1b()

    # Zheng depth-to-depth calibration: 0.5/1/2x -> 10x all-finite V-fit.
    nested = nested_ridge_loomo(zheng, A_FEATURES, "y10")
    predictions = pd.DataFrame(
        {
            "medium": zheng["medium"],
            "growth_rate": zheng["growth_rate"],
            "observed_10x_vfit": zheng["y10"],
            "baseline_0.5x": zheng["A_0.5_mean"],
            "baseline_1x": zheng["A_1_mean"],
            "baseline_2x": zheng["A_2_mean"],
            "ridge_fixed_alpha1": fixed_ridge_loomo(zheng, A_FEATURES, "y10", 1.0),
            "ridge_nested": nested["prediction"],
            "ridge_alpha": nested["alpha"],
            "ridge_residual_sd": nested["residual_sd"],
            "ridge_ci_low": nested["ci_low"],
            "ridge_ci_high": nested["ci_high"],
        }
    )
    predictions.to_csv(OUT / "zheng_loomo_predictions.tsv", sep="\t", index=False)

    metric_rows = []
    for model in [
        "baseline_0.5x",
        "baseline_1x",
        "baseline_2x",
        "ridge_fixed_alpha1",
        "ridge_nested",
    ]:
        metric_rows.append(
            metrics_row(
                "high_depth_10x_vfit",
                model,
                zheng["y10"].to_numpy(float),
                predictions[model].to_numpy(float),
            )
        )
        metric_rows.append(
            metrics_row(
                "measured_growth_rate_posthoc",
                model,
                zheng["growth_rate"].to_numpy(float),
                predictions[model].to_numpy(float),
            )
        )
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(OUT / "zheng_metrics.tsv", sep="\t", index=False)

    # External C1b transfer. Alpha is selected only from Zheng nested LOOMO.
    transfer_alpha = nested_alpha_mode(zheng, TRANSFER_FEATURES, "y10")
    x_train = zheng[TRANSFER_FEATURES].to_numpy(float)
    y_train = zheng["y10"].to_numpy(float)
    transfer_model = ridge_fit(x_train, y_train, transfer_alpha)
    x_c1b = c1b[TRANSFER_FEATURES].to_numpy(float)
    c1b_predictions = pd.DataFrame(
        {
            "run": c1b["run"],
            "sample": c1b["sample"],
            "species": c1b["species"],
            "excluded_negative": c1b["excluded_negative"],
            "observed_10x_vfit": c1b["y10"],
            "growth_rate": c1b["growth_rate"],
            "baseline_2x": c1b["A_2_mean"],
            "ridge_transfer": ridge_predict(transfer_model, x_c1b),
            "transfer_alpha": transfer_alpha,
        }
    )
    c1b_predictions.to_csv(OUT / "c1b_transfer_predictions.tsv", sep="\t", index=False)

    c1b_metric_rows = [
        metrics_row(
            "high_depth_10x_vfit",
            "baseline_2x",
            c1b["y10"].to_numpy(float),
            c1b_predictions["baseline_2x"].to_numpy(float),
        ),
        metrics_row(
            "high_depth_10x_vfit",
            f"zheng_ridge_alpha_{transfer_alpha:g}",
            c1b["y10"].to_numpy(float),
            c1b_predictions["ridge_transfer"].to_numpy(float),
        ),
    ]
    nonnegative = ~c1b["excluded_negative"].to_numpy(bool)
    c1b_metric_rows.extend(
        [
            metrics_row(
                "mu_per_h_nonnegative_posthoc",
                "baseline_2x",
                c1b.loc[nonnegative, "growth_rate"].to_numpy(float),
                c1b_predictions.loc[nonnegative, "baseline_2x"].to_numpy(float),
            ),
            metrics_row(
                "mu_per_h_nonnegative_posthoc",
                f"zheng_ridge_alpha_{transfer_alpha:g}",
                c1b.loc[nonnegative, "growth_rate"].to_numpy(float),
                c1b_predictions.loc[nonnegative, "ridge_transfer"].to_numpy(float),
            ),
        ]
    )
    c1b_metrics = pd.DataFrame(c1b_metric_rows)
    c1b_metrics.to_csv(OUT / "c1b_transfer_metrics.tsv", sep="\t", index=False)

    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:
        commit = "unavailable"

    zheng_y10 = metrics[metrics.evaluation_target == "high_depth_10x_vfit"]
    c1b_y10 = c1b_metrics[c1b_metrics.evaluation_target == "high_depth_10x_vfit"]

    def fmt(model: str, table: pd.DataFrame) -> str:
        row = table[table.model == model].iloc[0]
        return (
            f"r={row.pearson_r:.4f} "
            f"(95% bootstrap CI {row.pearson_ci_low:.4f}–{row.pearson_ci_high:.4f}), "
            f"RMSE={row.rmse:.4f} ({row.rmse_ci_low:.4f}–{row.rmse_ci_high:.4f}), "
            f"MAE={row.mae:.4f}, R²={row.r2:.4f}"
        )

    report = f"""# Low-depth V-fit calibration prototype

Generated: 2026-10-03  
Git commit: `{commit}`  
Script: `benches/prototype_20261003/run_vfit_calibration.py`

## Purpose

This is a deliberately minimal prototype for learning a depth-to-depth
calibration from Zheng 0.5/1/2× profiles to the mean 10× all-finite V-fit
estimate. It is **not** a production estimator and does not claim true PTR
accuracy.

## Design

- Training/inference unit: one Zheng growing medium.
- Technical seeds s0/s1/s2 were aggregated within each medium; they were not
  treated as independent biological replicates.
- RUN_OUT was excluded, matching the manuscript's correlation benchmark.
- Outer validation: leave-one-medium-out (LOOMO), n=16.
- Model: ridge regression with standardized A-arm features from 0.5/1/2× plus
  low-depth trend features.
- Alpha was selected inside each outer training fold by nested LOOMO.
- Target: mean A-arm log2PTR at 10×. This is a stable-reference prediction,
  not independent biological truth.

## Zheng internal LOOMO results

| model | result |
|---|---|
| 0.5× baseline | {fmt('baseline_0.5x', zheng_y10)} |
| 1× baseline | {fmt('baseline_1x', zheng_y10)} |
| 2× baseline | {fmt('baseline_2x', zheng_y10)} |
| Ridge, fixed α=1 | {fmt('ridge_fixed_alpha1', zheng_y10)} |
| Ridge, nested α | {fmt('ridge_nested', zheng_y10)} |

Under the strict no-5×-feature rule, the fixed α=1 ridge improves squared-error
calibration relative to the 2× baseline: RMSE falls from 0.1955 to 0.1356 and
MAE from 0.1632 to 0.1022, while R² rises from 0.784 to 0.896. However, its
Pearson r is slightly lower (0.950 versus 0.959), and the nested-alpha choice
does not improve rank accuracy. Therefore this prototype supports only modest
error calibration, not a claim of universally better low-depth ranking. The
earlier 5×→10× setting is already near the estimator ceiling and leaves little
room for improvement.

## C1b external transfer warning

A transfer-compatible feature set (omitting seed SD and estimated coverage,
which are unavailable per C1b run) was trained on Zheng. Alpha was selected
only from Zheng nested LOOMO (α={transfer_alpha:g}), then applied without
further fitting to C1b.

| model | C1b 10× result |
|---|---|
| 2× baseline | {fmt('baseline_2x', c1b_y10)} |
| Zheng-trained ridge | {fmt(f'zheng_ridge_alpha_{transfer_alpha:g}', c1b_y10)} |

The Zheng-calibrated ridge does **not** transfer to the four C1b species. This
is the expected warning result: the current E. coli prototype learns a
condition- and species-specific depth correction, not a universal law.

## Interpretation limits

1. n=16 media is sufficient only for a prototype.
2. Zheng subsample seeds are technical, not biological replicates.
3. The 10× target is an estimator reference, not independent PTR truth.
4. A linear ridge model uses shallow window summaries, not raw anchor counts.
5. C1b transfer failure shows that cross-species training data are required.
6. The model must report calibrated intervals and QC before deployment use.

## Files

- `zheng_loomo_predictions.tsv`
- `zheng_metrics.tsv`
- `c1b_transfer_predictions.tsv`
- `c1b_transfer_metrics.tsv`
- `vfit_lowdepth_prototype.png`
- `vfit_lowdepth_prototype.pdf`
"""
    (OUT / "README.md").write_text(report, encoding="utf-8")

    make_figure(zheng, predictions, c1b, c1b_predictions)
    print(OUT)
    print(metrics[(metrics.evaluation_target == "high_depth_10x_vfit")][
        ["model", "pearson_r", "rmse", "mae", "r2"]
    ].to_string(index=False))
    print("\nC1b transfer")
    print(c1b_y10[["model", "pearson_r", "rmse", "mae", "r2"]].to_string(index=False))


if __name__ == "__main__":
    main()
