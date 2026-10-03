# Low-depth V-fit calibration prototype

Generated: 2026-10-03  
Git commit: `56d0e36df21747424f8aeb4ffcbfe946e4709750`  
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
| 0.5× baseline | r=0.8864 (95% bootstrap CI 0.7122–0.9773), RMSE=0.4382 (0.3246–0.5349), MAE=0.3786, R²=-0.0878 |
| 1× baseline | r=0.9489 (95% bootstrap CI 0.8855–0.9832), RMSE=0.3146 (0.2454–0.3778), MAE=0.2849, R²=0.4393 |
| 2× baseline | r=0.9586 (95% bootstrap CI 0.9142–0.9844), RMSE=0.1955 (0.1359–0.2461), MAE=0.1632, R²=0.7835 |
| Ridge, fixed α=1 | r=0.9504 (95% bootstrap CI 0.8568–0.9876), RMSE=0.1356 (0.0770–0.1952), MAE=0.1022, R²=0.8959 |
| Ridge, nested α | r=0.9368 (95% bootstrap CI 0.8311–0.9844), RMSE=0.1473 (0.1003–0.1979), MAE=0.1227, R²=0.8771 |

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
only from Zheng nested LOOMO (α=3), then applied without
further fitting to C1b.

| model | C1b 10× result |
|---|---|
| 2× baseline | r=0.9040 (95% bootstrap CI 0.7918–0.9618), RMSE=0.1543 (0.1130–0.1911), MAE=0.1253, R²=0.4637 |
| Zheng-trained ridge | r=0.8595 (95% bootstrap CI 0.7697–0.9321), RMSE=0.1736 (0.1368–0.2076), MAE=0.1529, R²=0.3209 |

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
