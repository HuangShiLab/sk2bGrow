# Supplementary material

## S1 Software provenance

| item | value |
|---|---|
| repository | sk2bGrow, current implementation after signed fixed-origin correction |
| Zheng source archive | `sk2bgrow-source-20261001.tar.gz` |
| source SHA256 | `8a157b800887df43d5417b901e229fc302554418c119a682499f129ea01c92a7` |
| Zheng result archive | `sk2bgrow-refresh-results-20260930.tar.gz` |
| result SHA256 | `4b5840dff3681aa0ca3c639e48d26c0cbcf4de48bc20eae6d7f3404da25315d1` |
| Pilea | vendor `pilea138` environment, v1.3.8 |
| figure generator | `benches/refresh_20260930/make_paper_figures.py` |
| review diagnostics | `benches/refresh_20260930/review_checks.py` |

The final mixed-strain validation also records the current Git commit and the
16 reference FASTA SHA256 hashes in
`benches/mixedstrain_20261001/results/provenance.tsv`.

Real-community results are separately labelled below. They were run before the
current implementation refresh and are not described as current-code reruns.

## S2 Zheng benchmark design

* BioProject: PRJNA615952.
* Conditions: 16 growing media plus RUN_OUT.
* Reads: paired sequencing subsamples.
* Subsample seeds: `s0` archived C1 subsample; `s1`, `s2` deterministic paired
  resamples from the archived C1 20x pool.
* Depths: 0.5x, 1x, 2x, 5x, 10x.
* Completed profile outputs: 1,275.

### S2.1 Arms

| arm | sketch | estimator |
|---|---|---|
| A | 16-enzyme anchors | adaptive windows + coordinate V-fit |
| E | FracMinHash | adaptive windows + coordinate V-fit |
| B | 16-enzyme anchors | 25 kb windows + sorted/RANSAC |
| C relaxed | FracMinHash | Pilea gates-off |
| C default | FracMinHash | Pilea shipped defaults |

## S3 Zheng statistics

Metrics average the three seed estimates within each growing medium and then
compare across media. The RUN_OUT condition is used as a negative control and
is excluded from correlations. Pearson and Spearman correlations use measured
growth rate. RMSE uses the Zheng predicted log2 PTR and is not independent
absolute accuracy.

### S3.1 Median bootstrap intervals

`results/zheng_bootstrap_media.tsv` contains point estimates and 10,000-resample
medium-level 95% intervals. The bootstrap resamples media; seed estimates are
averaged within each selected medium. Therefore the interval does not represent
biological replication.

### S3.2 QC instability

`results/zheng_qc_pass_counts.tsv` reports exact pass counts per seed and the
number of media passing in all or any seed. Mean pass counts are not equivalent
to reproducible deployment.

### S3.3 Paired A-minus-E diagnostics

`results/zheng_arm_delta_bootstrap.tsv` reports deltas for Pearson, Spearman,
RMSE and slope. Exploratory bootstrap p values are unadjusted. The primary
inference uses intervals and effect direction, not significance testing.

## S4 A4 window-rate statistics

A4 simulations used real E. coli K-12 coordinates and known tent amplitudes.
Depths were 0.5/1/2/5x; true log2 PTR values were 0/0.5/1/1.5/2; there were
three seeds. The Poisson arm used integer Poisson counts. The NB arm applied
per-kilobase lognormal efficiencies with sigma 0.5. Exact expected counts were
reconstructed from the simulator and joined to production windows.

`results/a4_slope_bootstrap.tsv` contains slope and b=0 intervals from 10,000
replicate bootstrap resamples. Replicates are technical simulation seeds.

## S5 R3 coordinate-QC statistics

R3 evaluated 27 references across complete, ordered and scrambled coordinate
families, 12 read cells and 324 profiles. The primary statistic used two
branches:

* slope branch: within-contig amplitude versus global fitted amplitude;
* jump branch: excess squared rate jump across contig boundaries.

The primary threshold used `z_short > 2.0`. The review sensitivity table varied
this threshold from 2.0 to 3.0. A post-hoc threshold of 2.5 removed the observed
stationary false fire while retaining all strong 5–10x detections. Independent
validation is required before production use.

## S6 Mixed-strain simulation

The mixed-strain validation used 16 bacterial reference genomes. Communities
contained 4, 8 or 16 strains. Selected genomes received 1, 2 or 4x coverage per
genome, and each cell had two simulation replicates. Strain-specific true log2
PTR values were drawn from Uniform(0,2). Reads followed a V-shaped origin-terminus
profile.

Reported metrics:

* recall: selected true genomes with a finite estimate / all selected genomes;
* RMSE: root mean squared error of matched log2 PTR estimates;
* bias: mean signed matched error;
* spurious: reported genomes absent from truth.

## S7 Real-community benchmarks

### S7.1 C1b cross-species isolates

* BioProject: PRJNA1280254.
* Design: four wastewater isolate species (_Bacillus subtilis_,
  _Klebsiella pneumoniae_, _Morganella morganii_ and _Pseudomonas putida_)
  across five nutrient concentrations.
* Subsampling: 0.5, 1, 2, 5 and 10x species-assembly coverage.
* Protocols: sk2bGrow fixed coordinate estimator; Pilea v1.3.8 default;
  Pilea gates-off (`-x 0 -z 0 -c 0`).
* Truth: digitized from the source article's growth-rate figure and validated
  against printed panel correlations to within 2e-4.
* Caveat: RMSE against the paper's predicted log2 PTR is a same-reads
  comparison, not independent absolute accuracy.

### S7.2 C4 marine metagenomes

* BioProject: PRJNA551656.
* Design: 20 marine surface-water metagenomes and 101 marine MAGs.
* Protocols: sk2bGrow and Pilea default; Pilea gates-off was recorded but had
  16/20 sample crashes and no balanced protocol correlation set.
* Interpretation: a MAG is protocol-comparable only when both protocols have
  finite matched samples. Growth-rate truth from Long et al. was derived from
  the same reads, so C4 is agreement analysis rather than absolute accuracy.

### S7.3 C5 rotating biological contactor application

* BioProject: PRJNA974210.
* Design: nine metagenomes and 522 bacterial MAGs (4,698 expected MAG-sample
  cells).
* Outcomes: estimate recall, default-QC passes, suspicious absent-genome
  confident estimates, cross-sample consistency and wall time.
* Interpretation: C5 tests application-scale recall and workflow coverage. It
  has no independent growth-rate truth and does not establish accuracy.

### S7.4 Real-community provenance caveat

The C1b, C4 and C5 runs used sk2bGrow commit
`275778f350b87e10c6366bf90884d965fbab45a6` plus working-tree fixes to
`python/sk2bgrow/fit.py` and `python/sk2bgrow/ztp.py`, with Pilea v1.3.8. They
were not rerun at the current implementation commit.

## S8 Data availability

* Zheng 2020 reads: BioProject PRJNA615952.
* C1b reads: BioProject PRJNA1280254.
* C4 reads: BioProject PRJNA551656.
* C5 reads: BioProject PRJNA974210.
* Reference panel: 16 bacterial FASTA files included with benchmark scripts and
  checksummed in the mixed-strain provenance table.
* Implementation: sk2bGrow repository.
* Result archives: hashes listed in S1.
* Figure source: `benches/refresh_20260930/make_paper_figures.py`.
* Mixed-strain source: `benches/mixedstrain_20261001/run_validation.sh`.
* Review diagnostics: `benches/refresh_20260930/review_checks.py`.

A citable archival DOI should be generated for the final source and result
archives on Zenodo or another DOI-issuing repository before submission.

## S9 Code availability

The code is available under the repository LICENSE. The exact source archive
hash and result archive hash are part of this supplement. The final commit hash
should be recorded beside these hashes when the manuscript is submitted.

## S10 Claim-consistency rules

1. Use "all finite" for estimator signal and "default QC passed" for deployment.
2. Do not describe the three sequencing seeds as biological replicates.
3. Do not describe a Pilea comparison as a uniform accuracy win.
4. Do not describe current PTR as EM-reassigned.
5. Do not describe WCG as production QC.
6. State Pilea v1.3.8, archive hashes and the exact benchmark arms.
7. Do not infer real-community accuracy from C5 because C5 has no independent
   growth-rate truth.
8. Keep C1b and C4 negative or boundary comparisons visible in abstract,
   results and conclusions.

## S11 Supplementary tables

Machine-readable TSVs are preferred for analysis. Human-readable versions of
the key result tables are in `docs/paper/supplementary_tables.md`.
