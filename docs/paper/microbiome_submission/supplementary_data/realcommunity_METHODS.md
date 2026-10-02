# C1b — cross-species isolate accuracy: METHODS (running log)

Experiment C1b of the HPC phase: PRJNA1280254 (Pilea paper's own isolate
dataset; Chen et al., Microbiome 2026, PMC13123247). 20 runs, 4 species x 5 LB
nutrient concentrations (2x, 1x, 1/3x, 1/5x, 1/10x), subsampled to
0.5/1/2/5/10x per-species assembly coverage. sk2bGrow (fixed estimator) vs
Pilea v1.3.8 (default and gates-off `-x 0 -z 0 -c 0`).

## Versions / commits

- sk2bgrow: commit `275778f350b87e10c6366bf90884d965fbab45a6` + uncommitted
  working-tree fixes to `python/sk2bgrow/{fit.py,ztp.py}` (the FIXED estimator,
  same state as C1): `$BASE/src/target/release/sk2bgrow`, stats via
  `$BASE/micromamba/envs/sk2bgrow/bin/python -m sk2bgrow.cli`.
- Pilea: v1.3.8, `$BASE/bench/C1/vendor/env-pilea138/bin/pilea` (bioconda
  micromamba env, same as C1).
- seqkit v2.10.0 (`$BASE/tools/seqkit`).

## Species, references, accessions

run -> species from ENA `scientific_name` (filereport, 2026-08-26):
K-* = Klebsiella pneumoniae, B-* = Bacillus subtilis, M-* = Morganella morganii,
P-* = Pseudomonas putida (`run_species.tsv`).

| species | ref file ($BASE/refs) | accession | why |
|---|---|---|---|
| B. subtilis | Bacillus_subtilis_168.fna | GCF_000009045.1 (str. 168) | the lab-strain reference; same file as local benches/genomes (md5 81e24722...) |
| K. pneumoniae | Klebsiella_pneumoniae.fna | GCF_000240185.1 (HS11286) | RefSeq reference; local benches/genomes file (md5 9f1efc72...). 7 contigs: chromosome + 6 plasmids, 5,682,322 bp total |
| M. morganii | Morganella_morganii.fna | GCF_006094455.1 (ATCC 25830) | RefSeq `refseq_category: reference genome`, complete, single contig 3,889,903 bp. Fetched from NCBI 2026-08-26 (md5 144fa28e...) |
| P. putida | Pseudomonas_putida_KT2440.fna | GCF_000007565.2 (KT2440) | the KT2440 reference, complete, single contig 6,181,873 bp. Fetched from NCBI 2026-08-26 (md5 a16461d4...) |

Isolates are wastewater strains, not these type strains — divergence is part
of the experiment (ZTP truncation handles absent anchors; detected-fraction is
reported by QC).

## FASTQ state warning (deviation from the C1b brief)

The brief said fastq/PRJNA1280254 was "complete + verified". At start
(2026-08-26 08:20 HKT) a re-verification/re-download loop was STILL RUNNING on
io2: 12 files in log/failed.log, .part files actively growing, finals being
atomically replaced by ENA-original downloads. 35 of 40 final files did NOT
byte-match the ENA manifest (fastq_bytes); they are ~+4.5% larger, consistent
with SRA-normalized re-encodes (headers `@SRR....N <orig> length=150`, no /1
/2) rather than ENA originals (`@...DP.../1`). Read CONTENT is equivalent
(same reads, same order; md5 mismatch is a gzip-encoding difference, not data).
Decision: proceed with the on-disk files; per-run input read-count
verification against the counts manifest is folded into 01_subsample.sh
(fq/{run}.input_check.tsv). Atomic mv replacement makes mid-run reads
consistent per open file. Any cell whose mate-sync check FAILs is re-run after
the downloads settle. Recorded here so the provenance is explicit.

## Databases

- per species: `sk2bgrow index refs/<sp>.fna -o db/<sp> --enzymes all` (full
  16-enzyme panel) and `pilea index <ref> -o pileadb/<sp> -t 8` (defaults
  k=31 s=250 w=25000). Per-species anchor counts: [TBD after 00_index].

## Subsampling

- Depths 0.5/1/2/5/10 (x species assembly length, plasmids included for K).
- fraction = depth * glen / base_count (base_count from
  manifests/PRJNA1280254_counts.tsv, both mates). Depths with frac >= 1 are
  skipped (none expected: min run ~1.14 Gb vs 10x of the largest assembly
  6.18 Mb = 62 Mb).
- seed = depth_code*1000 + ridx; depth_code = depth*100 (50..1000); ridx =
  1..20 sorted accession order. Same seed both mates (seqkit sample keeps
  record if uniform < p; nested, deterministic — verified in C1).
- Outputs fq/{run}.{depth}x_{1,2}.fq.gz + {run}.{depth}x.info.tsv (seed,
  fraction, mate counts, sync OK/FAIL).

## Pipeline stages

- count: `sk2bgrow profile fq1 fq2 -d db/<sp> -o counts/{name} --no-stats
  --threads 8 --quiet`
- stats: `python -m sk2bgrow.cli profile <counts.tsv> --db db/<sp> --output
  counts/{name} --windows counts/{name}/windows.tsv` (fixed estimator)
- pilea default: `pilea profile fq1 fq2 -d pileadb/<sp> -o
  pilea/default/{name} -t 8`; relaxed: add `-x 0 -z 0 -c 0`. Output file is
  output.tsv.
- /usr/bin/time -v per stage per cell; logs $BASE/logs/C1b/.

## SLURM

- Partition intel, arrays 0-19 (per run, looping 5 depths), cell-level skip
  checks make reruns incremental. QOS: MaxSubmitJobsPerAccount=85 /
  MaxJobsPerAccount=75 (normal) — submitted with %25 throttle and sequential
  stage submission.

## Ground truth (Part 2.2 labelling)

Per-run measured growth rates are NOT published numerically. They exist only
as points in Supplementary Fig S4 (x-axis: mu from OD595, t1=5.75h/t2=6h
interval). We digitized the vector PDF (scripts/digitize_s4.py): axis
calibration from tick marks, points from bezier circles; reproduces the
figure's printed per-panel Pearson r to <=2e-4 (0.3937/0.7870/0.7918/0.7318 vs
0.3939/0.7871/0.7918/0.7318) -> digitization error negligible. Two samples
excluded by the paper (negative measured mu): K at 1/10x, P at 1/3x.
run -> concentration is not in SRA metadata; resolved empirically by matching
our Pilea estimates to the figure's y-values (scripts/resolve_conditions.py,
hypotheses fwd/rev + unconstrained best assignment cross-check).
- correlation: vs digitized mu (and vs nutrient concentration, ordinal
  robustness column).
- RMSE/bias: vs pred_log2ptr_CH = mu*C/ln2 with the PAPER'S OWN
  C ~= 0.6032 + 0.2948*tau — an AFFINE function of mu, not an independent
  marker-frequency target like C1's. Labelled accordingly.
- slope: OLS of estimate on measured mu (table2-style).

## Progress log

- [setup] 4 refs pushed to $BASE/refs, md5 vs local/NCBI verified.
- [setup] Fig S4 digitized + validated (r match <=2e-4). growth_rates_fig.tsv.
- [setup] scripts 00-04 + resolve_conditions.py + analyze_c1b.py staged.
- [TBD] 00_index submitted.
- [index] job 3944837 done in <2 min. sk2bgrow anchors (16-enzyme panel):
  B. subtilis 35,271 (8,366/Mb); K. pneumoniae 53,428 (9,406/Mb, incl. 6
  plasmid contigs); M. morganii 32,489 (8,353/Mb); P. putida 60,048
  (9,717/Mb). E. coli K-12 reference (C1): 43,735 (9,425/Mb). Density spread
  across the four species is only +/-6% around the E. coli value -- a thin
  range; any species effect seen here is not a density effect of this panel.
  Pilea sketch DBs built per species (k=31 s=250 w=25000).
- [subsample] array 3944838 (0-19%25) complete: 100/100 cells, all mate-sync
  OK, 1.7 GB. Input verification: all 20 runs x both mates match the counts
  manifest read counts EXACTLY (fq/*.input_check.tsv all OK), including the
  SRA-normalized re-encoded files -> read content matches the manifest; the
  byte mismatches are gzip-encoding differences only.
- [grid] count array 3944887 submitted; pilea array 3944888 queued
  (afterany:3944887). stats submission hit QOSMaxSubmitJobPerUserLimit
  (contention with fqconv/m2b pipelines); retry pending.
- [grid] count 3944887 + pilea 3944888 + stats 3944937 all COMPLETED: 100/100
  cells per arm (sk2bgrow output.tsv, pilea default, pilea relaxed).
- [analyze 2026-08-26] resolve_conditions.py: species-name bug fixed (fig TSV
  uses Bacillus_subtilis / Pseudomonas_putida, not the ref names; B/P had
  silently matched n=0). Mapping decision: fwd (methods listing order) kept
  for ALL species as the a priori hypothesis; rev/best are cross-check only.
  B/K/M consistent with fwd (SSE 0.021/0.064/0.055 vs rev 0.158/0.169/0.075);
  P unidentified -- 10x Pilea estimates nearly flat across runs (range 0.28),
  fwd 0.420 vs rev 0.374 both terrible; choosing rev would fit noise.
  analyze_c1b.py excluded-block bug fixed (rr.name -> rr["run"]).
  Outputs: c1b_results_raw.tsv, c1b_summary.tsv, c1b_summary_pooled.tsv,
  c1b_cost_raw.tsv; full print in logs/C1b/analyze_final.txt.
