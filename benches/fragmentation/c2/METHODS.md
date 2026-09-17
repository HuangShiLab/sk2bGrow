# C2 — fragmentation, generalised + WCG null-margin validation: METHODS (running log)

Experiment C2 of the HPC phase (`benches/RESEARCH_PLAN.md` Part 6, and Part 4
R2/R3). Two tasks:

- **Task 1 (priority): WCG null-margin validation on real data.** The R3
  within-contig-gradient gate (`benches/fragmentation/qc_within_contig.py`,
  two branches: slope `A_loc` vs global fit, jump `J` at contig boundaries)
  was tuned on one genome / one seed. Its null margins (max J = 0.013 vs the
  0.03 fire line; max A_loc − b_hat = 0.18 on growing controls vs 0.5) are
  re-measured here on real reads at scale: every C1 cell (315, K-12 complete
  reference, Zheng reads, 0.25–20×) and every C1b cell (100, 4 species,
  divergent wastewater isolates vs type-strain references, 0.5–10×) is a
  NULL for WCG (correct coordinate). Scored offline from the existing
  `windows.rates.tsv` + `output.tsv`; no estimator re-runs.
- **Task 2: fragmentation, generalised** (RESEARCH_PLAN C2, in priority
  order): (1) genuinely incomplete MAGs (drop + splice, not just cut);
  (2) scaffolding distance ladder O157:H7 → Shigella → Salmonella → Vibrio
  (R2's failure boundary); (3) the 16-genome simulation grid under
  fragmentation; (4) the enzyme-panel sweep under fragmentation (A6
  blocker).

## Versions / commits

- sk2bgrow: commit `275778f350b87e10c6366bf90884d965fbab45a6` + uncommitted
  working-tree fixes to `python/sk2bgrow/{fit.py,ztp.py}` (the FIXED
  estimator, same state as C1/C1b): `$BASE/src/target/release/sk2bgrow`,
  stats via `$BASE/micromamba/envs/sk2bgrow/bin/python -m sk2bgrow.cli`.
- WCG statistic: `scripts/qc_within_contig.py`, byte-identical copy of
  `benches/fragmentation/qc_within_contig.py` (laptop). Thresholds as tuned
  in R3_QC.md: slope branch fires if usable (≥2 usable contigs, median ≥20
  windows/contig) AND z_loc>2 AND z_short>2 AND A_loc−b_hat>0.5; jump branch
  fires if ≥8 boundary pairs AND z_J>3.5 AND J>0.03.
- Simulator for Task 2.3: `scripts/simulate_c3.py`, byte-identical copy of
  the C3 script (assembly-wide V-profile, split pick/sample seeds).
- Reference genomes: local `benches/genomes/*.fna` pushed to
  `$BASE/refs/genomes16/` (md5-verified against local, e.g. K-12
  92c997bcd88e983ffdb21b2712ed3736; K-12 identical to the C1 reference).

## Deviations / decisions stated up front

- **Cell counts.** The C2 brief said "C1: 945 cells". `bench/C1/counts/`
  holds 315 sk2bgrow cells (45 runs × 7 depths); 945 counts the three arms
  (sk2bgrow + Pilea default + Pilea relaxed). Only sk2bgrow cells carry
  `windows.rates.tsv`, so the WCG null set is 315 + 100 = 415 cells.
- **Single-contig references are vacuous for WCG — by design.** K-12, B.
  subtilis, M. morganii and P. putida are single-record references; both
  branches are undefined at N=1 (A_loc needs ≥2 usable contigs, J needs
  boundary pairs), so those cells contribute "cannot fire", not a measured
  margin. The only real-data multi-contig correct reference in C1/C1b is
  K. pneumoniae (chromosome + 6 plasmids) — see the 2026-08-26 finding
  below. Measured J/A_loc margins on multi-contig correct *fragment*
  references come from the Task-2 controls (ordered/scaffolded drafts).
- **Reads for the real-read grid are the C1 subsamples** (`bench/C1/fq`,
  17 laptop-pick runs: 16 media + RUN_OUT rep SRR11558950, depths
  0.5/1/2/5/10×). No re-subsampling; ground truth
  `bench/C1/growth_rates.tsv` via `runs.tsv`. The complete-reference k=16
  arm is NOT rerun — C1 results for the same cells are the control.
- **Panel subsets (Task 2.4).** The laptop panel-sweep subset definitions
  were not recoverable from the repos. Nested density-ranked subsets
  (usable anchors on K-12) are used and recorded exactly:
  k2 = CjeI,CjePI (this pair reproduces the laptop k=2 index size of
  17,055 anchors exactly); k4 = +HaeIV,Hin4I; k8 = +BcgI,BslFI,AlfI,Bsp24I;
  k16 = all. Applied to both frag100 and the complete K-12 reference.
- **Incomplete-MAG variants (Task 2.1)** from the seed-0 frag100 draft
  (identical to the laptop draft): comp90/comp75/comp50 (drop 10/25/50
  random contigs, rng seed 0, independent draws); comp75_c10_O157H7,
  comp50_c10_O157H7, comp75_c10_Salmonella (splice in 10 foreign contigs
  cut from the named genome with the identical fragment.py protocol,
  seed 0). Composition in `refs/refs_manifest.tsv`.
- **frag16 simulation references (Task 2.3)**: each of the 16 genomes cut
  into 100 contigs (longest record only; plasmids/other records appended
  unchanged), original file stems kept so index genome names match the
  simulation truth. Reads are simulated on the COMPLETE genomes (truth on
  the real coordinate) and profiled against the frag16 DB; complete16 on
  identical reads is the control arm. Laptop seeds: pick = s*1000+c*10+r,
  sample = pick*100+1.

## SLURM

- Partition intel, job names C2_*. Logs `$BASE/logs/C2/`.
- `C2_10_wcgnull.sbatch`: 1 task, 4 CPU, 8 GB — Task 1 scorer.
- `C2_20_build_refs.sh`: 1 task, 16 CPU, 32 GB — all Task-2 references +
  anchor DBs + `refs.tsv`/`runs.tsv` + anchor-count verification.
- `C2_30_profile.sh`: array 0-84 (17 runs × 5 depths) × 17 reference
  conditions from refs.tsv, count + FIXED stats per cell, skip-if-exists.
- `C2_40_simgrid.sh`: array 0-24-task grid s{4,8,16}×c{1,2,4,8}×r{1,2},
  two arms (frag16, complete16) per cell. Ready; submission pending slots.

## Progress log

- [setup 2026-08-26] scripts pushed to `$BASE/bench/C2/scripts/` (fragment,
  rescaffold, qc_within_contig, spread/score_spread, r3_*, analyze/sweep +
  new c2_* drivers); 16 genomes pushed to `$BASE/refs/genomes16/` (md5
  verified).
- [wcgnull 2026-08-26] job 3945186: first pass scored only 305/415 cells —
  cell-name parser bug (`rsplit` broke on `0.25x`/`0.5x` tags, dropping all
  90 C1 + 20 C1b sub-1× cells). Fixed (`split(".",1)`), rerun as 3945217
  (415/415 scored; `c2_wcg_null.tsv`).
- [wcgnull RESULT 2026-08-26] Task-1 numbers (415 null cells):
  - **Fire rate: 5/415 (1.2%), all five on the jump branch, all
    K. pneumoniae** (SRR34095773.5x, SRR34095774.{5,10}x,
    SRR34095775.{1,5}x; J = 0.29–2.17 vs the 0.03 fire line, z_J 3.88–4.62
    vs 3.5). Expected ~0 — the laptop null margin (max J = 0.013 on
    correctly-ordered fragment drafts) does NOT transfer to multi-replicon
    references.
  - Slope branch: 0 fires anywhere; max A_loc−b_hat on growing cells =
    **0.0088** vs the 0.5 fire line (margin ×57); max z_short = 0.18.
    The A_loc side of the gate is comfortable on real data.
  - C1 single-contig cells (315): both branches undefined at N=1 by
    construction → quiet trivially. They are not a measured margin.
  - Mechanism of the K. pneumoniae fires (job 3945232,
    `c2_kpneumo_plasmid_check.tsv`): EVERY boundary pair in all 25 K cells
    touches a plasmid — restricting J to the primary contig leaves 0
    boundary pairs (chromosome is one contig). Plasmid median centred rates
    sit +0.65 to +1.7 log2 above the chromosome (copy-number offset,
    measured on SRR34095774.10x), so the boundary step is real coverage
    biology, not coordinate scrambling. J>0.03 on 25/25 K cells; only
    5/25 also pass z_J>3.5 (the z term is doing the gating). J also
    inflates at 0.5–1× (0.9–2.2) — the A4 low-count SE miscalibration
    reaching the jump branch's noise floor.
  - Consequence for R3: the jump branch as tuned assumes contigs are
    fragments of ONE replicon. On multi-replicon references (chromosome +
    plasmids) it must not run (or must be restricted to primary-replicon
    boundaries). The "correctly-ordered multi-contig draft" null margin
    still has no real-data measurement — it comes from the Task-2
    scaffolded/ordered controls.
  - The R3_QC.md "slope branch only at ≥2×" recommendation is NOT tested
    by C1: at N=1 the slope branch is undefined at every depth, so the
    known 1× A_loc inflation cannot appear here. K cells at 0.5–1× show
    A_loc medians 0.048–0.083 (no inflation, but n=10 and multi-replicon).
    Verdict: recommendation neither confirmed nor refuted by real data yet;
    the frag100/ordered Task-2 cells at 0.5–1× are the test bed.
- [refs 2026-08-26] job 3945218 (C2_refs) COMPLETE. frag100 reproduces the
  laptop draft exactly (N50 77,939; 43,707 anchors vs laptop 43,707 OK).
  complete_k2 = 17,055 anchors — matches the laptop k=2 panel size exactly,
  validating the density-ranked subset choice for that anchor point.
  Incomplete variants: comp90/75/50 = 89.6/73.6/48.2% of bp; contaminated
  variants in `refs/refs_manifest.tsv`.
- [scaffold ladder RESULT 2026-08-26] (Task 2.2, R2 boundary; score files
  `refs/scaf_*.score.txt`): the laptop O157:H7 result reproduces digit-for-
  digit (99/100 placed, orientation 99/99, median start error 712,324 bp,
  73,276 bp after rotation removal, order Spearman 1.0000). Beyond that:
  - Shigella dysenteriae: 81/100 placed (64.9% of bp), orientation 65.4%,
    Spearman 0.9779 — degradation begins at the genus edge.
  - Salmonella enterica LT2: 59/100 placed (77.9% of bp), orientation
    93.2%, Spearman 0.9994 among placed — order survives, coverage of the
    draft does not.
  - Vibrio cholerae (different order): 2/100 placed — scaffolding fails.
    The R2 failure boundary sits between family (Enterobacterales, partial)
    and order (Vibrionales, fail) for this draft.
  Downstream accuracy (r/slope/RMSE per depth, and whether WCG flags the
  degraded scaffolds) comes from the C2_prof grid.
- [wcgnull verified 2026-08-26] `c2_wcg_null.tsv` final: 415 data rows
  (315 C1 + 100 C1b), all depth bands populated, 5 flag=True rows (the K.
  pneumoniae cells listed above). Task 1 deliverable complete.
- [submit 2026-08-26] C2_30_profile.sh (array 0-84) and C2_40_simgrid.sh
  (array 0-23) submissions REJECTED by QOS (MaxSubmitJobsPerAccount /
  QOSMaxSubmitJobPerUserLimit). Measured: SLURM counts ARRAY TASKS
  individually against MaxSubmitJobsPerUser=50 (C3 + C4 arrays hold ~48),
  so the 85-task array can never fit whole. Switched to chunked submission
  with a polling retry loop (600 s), priority order:
    cd $BASE/bench/C2/scripts
    sbatch --array=0-16%17  C2_30_profile.sh   # chunk 0
    sbatch --array=17-33%17 C2_30_profile.sh   # chunk 1
    sbatch --array=34-50%17 C2_30_profile.sh   # chunk 2
    sbatch --array=51-67%17 C2_30_profile.sh   # chunk 3
    sbatch --array=68-84%17 C2_30_profile.sh   # chunk 4
    sbatch --array=0-23%12  C2_40_simgrid.sh   # Task 2.3, last
  Skip-if-exists per cell makes the chunks independent and reruns safe.
  Scoring when done: `$PY scripts/c2_score.py` (real-read grid →
  c2_cells.tsv + c2_summary.tsv) and `$PY scripts/c2_score_sim.py`
  (sim grid → sim/c2_sim_summary.tsv).

## Task 1 — final summary (null margins on real data, job 3945217)

415/415 null cells scored (`c2_wcg_null.tsv`). Per dataset × depth band, the
null-side maxima against the fire lines (J fires > 0.03, z_J > 3.5,
A_loc−b_hat > 0.5). C1 (K-12, single contig) and the three single-contig
C1b species are vacuous: both branches undefined at N=1, quiet by
construction — they contribute "cannot fire", not a measured margin. The
only measured multi-contig null is C1b K. pneumoniae (chromosome + 6
plasmids); maxima below are over the cells where the statistic is defined.

| dataset | band | n | max J | max z_J | max A_loc−b_hat | fires |
|---|---|---|---|---|---|---|
| C1 | 0.25–1× | 135 (45 runs × 3 depths) | undefined (N=1) | — | — | 0 |
| C1 | 2–5× | 90 | undefined | — | — | 0 |
| C1 | 10–20× | 90 | undefined | — | — | 0 |
| C1b | 0.5–1× | 40 (10 K.) | 2.2418 | 3.88 | 0.0088 | 1 |
| C1b | 2–5× | 40 (10 K.) | 0.7024 | 4.62 | 0.0079 | 3 |
| C1b | 10× | 20 (5 K.) | 0.5008 | 3.98 | −0.0997 | 1 |

- **Null fire rate: 5/415 (1.2%)** — expected ~0. All five are the jump
  branch on K. pneumoniae; slope branch fires 0 anywhere. Slope-branch
  margins are comfortable: max A_loc−b_hat on growing cells = 0.0088 vs the
  0.5 line (~57× headroom), max z_short = 0.18 vs 2.
- **Firing cells** (all C1b K. pneumoniae, growing):
  SRR34095773.5x (J=0.386, z_J=4.53; qc: only 11/16 enzymes fit);
  SRR34095774.10x (J=0.309, z_J=3.98; qc clean);
  SRR34095774.5x (J=0.287, z_J=4.37; qc clean);
  SRR34095775.1x (J=2.169, z_J=3.88; qc: fraction 0.40 < 0.75, containment
  0.38 < 0.5, only 6/16 enzymes fit);
  SRR34095775.5x (J=0.407, z_J=4.62; qc clean).
- **Interpretation.** These are null-side fires on a divergent
  multi-REPlicon reference, and the driver is measured, not guessed
  (job 3945232, `c2_kpneumo_plasmid_check.tsv`): every boundary pair in all
  25 K. pneumoniae cells touches a plasmid — restricting J to the primary
  contig leaves zero boundary pairs on 25/25 cells. The isolate's plasmids
  sit at a different copy number than its chromosome (median centred rates
  +0.65 to +1.7 log2 above the chromosome on SRR34095774.10x), so each
  chromosome→plasmid boundary carries a genuine coverage step of the size J
  reports. That is a real step, not a scrambled coordinate: the jump
  branch's design assumption — contigs are fragments of ONE replicon — is
  violated by a reference whose contigs are separate replicons. Two
  aggravators: the isolate is highly diverged from the type strain
  (detected fraction 0.40 at 1×, so plasmid windows are sparse and noisy),
  and the J noise floor inflates at ≤1× (J = 0.9–2.2 across all K cells at
  0.5–1× — the A4 window-SE miscalibration reaching the jump branch; the
  z_J > 3.5 term, not the J > 0.03 term, is what kept the fire count at 5).
  Fix direction: do not run the jump branch on multi-replicon references
  (or restrict it to primary-replicon boundary pairs); the tuned thresholds
  themselves need no change for the fragmented-draft class they were
  designed for, but that class still lacks a real-data null measurement —
  it comes from the C2 scaffolded/ordered controls.
- **Depth-band recommendation.** R3_QC.md's "slope branch only at ≥2×" is
  NOT tested by C1: at N=1 the slope branch is undefined at every depth, so
  the known 1× A_loc inflation cannot appear. K cells at 0.5–1× show no
  A_loc inflation (medians 0.048–0.083) but n=10 and multi-replicon. The
  frag100/ordered Task-2 cells at 0.5–1× are the real test bed. The jump
  branch's ≤1× J inflation (above) argues for extending the same ≥2×
  restriction to the jump branch.

## Submission state (2026-08-26, final)

- SLURM counts array TASKS individually against MaxSubmitJobsPerUser=50;
  C3+C4 arrays held ~48 slots all day, so the grids are submitted in
  8-task chunks by a polling submitter (C5 02_submit_grid.sh pattern):
  **job 3945810 (`C2_submit`, 48 h limit, RUNNING)** — submits
  C2_30_profile.sh chunks 9-16 … 72-79 in order at ≤40 queued tasks, then
  C2_40_simgrid.sh 0-23%12 at ≤35, with a task-0 safety net. Markers in
  `bench/C2/.submit/`.
- Already queued before the submitter: pilot task 0 (3945569), chunk 1-8
  (3945712), chunk 80-84 (3945655) — all PENDING (Priority) behind C3/C4.
- Scoring once the grids finish:
    $PY $BASE/bench/C2/scripts/c2_score.py       # real-read grid →
                                                 # c2_cells.tsv, c2_summary.tsv
    $PY $BASE/bench/C2/scripts/c2_score_sim.py   # sim grid →
                                                 # sim/c2_sim_summary.tsv
