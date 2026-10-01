# Benchmark evaluation and scientific-significance assessment

This assessment uses the five-pillar benchmark framework and applies it to the
current sk2bGrow implementation and the 2026-09-30/2026-10-01 refresh.

## Overall verdict

**Accept with revisions.** The work is scientifically useful because it turns a
motif-defined restriction digest into a coordinate-aware, stratified,
reproducible PTR benchmark and companion estimator. The strongest evidence is
the controlled Zheng comparison and the low-depth simulation. The main
submission risk is scope: the real-data validation is still one isolate, and
the three sequencing seeds are not biological replicates.

## Five-pillar assessment

| pillar | score /5 | assessment |
|---|---:|---|
| Research gap | 4.5 | Existing sketch-based PTR estimators discard motif coordinates and use arbitrary windows. The need for deterministic loci, coordinate fitting and stratified QC is concrete. |
| Construction pipeline | 4.0 | Enzyme definitions, digest, TGT/anchor database, counting and statistics are implemented and versioned. EM and ZTP boundaries are now explicit. |
| Evaluation framework | 4.0 | The refresh separates all-finite from default-QC estimates, reports media bootstrap CIs, seed instability, negative controls and paired deltas. |
| Empirical findings | 4.0 | Current data show low-depth estimator signal, QC immaturity, sorted-regression failure, and NB-specific residual compression. |
| Companion method | 4.0 | sk2bGrow is both the benchmark consumer and a complete estimator; however, EM is not used for current PTR and WCG remains prototype-level. |

## Scientific significance

### High-value contributions

1. **Coordinate preservation.** Deterministic enzyme anchors retain motif and
   genomic position. This enables a direct V-fit and makes worst-case anchor
   gaps auditable before profiling.
2. **Stratified measurement.** Sixteen enzyme strata give repeated, enzyme-
   specific views of the same replication gradient. They also expose
   between-enzyme disagreement that a single random sketch can hide.
3. **Failure-mode diagnosis.** The benchmark distinguishes estimator signal
   from shipped-QC deployability, sorted-regression failure, stationary floor
   and reference fragmentation.
4. **Reproducibility discipline.** Archive hashes, fixed seeds, paired-delta
   diagnostics and both all-finite/default-QC views reduce selective
   interpretation.

### Boundaries on significance

1. The real-data validation is one *E. coli* isolate, not a community.
2. The three subsamples quantify technical read-sampling variability; they are
   not biological replicates.
3. Pilea gates-off remains competitive under the same all-finite rule, so the
   work should not claim a uniform accuracy win.
4. Low-depth sk2bGrow estimates are present but fail shipped default QC at
   0.5–1×.
5. R3 is prototype evidence and must not be described as production QC.

## Fatal-flaw audit

| candidate flaw | verdict | rationale |
|---|---|---|
| No real-data validation | not fatal | Zheng 2020 provides an independent measured growth-rate target and 1,275 current profiles. |
| No community validation | major, not fatal | Mixed-strain simulation is included, but no real community benchmark is complete. |
| Low-coverage claim overreach | resolved | Results now distinguish all-finite signal from default-QC deployment. |
| EM overclaim | resolved | PTR uses raw integer counts and excludes shared anchors; EM is a diagnostic/assignment interface. |

## Quality-improvement actions completed

1. Added formal manuscript draft, figures and statistical/availability
   supplement.
2. Added media bootstrap CIs and paired A-minus-E deltas.
3. Added seed-level QC pass counts and all-seed/any-seed deployment counts.
4. Added A4 replicate bootstrap CIs for recovered slope and b=0 floor.
5. Added R3 threshold sensitivity; the post-hoc `z_short > 2.5` candidate
   removes the observed stationary false fire but requires independent
   validation.
6. Explicitly scoped EM output away from current PTR estimation.
7. Replaced historical/current number mixing with a provenance gate.
