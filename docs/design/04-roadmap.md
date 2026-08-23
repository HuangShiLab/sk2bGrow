# Roadmap

Two interleaved tracks: **M1–M4** from the architecture document (what to build)
and **P0–P4** from the algorithm report (what to prove). Building without the
matching proof is how a method ships that nobody can trust.

## M1 — minimal closed loop ✅

`digest` + `anchor_db` + `count` + `ztp` + `fit`, running end to end and
recovering planted PTR values in simulation.

Delivered, plus more than planned: the V-shape coordinate fit, cross-enzyme
fusion and the ori search all exist, because implementing the parity path alone
would have left the central claim untested.

Evidence in the repository:

* `cargo test --workspace` — unit + integration, including a synthetic genome
  whose planted replication gradient is recovered from simulated reads.
* `pytest` — the statistical core, including PTR recovery at 0.5/1.0/2.0 through
  the real CLI.
* `python -m sk2bgrow.cli simulate a|b` — reproduces the report's §5 ordering.

## P0 — in silico benchmark

Replicate Pilea's 400-sample simulated-community design (4–32 strains × 4–32×,
log2 PTR ~ U[0,2]) with the same estimator over three anchor sets: random sketch,
single enzyme, 16-enzyme union.

*Gate:* union L2 distance < random sketch, with the gain quantified in the 1–2×
band. `simulate.route_a` is the harness; feed it real digested coordinates via
`anchors_from_digest` rather than the synthetic sets.

## M2 / P1 — A/B against Pilea on real data

Zheng E. coli, 16 growth conditions (PRJNA615952, >300×, Pilea reaches
r = 0.9764). Subsample to 0.5–10× and run both tools through one script.

*Gate:* correlation with measured growth rate ≥ Pilea at every coverage, and a
clear win at 1×.

**This is the honest gate.** If sk2bGrow cannot beat the Pilea baseline at 1×
subsampling, the premise needs rethinking rather than more engineering — the
whole argument for a deterministic sketch is the low-coverage band. Protocol in
[`../../benches/README.md`](../../benches/README.md).

## M3 / P2 — the innovation layer, and the wet-lab calibration

Code: ori-aware V-shape MLE, cross-enzyme fusion + QC, `--mode 2brad`. All three
exist; M3 is about validating them rather than writing them.

Experiment: the 4-species LB gradient (PRJNA1280254), one BcgI library plus a 3–4
enzyme combination per sample, three technical replicates.

*Gate:* per-enzyme PTR consistency; a **measured** σ_eff; correlation with a
simultaneous WMS-PTR gold standard.

σ_eff is the single largest open risk (report R1). Every route-B claim rests on
assumed values of 0.3–0.6, and no published per-site depth CV exists for 2bRAD.
If σ_eff turns out above ~0.8, the route-B gain estimates need revising downward.

## M4 / P3–P4 — scale and community data

* Two-tier architecture: containment pre-filter, then anchor-level counting, to
  keep GTDB-scale profiling tractable (~190 GB raw for 137 k representatives).
* `scaffold` integrated into the profile path for fragmented MAGs.
* `dynamics` on the Long marine time series (PRJNA551656; 98 MAGs with measured
  cell-count growth rates).
* Mock communities (Mock-CAS / MSA1002) for abundance–PTR decoupling.

*Gate:* end-to-end time on a 100 Gbp dataset comparable to Pilea's 106 s; a
single-sample four-dimensional profile (taxonomy + ANI + SV + PTR) demonstrated.

## Open questions

| | question | resolution path |
|---|---|---|
| R1 | actual magnitude of per-site efficiency noise | P2 technical replicates |
| R2 | shared anchors across co-occurring strains | strain-specific anchors only; accept the sensitivity loss |
| R3 | multi-fork profile non-linearity | segmented model exists; needs validation at PTR > 2 |
| R4 | ori annotation coverage for MAGs | joint search + `ori_confidence`, cross-checked against Ori-Finder |
| — | `BslFI` density vs its `GGGAC` motif | reconcile against `bsyn::enzyme` ([`../enzymes.md`](../enzymes.md)) |
| — | HaeIV ⊂ Hin4I breaks strata independence | handled in counting; quantify the effect on the χ² test |
