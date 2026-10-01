# Pre-submission quality review

This review used the benchmark-paper structure, pre-submission severity
taxonomy and figure-quality rules from the installed Supervisor-Skills. It
assesses `docs/paper/manuscript.md` after the mixed-strain rerun.

## Summary

- CRITICAL: 0
- MAJOR: 3
- MINOR: 4
- Top fixes: complete tool citations; add real-community or biological-replicate
  validation; move the manuscript into the target venue template.

## Dimension 1: Macro logic

| # | finding | severity | status / fix |
|---|---|---|---|
| 1 | The manuscript states that “three sequencing subsamples ... are technical read-sampling replicates, not biological replicates.” | MAJOR | Correctly scoped, but a reviewer will still request biological replicates or a real community. |
| 2 | The mixed-strain track is simulated. The text says, “This controlled test measures attribution in a known community. It does not replace real metagenomic validation.” | MAJOR | Keep as a limitation; real-community validation should precede a broad method-paper claim. |
| 3 | The comparison focuses on Pilea. The limitations state that broader comparison with iRep, GRiD, CoPTR and DEMIC is required. | MAJOR | Acceptable for a method-development paper, but not enough for a broad PTR-tool benchmark claim. |
| 4 | The Introduction follows background, three limitations, research questions, design considerations, proposal and contributions. | resolved | No further change needed. |
| 5 | Results separate all-finite and default-QC views throughout. | resolved | Maintain this wording in the venue template. |

## Dimension 2: Writing and claims

| # | finding | severity | status / fix |
|---|---|---|---|
| 1 | The phrase “does not support a uniform accuracy-superiority claim” is negative but evidence-calibrated. | resolved | Keep; do not replace with a positive superiority claim. |
| 2 | The text explicitly says, “current ZTP PTR model does not use EM fractional weights.” | resolved | Keep this scope statement in the main text. |
| 3 | WCG is described as prototype and the `z_short > 2.5` result as post-hoc. | resolved | Keep; do not install WCG as production QC yet. |
| 4 | Some paragraphs in Methods exceed six lines. | MINOR | Split during venue formatting if the journal requests short paragraphs. |

## Dimension 3: English grammar

No grammar error was identified in the sampled full text. Verb tense is
consistent: present tense for findings, past tense for performed experiments.

## Dimension 4: LaTeX and venue format

| # | finding | severity | fix |
|---|---|---|---|
| 1 | The current manuscript is Markdown, not a journal template. | MAJOR | Convert to the target venue LaTeX/Word template after tool citations are completed. |
| 2 | Reference entries 2 and 4 require final version-specific bibliographic details. | MAJOR | Verify and export tool citations before submission. |
| 3 | Supplementary tables are Markdown rather than journal supplementary files. | MINOR | Export the same tables to XLSX or PDF when preparing submission files. |

## Dimension 5: Figures

| # | finding | severity | status / fix |
|---|---|---|---|
| 1 | Figures are vector PDF and raster PNG. | resolved | Use PDF for submission. |
| 2 | Figure captions state the finding first. | resolved | No further change needed. |
| 3 | Figure 1 is a simple pipeline diagram. | MINOR | Replace with a designer-drawn schematic if submitting to a high-visibility venue. |
| 4 | Figure 4 recall is intentionally a constant line at 1.0. | resolved | This is informative because recall, not error, is the quantity. |

## Banned-vocabulary and em-dash scan

The full manuscript was scanned. No prohibited em-dash sentence connector and
no banned intensifier was found. The word root “superiority” occurs once in a
negative calibration statement and is appropriate.

## Numerical consistency

`benches/refresh_20260930/validate_manuscript_numbers.py` compared 97 headline
values and QC counts against the archived TSVs. All checks passed.

## Final score

**7 / 10** for a method-development manuscript at pre-submission stage.

## Recommendation

**Needs major revision before submission.** The computational claims, figures
and statistical reporting are now much stronger, but the manuscript still needs
(1) final tool citations and venue formatting, (2) either real-community/
biological-replicate validation or a narrower title and claim set, and
(3) submission-quality supplementary files.
