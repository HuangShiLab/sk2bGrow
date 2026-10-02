# Pre-submission quality review

This review assessed the revised manuscript after integrating the C1b,
C4 and C5 real-community results and generating the Microbiome/BMC package.

## Summary

- CRITICAL: 0
- MAJOR: 0
- MINOR: 3 administrative items
- Numerical audit: 136 claims PASS
- Package: `docs/paper/microbiome_submission/`

## Claim calibration

The manuscript now separates three result classes:

1. Controlled estimator signal: the Zheng all-finite, A4 overdispersion and
   mixed-strain results.
2. Deployment readiness: default-QC pass rates and the prototype status of
   within-contig gradient QC.
3. Real-data limits: C1b cross-species weakness, C4 marine protocol-coverage
   weakness and C5 RBC recall/runtime tradeoff.

No blanket superiority claim over Pilea remains. The C1b statement explicitly
reports that Pilea gates-off is stronger at 1–10×. C4 explicitly reports the
small matched-MAG set and noisy same-reads truth. C5 explicitly states that it
demonstrates recall and workflow coverage, not independent growth-rate
accuracy.

## Venue format

The manuscript uses Microbiome/BMC Research Article structure: title,
unstructured abstract, keywords, Background, Results, Discussion, Conclusions,
Limitations, Methods, Availability of data and materials, Abbreviations,
Declarations and numbered References.

The package contains the main DOCX, title-page DOCX, two additional-file DOCX
files, six figures in PDF/PNG, machine-readable real-community TSVs and BibTeX
and RIS reference exports.

## Checks completed

| check | result |
|---|---|
| real-community claim checker | 136 claims PASS |
| Figure 6 generation | PASS; PNG/PDF generated |
| DOCX ZIP integrity | PASS for all four DOCX files |
| main DOCX image embedding | six figures embedded |
| heading structure | BMC sections present |
| em-dash sentence connector scan | none found |

## Remaining administrative items

1. Enter approved authors, affiliations and corresponding-author details.
2. Complete funding, acknowledgements and CRediT author contributions.
3. Confirm preprint policy and, if applicable, add the preprint DOI.

## Recommendation

**Ready for internal submission review after the three administrative items.**
The scientific text, real-data boundary claims, figures, supplementary tables
and journal package are now integrated.
