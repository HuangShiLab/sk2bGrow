#!/usr/bin/env python3
"""Build a Microbiome/BMC submission package from the maintained Markdown."""
from __future__ import annotations

import datetime as dt
import hashlib
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "docs/paper"
OUT = PAPER / "microbiome_submission"
FIGURES = PAPER / "figures"
REAL = ROOT / "benches/realcommunity_20261001/remote_summaries"

TITLE_PAGE = """---
title: "Microbiome submission title page"
---

# Deterministic 2bRAD anchors enable coordinate-aware bacterial replication-rate estimation

## Authors and affiliations

To be completed in the Microbiome submission system. Do not add authors who have
not approved the manuscript and this submission.

## Corresponding author

To be completed with name, department, institution, postal address, telephone
number, email address and ORCID.

## Running title

Deterministic 2bRAD anchors for PTR benchmarking

## Keywords

bacterial growth rate; peak-to-trough ratio; replication rate; 2bRAD;
metagenomics; sketching; benchmarking; microbiome
"""

REFERENCES_BIB = r"""@article{zheng2020general,
  author = {Zheng, Hai and Bai, Yang and Jiang, Meiling and Tokuyasu, Taku A. and Huang, Xiongliang and Zhong, Fajun and Wu, Yuqian and Fu, Xiongfei and Kleckner, Nancy and Hwa, Terence and Liu, Chenli},
  title = {General quantitative relations linking cell growth and the cell cycle in {\em Escherichia coli}},
  journal = {Nature Microbiology},
  year = {2020},
  volume = {5},
  number = {8},
  pages = {995--1001},
  doi = {10.1038/s41564-020-0717-x}
}

@article{chen2026pilea,
  author = {Chen, Xi and Xu, Xiaoqing and Lin, Yunqi and Shi, Xianghui and Wang, Dou and Zhang, Tong},
  title = {Pilea: profiling bacterial growth dynamics from metagenomes with sketching},
  journal = {Microbiome},
  year = {2026},
  volume = {14},
  number = {1},
  pages = {128},
  doi = {10.1186/s40168-026-02374-0}
}

@article{brown2016measurement,
  author = {Brown, Christopher T. and Olm, Matthew R. and Thomas, Brian C. and Banfield, Jillian F.},
  title = {Measurement of bacterial replication rates in microbial communities},
  journal = {Nature Biotechnology},
  year = {2016},
  volume = {34},
  number = {12},
  pages = {1256--1263},
  doi = {10.1038/nbt.3704}
}

@article{emiola2018high,
  author = {Emiola, Akintunde and Oh, Julia},
  title = {High throughput in situ metagenomic measurement of bacterial replication at ultra-low sequencing coverage},
  journal = {Nature Communications},
  year = {2018},
  volume = {9},
  pages = {4956},
  doi = {10.1038/s41467-018-07240-8}
}

@article{joseph2022accurate,
  author = {Joseph, Tal A. and Chlenski, Phillip and Litman, Aaron and Korem, Tal and Pe'er, Itsik},
  title = {Accurate and robust inference of microbial growth dynamics from metagenomic sequencing reveals personalized growth rates},
  journal = {Genome Research},
  year = {2022},
  volume = {32},
  number = {3},
  pages = {558--568},
  doi = {10.1101/gr.275533.121}
}

@article{long2021benchmarking,
  author = {Long, Andrew M. and Hou, Shengwei and Ignacio-Espinoza, J. Cesar and Fuhrman, Jed A.},
  title = {Benchmarking microbial growth rate predictions from metagenomes},
  journal = {The ISME Journal},
  year = {2021},
  volume = {15},
  number = {1},
  pages = {183--195},
  doi = {10.1038/s41396-020-00773-1}
}

@article{sun2022species,
  author = {Sun, Zheng and Huang, Shi and Zhu, Pengfei and Tzehau, Lam and Zhao, Helen and Lv, Jia and Zhang, Rongchao and Zhou, Lisha and Niu, Qianya and Wang, Xiuping and Zhang, Meng and Jing, Gongchao and Bao, Zhenmin and Liu, Jiquan and Wang, Shi and Xu, Jian},
  title = {Species-resolved sequencing of low-biomass or degraded microbiomes using 2bRAD-M},
  journal = {Genome Biology},
  year = {2022},
  volume = {23},
  number = {1},
  pages = {36},
  doi = {10.1186/s13059-021-02576-9}
}
"""

REFERENCES_RIS = """TY  - JOUR
AU  - Zheng H
AU  - Bai Y
AU  - Jiang M
TI  - General quantitative relations linking cell growth and the cell cycle in Escherichia coli
JO  - Nature Microbiology
VL  - 5
IS  - 8
SP  - 995
EP  - 1001
PY  - 2020
DO  - 10.1038/s41564-020-0717-x
ER  -

TY  - JOUR
AU  - Chen X
AU  - Xu X
AU  - Lin Y
AU  - Shi X
AU  - Wang D
AU  - Zhang T
TI  - Pilea: profiling bacterial growth dynamics from metagenomes with sketching
JO  - Microbiome
VL  - 14
IS  - 1
SP  - 128
PY  - 2026
DO  - 10.1186/s40168-026-02374-0
ER  -

TY  - JOUR
AU  - Brown CT
AU  - Olm MR
AU  - Thomas BC
AU  - Banfield JF
TI  - Measurement of bacterial replication rates in microbial communities
JO  - Nature Biotechnology
VL  - 34
IS  - 12
SP  - 1256
EP  - 1263
PY  - 2016
DO  - 10.1038/nbt.3704
ER  -

TY  - JOUR
AU  - Emiola A
AU  - Oh J
TI  - High throughput in situ metagenomic measurement of bacterial replication at ultra-low sequencing coverage
JO  - Nature Communications
VL  - 9
SP  - 4956
PY  - 2018
DO  - 10.1038/s41467-018-07240-8
ER  -

TY  - JOUR
AU  - Joseph TA
AU  - Chlenski P
AU  - Litman A
AU  - Korem T
AU  - Pe'er I
TI  - Accurate and robust inference of microbial growth dynamics from metagenomic sequencing reveals personalized growth rates
JO  - Genome Research
VL  - 32
IS  - 3
SP  - 558
EP  - 568
PY  - 2022
DO  - 10.1101/gr.275533.121
ER  -

TY  - JOUR
AU  - Long AM
AU  - Hou Z
AU  - Ignacio-Espinoza JC
AU  - Fuhrman JA
TI  - Benchmarking microbial growth rate predictions from metagenomes
JO  - The ISME Journal
VL  - 15
IS  - 1
SP  - 183
EP  - 195
PY  - 2021
DO  - 10.1038/s41396-020-00773-1
ER  -

TY  - JOUR
AU  - Sun Z
AU  - Huang S
AU  - Zhu P
TI  - Species-resolved sequencing of low-biomass or degraded microbiomes using 2bRAD-M
JO  - Genome Biology
VL  - 23
IS  - 1
SP  - 36
PY  - 2022
DO  - 10.1186/s13059-021-02576-9
ER  -
"""


def pandoc(source: Path, target: Path) -> None:
    subprocess.run(
        [
            "pandoc", str(source), "--from", "markdown", "--to", "docx",
            "--resource-path", str(PAPER), "--output", str(target),
        ],
        cwd=ROOT,
        check=True,
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if OUT.exists():
        shutil.rmtree(OUT)
    (OUT / "figures").mkdir(parents=True)
    (OUT / "supplementary_data").mkdir(parents=True)

    title_path = OUT / "title_page.md"
    title_path.write_text(TITLE_PAGE, encoding="utf-8")
    pandoc(title_path, OUT / "title_page.docx")

    pandoc(PAPER / "manuscript.md", OUT / "sk2bgrow_microbiome_research_article.docx")
    pandoc(PAPER / "supplement.md", OUT / "Additional_file_1_supplement.docx")
    pandoc(PAPER / "supplementary_tables.md", OUT / "Additional_file_2_supplementary_tables.docx")

    for stem in [
        "fig1_pipeline", "fig2_zheng_performance", "fig3_a4_compression",
        "fig4_r3_qc", "fig5_mixed_strain", "fig6_real_community",
    ]:
        for suffix in ("pdf", "png"):
            shutil.copy2(FIGURES / f"{stem}.{suffix}", OUT / "figures" / f"{stem}.{suffix}")

    for path in sorted(REAL.glob("*.tsv")):
        shutil.copy2(path, OUT / "supplementary_data" / path.name)
    shutil.copy2(REAL / "METHODS.md", OUT / "supplementary_data" / "realcommunity_METHODS.md")

    (OUT / "references.bib").write_text(REFERENCES_BIB, encoding="utf-8")
    (OUT / "references.ris").write_text(REFERENCES_RIS, encoding="utf-8")

    generated = dt.datetime.now(dt.timezone.utc).date().isoformat()
    manifest_lines = ["# Microbiome submission package", "", f"Generated: {generated}", "", "## Files", ""]
    for path in sorted(p for p in OUT.rglob("*") if p.is_file() and p.name != "README.md"):
        rel = path.relative_to(OUT).as_posix()
        manifest_lines.append(f"- `{rel}` — sha256 `{sha256(path)}`")
    manifest_lines += [
        "",
        "## Before submission",
        "",
        "1. Complete approved authors, affiliations, corresponding author and CRediT contributions.",
        "2. Complete funding and acknowledgements.",
        "3. Confirm whether a prior preprint DOI must be cited.",
        "4. Upload figures as PDF or TIFF and the additional files listed above.",
        "5. Upload title-page information in the journal portal; do not add unapproved authors.",
        "",
        "## Provenance caution",
        "",
        "Zheng, A4, R3 and mixed-strain analyses use the current implementation refresh.",
        "C1b, C4 and C5 real-community runs used an earlier sk2bGrow commit plus recorded",
        "working-tree fixes, as stated in Methods and Additional file 1.",
        "",
    ]
    (OUT / "README.md").write_text("\n".join(manifest_lines), encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
