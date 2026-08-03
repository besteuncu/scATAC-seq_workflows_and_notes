# scATAC-seq workflows and notes

Curated notes, reproducible workflows, tool comparisons, and reading lists for single-cell ATAC-seq analysis.

**Author:** Beste Uncu, PhD in Statistics (2026)

**Focus:** Reproducible end-to-end scATAC-seq pipelines in both R and Python, alongside topic notes that view current-literature methods through a statistical lens.

## Why this repo

Single-cell RNA sequencing and spatial transcriptomics have opened remarkable windows into cellular identity and tissue architecture. But cellular state is not expression alone; it is also chromatin accessibility, transcription factor binding, and the wider regulatory landscape. Reading these layers together, particularly in a disease context, surfaces mechanisms that expression data alone cannot resolve. scATAC-seq brings the regulatory potential of each cell into focus.

This repo is a working reference: notes and reproducible pipelines built by a statistician learning single-cell chromatin accessibility from the ground up. The framing is intentional. Statistical clarity on why TF-IDF, why LSI, why the first component is dropped, why pseudobulk beats cell-level testing. Concepts read through that lens first, code follows.

Inside are four things: topic notes covering biology, statistics, and pipeline conventions; end-to-end reproducible workflows in R (Signac) and Python (SnapATAC2) on the same public dataset; a tool comparison across Signac, ArchR, and SnapATAC2; and a curated reading list. The workflows are the anchor. The notes explain the choices they encode.

## Repository layout

```
scATAC-seq_workflows_and_notes/
├── docs/            Topic notes (biology, stats, QC, dim red, ChromVAR, etc.)
├── workflows/       Reproducible pipelines (Signac PBMC 10k, SnapATAC2 parallel)
├── tools/           Tool comparison (Signac vs ArchR vs SnapATAC2)
├── papers/          Curated reading list with one-line summaries
└── benchmarks/      Multi-tool comparisons on the same data (roadmap)
```

## Table of contents

Sections are in progress; this table will populate as content lands.

### Concepts and statistics

1. Biology and assay landscape (planned)
2. Data structure and statistical model (planned)
3. Quality control (planned)
4. Normalization and dimensionality reduction (planned)
5. Clustering (planned)
6. Cell type annotation (planned)
7. Peak calling (planned)
8. TF motif activity (ChromVAR) (planned)
9. Peak-to-gene linking (planned)
10. TF footprinting (planned)
11. Batch integration (planned)
12. GWAS variant overlay (planned)
13. Common pitfalls (planned)

### Reproducible workflows

- Signac PBMC 10k pipeline (R) (planned)
- SnapATAC2 PBMC 10k pipeline (Python) (planned)

### Reference material

- Tool comparison (planned)
- Papers with one-line summaries (planned)
- Benchmarks (roadmap)

## Running the workflow

Workflows land under `workflows/`. Each ships with its own README covering data download, environment setup, and script execution order.

Preview of the R pipeline (Signac):

```bash
cd workflows/signac_pbmc_10k
# See workflow README for details
Rscript 01_qc.R
Rscript 02_normalize_dimred.R
Rscript 03_cluster_annotate.R
Rscript 04_chromvar.R
Rscript 05_linkpeaks.R
```

Preview of the Python pipeline (SnapATAC2):

```bash
cd workflows/snapatac2_pbmc_10k
# See workflow README for details
conda env create -f environment.yml
conda activate scatac-py
python 01_qc.py
# ...
```

Expected runtime on a modern laptop: about 25 minutes for 10k PBMC on 16 GB RAM.

## Environment

Primary stack:

- **R 4.5+** with Signac 1.14+, Seurat 5.x, GenomicRanges, EnsDb.Hsapiens.v86, BSgenome.Hsapiens.UCSC.hg38, chromVAR, motifmatchr, JASPAR2020.
- **Python 3.11** (isolated conda env) with SnapATAC2, scanpy, anndata, scvi-tools.
- Package versions locked per workflow via `renv.lock` (R) or `environment.yml` (Python).

## References that shaped this repo

- Stuart et al. 2021. Single-cell chromatin state analysis with Signac. *Nat Methods*.
- Granja et al. 2021. ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis. *Nat Genet*.
- Schep et al. 2017. chromVAR: inferring transcription-factor-associated accessibility from single-cell epigenomic data. *Nat Methods*.
- Zhang et al. 2024. SnapATAC2: a fast, scalable and versatile tool for analysis of single-cell omics data. *Nat Methods*.
- Cusanovich et al. 2018. A Single-Cell Atlas of In Vivo Mammalian Chromatin Accessibility. *Cell*.
- Buenrostro et al. 2015. Single-cell chromatin accessibility reveals principles of regulatory variation. *Nature*.

Full curated list with one-line summaries will live in [papers/](papers/).

## Similar community reference repos

Structure and coverage inspired by:

- [crazyhottommy/scATACseq-analysis-notes](https://github.com/crazyhottommy/scATACseq-analysis-notes)
- [mdozmorov/scATAC-seq_notes](https://github.com/mdozmorov/scATAC-seq_notes)
- [Liu Lab, bioinformatics-in-combio, scATAC chapter](https://liulab-dfci.github.io/bioinfo-combio/scatac.html)

## License

MIT. See [LICENSE](LICENSE).

## Contact

Beste Uncu, besteuncu@gmail.com, [github.com/besteuncu](https://github.com/besteuncu)
