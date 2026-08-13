# 01. Multiome introduction — why paired scRNA + scATAC matters

## What Multiome is

The 10x Genomics Multiome ATAC + Gene Expression kit produces **two libraries from the same nucleus**:
1. **scRNA-seq library** — poly-A selection to capture mRNA
2. **scATAC-seq library** — Tn5 tagmentation to capture open chromatin

Both libraries share the **same barcode set**. Every barcode (for example `AAACGAAAGAAACGAA-1`) carries both RNA counts and ATAC peak counts.

Cellranger-ARC (ARC = Assay for RNA + Chromatin) processes both jointly and outputs:
- `filtered_feature_bc_matrix.h5` — RNA and ATAC counts in one file (two assays)
- `atac_fragments.tsv.gz` — ATAC fragments
- `per_barcode_metrics.csv` — per-cell QC (RNA and ATAC metrics)

## Multiome versus unpaired

The classic approach is to run scRNA and scATAC on separate samples and integrate computationally.

| Dimension | Unpaired (separate scRNA + scATAC) | Multiome (paired) |
|-----------|-----------------------------------|-------------------|
| Cell-cell matching | Inferred via anchors (potentially erroneous) | Direct (same barcode) |
| Peak-to-gene inference | Population-level correlation proxy | Direct calculation on paired values |
| Cell type annotation | ATAC uses gene activity proxy | RNA labels directly |
| Priming state | Hard to detect (different cells) | Easy (both modalities in same cell) |
| Cost per cell | Two libraries (~$0.50 × 2) | One prep (~$0.70) |
| Depth per modality | Full for each | Slightly split |
| Reference framework | Bridge integration required | WNN built in |

**Trade-off summary:**
- Multiome wins for regulatory inference, priming, and annotation.
- Unpaired wins when maximum per-modality depth is needed and sample optimization is separately possible.

## Why Multiome matters for regulatory inference

scRNA-seq tells you: "In cluster X, CD8A mRNA is high."
scATAC-seq tells you: "In cluster X, the CD8A promoter is open."

When both measurements come from the same cell, you can state: "In cluster X, CD8A is both accessible and expressed — active transcription."

The four possible states:
- RNA high + ATAC closed → rare; measurement error or declining state
- RNA absent + ATAC open → **priming** — chromatin ready but transcription has not started. Developmental or differentiation state indicator.
- RNA absent + ATAC closed → normal silent gene
- RNA high + ATAC open → actively transcribed

**Priming** is the key concept. It is what unpaired workflows cannot detect. Multiome exposes it directly in a single cell, on a per-cell basis. This is why differentiation, regeneration, and cancer plasticity studies increasingly favor Multiome.

## Peak-to-gene linkage — Multiome's main strength

In scATAC alone, linking an enhancer peak to the gene it regulates is inference:
- The distal enhancer is within 500 kb of the gene's promoter
- Both are open in the same cell type
- Conclusion: probable regulatory link (correlational)

In Multiome, the inference becomes "causal-adjacent":
- Peak accessibility and gene expression are paired at cell resolution
- Strong per-cell correlation → strong regulatory link
- Example: a cluster-specific open peak plus cluster-specific gene expression are directly connected

Functions:
- Signac: `LinkPeaks(pbmc, peak.assay = "ATAC", expression.assay = "SCT")`
- muon: Pearson correlation on paired observations

Peak-to-gene links form the foundation of regulatory network inference. This is Multiome's most substantial contribution over bulk ATAC.

## Cell type annotation — the easy path

In scATAC-only pipelines annotation is difficult:
- Gene activity is a peak-based proxy for expression
- Long genes are over-estimated
- Marker gene promoter accessibility does not always track expression
- Reference-based label transfer (Azimuth) requires reference-query domain match

In Multiome, annotation is trivial:
- RNA modality provides canonical markers (CD3E, MS4A1, CD14) directly
- Full compatibility with the Seurat ecosystem (SCT, PCA, standard workflows)
- WNN clusters reflect groups both modalities agree on

Result: annotation quality is higher than unpaired analysis, with less effort.

## Multiome in this repository

The repository contains two pipeline sets:
- **scATAC-only** (`workflows/signac_pbmc_10k/` + `workflows/snapatac2_pbmc_10k/`) on 10x ATACv2 PBMC 10k
- **Multiome** (`workflows/multiome_pbmc_10k_signac/` + `workflows/multiome_pbmc_10k_muon/`) on 10x Multiome PBMC 10k

The datasets differ: ATACv2 is ATAC-only; Multiome has both RNA and ATAC. Cell-type compositions are comparable:
- scATAC: ~9,800 cells, ~18 clusters (annotated via gene activity)
- Multiome: ~11,300 cells, ~19 clusters (annotated directly from RNA)

**Portfolio value:** two chemistries, two annotation philosophies (proxy versus direct RNA), two integration approaches (single modality versus multi-modal). This breadth covers the modern scATAC workflow space.

## When Multiome is not the right choice

- **Sample size constraint.** Multiome kits produce fewer cells per sample than dedicated protocols. For 100K+ cells, separate scRNA and scATAC preparations are more cost-effective.
- **Depth-focused analysis.** For deep regulatory landscape work (motif footprinting, GWAS variant overlay), scATAC-only provides higher coverage per modality.
- **Legacy data.** If a published dataset is unpaired, bridge integration is more practical than generating new Multiome data.

## References

- Chen et al. 2019, *Nat Biotechnol* — Multiome method paper (first demonstration)
- Hao et al. 2021, *Cell* — Seurat v4 WNN (Multiome integration algorithm)
- Bredikhin et al. 2022, *Genome Biol* — muon (Python multi-modal framework)
- 10x Multiome page: https://www.10xgenomics.com/products/single-cell-multiome-atac-plus-gene-expression
