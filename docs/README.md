# Documentation

Conceptual and methodological documentation for the pipelines in this repository. Written to be read alongside the code, not to replace it.

## scATAC-seq (single-modality)

Sequenced walkthrough of the standard scATAC-seq pipeline, from raw fragments to annotated clusters.

- [01. Biology and assays](scatac/01_biology_and_assays.md) — what Tn5, fragments, and peaks measure; scATAC in the context of scRNA, ChIP, and WGS
- [02. Pipeline overview](scatac/02_pipeline_overview.md) — the six stages, mapping biological question to mathematical operation to interpretation
- [03. Quality control](scatac/03_quality_control.md) — four QC metrics, threshold rationale, filtering decisions
- [04. Normalization and dimensionality reduction](scatac/04_normalization_and_dim_reduction.md) — TF-IDF, LSI, and why the first component is dropped
- [05. Clustering and annotation](scatac/05_clustering_and_annotation.md) — Leiden, gene activity, marker overlay, reference-based label transfer
- [06. R vs Python](scatac/06_r_vs_python.md) — Signac and SnapATAC2 side by side, framework trade-offs, migration paths

Reference workflow implementations:
- R (Signac + Seurat): [`workflows/signac_pbmc_10k/`](../workflows/signac_pbmc_10k/)
- Python (SnapATAC2 + scanpy): [`workflows/snapatac2_pbmc_10k/`](../workflows/snapatac2_pbmc_10k/)

## Multiome (paired scRNA + scATAC)

Extension to multi-modal Multiome data, adding WNN integration and joint clustering.

- [01. Introduction](multiome/01_introduction.md) — why paired > unpaired for regulatory inference, priming, and annotation
- [02. Pipeline overview](multiome/02_pipeline_overview.md) — the eight blocks that extend scATAC's six with WNN
- [03. QC and intersection filter](multiome/03_qc_and_intersection_filter.md) — per-modality metrics and intersection filtering
- [04. Modality processing](multiome/04_modality_processing.md) — SCTransform + PCA for RNA, TF-IDF + LSI for ATAC, ordered correctly
- [05. WNN and joint clustering](multiome/05_wnn_and_joint_clustering.md) — the weighted nearest neighbor algorithm and joint UMAP
- [06. Annotation and downstream](multiome/06_annotation_and_downstream.md) — Azimuth, CellTypist, peak-to-gene linkage, motif enrichment

Reference workflow implementations:
- R (Signac WNN + Seurat): [`workflows/multiome_pbmc_10k_signac/`](../workflows/multiome_pbmc_10k_signac/)
- Python (muon + scanpy): [`workflows/multiome_pbmc_10k_muon/`](../workflows/multiome_pbmc_10k_muon/)

## Suggested reading order

If new to single-cell chromatin analysis:
1. scATAC 01 (biology grounding)
2. scATAC 02 (pipeline overview)
3. scATAC 03 through 05 in order (QC through annotation)
4. scATAC 06 (framework comparison)
5. Multiome 01 through 06 in order (paired modality workflow)

For readers already familiar with scRNA-seq analysis, scATAC 01 is a good starting point; the differences from scRNA are made explicit throughout.

For readers focused on framework choice, scATAC 06 has the direct comparison.
