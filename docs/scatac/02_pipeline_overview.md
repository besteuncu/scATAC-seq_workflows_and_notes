# 02. Pipeline overview — biology, math, and meaning

## The six stages

A single-cell ATAC-seq analysis proceeds through six logical stages:

```
Raw fragments → Peak calling → Cell × peak matrix → QC filter
                                                        |
                            Cell annotation ← Clustering ← Dim reduction
```

Each stage answers a specific biological question with a specific mathematical operation:

| Stage | Biological question | Mathematical operation | What it tells us |
|-------|--------------------|-----------------------|------------------|
| 1. Peak calling | Where do fragments concentrate? | MACS2 density detection | Catalog of regulatory elements |
| 2. Quality control | Which cells are real nuclei? | TSS + fragment count thresholds | Trustworthy cell subset |
| 3. TF-IDF | How to correct depth and prevalence bias? | NLP-style weighting | Comparable cell vectors |
| 4. LSI (SVD) | How to compress the feature space? | SVD to 30 dims, drop dim 1 | Cell embedding (30D) |
| 5. Clustering | Which cells are similar? | Leiden on KNN graph | Discrete cell groups |
| 6. Annotation | What cell type is each cluster? | Gene activity plus reference | Biological identity |

## Stage 1: peak calling

**Question:** Where in the genome does Tn5 tagmentation signal concentrate?

**Input:** Fragment file with all cells pooled.

**Operation:** MACS2 computes per-position fragment density, tests for enrichment against background, and calls peaks with ATAC-specific parameters (`--nomodel --shift -75 --extsize 150`).

**Output:** `peaks.bed` with 50,000–200,000 peaks.

**Meaning:** A catalog of accessible regulatory elements (enhancers, promoters, TFBS) in this dataset. The cell-of-origin information is not yet resolved.

In this repository we rely on Cellranger-ATAC's built-in peak calling, which produces `filtered_peak_bc_matrix.h5` directly. Custom peak calling is left for future work.

## Stage 2: quality control

**Question:** Of Cellranger's roughly 10,000 candidate cells, which represent real nuclei?

**Input:** Cell × peak matrix plus per-cell metadata (`singlecell.csv`).

**Operation:** Compute per-cell metrics — TSS enrichment, nucleosome signal, fraction of reads in peaks (FRiP), total fragment count. Then apply threshold intersection.

**Output:** Filtered subset of about 8,500 to 9,500 cells.

**Meaning:** Cell integrity indicators. Low TSS enrichment implies ambient DNA or damaged nuclei. High nucleosome signal implies overtagmentation. Low FRiP implies weak biological signal relative to background.

The filter timing is critical: all downstream analyses use the filtered subset.

## Stage 3: TF-IDF normalization

**Question:** How can cells with different sequencing depths be compared fairly?

**Input:** Filtered sparse binary-like cell × peak matrix.

**Operation:** Term Frequency-Inverse Document Frequency, borrowed from information retrieval:
- **TF** (Term Frequency): `count[i,j] / total_count[i]` normalizes depth
- **IDF** (Inverse Document Frequency): `log(N_cells / N_cells_with_peak[j])` down-weights ubiquitous peaks
- **Combined:** `log(1 + TF × IDF × 1e4)` stabilizes variance

**Output:** Depth-normalized, informative-peak-weighted matrix.

**Meaning:** Each cell's chromatin profile is now fair to compare. Deep-sequenced cells no longer dominate. Ubiquitous promoter peaks (accessible in every cell) receive small weights; cell-type-specific enhancer peaks receive large weights.

The NLP analogy is not superficial. A sparse binary document × term matrix has the same structural properties as a sparse binary cell × peak matrix. The information retrieval community solved this class of problem in the 1980s (LSA + IDF weighting), and scATAC-seq adopted the recipe directly.

## Stage 4: LSI (Latent Semantic Indexing)

**Question:** How to compress a 100K-peak space into an interpretable low-dimensional embedding?

**Input:** TF-IDF normalized matrix.

**Operation:** Singular Value Decomposition (SVD). Compute 30 leading singular vectors and construct a per-cell embedding.

**Critical detail:** The largest singular vector (LSI-1) typically correlates strongly with sequencing depth. TF-IDF partially normalizes depth but cannot remove it fully. LSI-1 is dropped, and downstream analysis uses `dims = 2:30`.

**Output:** Cell × 30 embedding matrix (or 29, after dropping LSI-1).

**Meaning:** Each cell has biological latent coordinates. Neighbors in this 30D space are biologically similar cells.

The framework difference: Signac requires manual `dims = 2:30`; SnapATAC2 uses spectral embedding which handles depth confound internally.

## Stage 5: clustering

**Question:** Into how many groups do cells fall, and which cell belongs to which group?

**Input:** LSI embedding.

**Operation:**
1. **KNN graph:** each cell connects to its 20 nearest neighbors in the 30D space.
2. **Leiden algorithm:** community detection on the graph, maximizing modularity (dense within groups, sparse between).

**Output:** Cluster label per cell.

**Meaning:** Discrete cell groups. For 10x PBMC 10k, typical output is 12–20 clusters at resolution 0.5–1.0.

The resolution parameter controls granularity: lower values give fewer, broader clusters; higher values give finer subclustering.

## Stage 6: UMAP visualization

UMAP (Uniform Manifold Approximation and Projection) is not an analysis step. It is a two-dimensional projection for visualization only.

**Note:** Cluster decisions are made in the 30D LSI space. UMAP preserves neighborhood structure but not distance. Its 2D shape is stochastic and reproducibility depends on seed. Cross-run comparison is done via cluster centroid signatures, not UMAP overlay.

## Stage 7: cell type annotation

**Question:** Which cell type does each cluster correspond to?

**Approach 1 — Gene activity plus canonical markers:**
1. Sum peak accessibility over each gene's body plus promoter (upstream 2 kb) to create a proxy expression matrix.
2. Normalize and log-transform.
3. Overlay canonical marker genes (CD3E, MS4A1, CD14) on the UMAP.
4. Assign labels based on marker enrichment per cluster.

Gene activity is a proxy, not real expression: longer genes are over-estimated, distal enhancer information is lost. For immune markers with strong promoter accessibility (like CD3E, CD14), the proxy works well enough.

**Approach 2 — Reference-based label transfer (recommended):**
Use Azimuth (R) or CellTypist (Python) to transfer labels from a reference atlas (Hao 2021 PBMC atlas, Domínguez 2022 immune atlas). Faster, reproducible, and provides confidence scores.

## Information flow across stages

Each stage builds on the previous and loses raw information while gaining meaning:

```
Fragments (billions of signals, no interpretation)
    |
Peaks (~100K regions, "regulatory potential" interpretation)
    |
Clusters (~18 groups, "different cell states" interpretation)
    |
Cell types (7-10 populations, "biological composition" interpretation)
```

Cell counts stay roughly constant after Stage 4; only labels are added. Compression is aggregation of interpretation, not loss of information.

## Beyond this pipeline

The six-stage core is only the beginning of what scATAC-seq analysis enables. Advanced analyses include:
- **Differential accessibility** between clusters
- **Motif enrichment** in cluster-marker peaks
- **ChromVAR** per-cell TF activity
- **Peak-to-gene linking** for enhancer-gene inference
- **TF footprinting** for binding inference (a ChIP-seq alternative)
- **Trajectory analysis** for cell state transitions

Those are the payoff for a robust core pipeline.

## R versus Python

Both frameworks implement the same six stages under different names:

| Stage | Signac (R) | SnapATAC2 (Python) |
|-------|------------|--------------------|
| Load | `Read10X_h5` + `CreateChromatinAssay` | `snap.pp.import_fragments` |
| QC | `NucleosomeSignal` + `TSSEnrichment` | `snap.metrics.tsse` |
| Filter | `subset(..., subset = ...)` | `snap.pp.filter_cells` |
| Normalize | `RunTFIDF` | `snap.pp.select_features` (implicit) |
| Dim red | `RunSVD` → LSI dims 2:30 | `snap.tl.spectral` |
| Cluster | `FindNeighbors` + `FindClusters` | `snap.pp.knn` + `snap.tl.leiden` |
| UMAP | `RunUMAP` | `snap.tl.umap` |
| Annotate | `GeneActivity` + `FeaturePlot` | `snap.pp.make_gene_matrix` + `sc.pl.umap` |

For a detailed side-by-side comparison see [06_r_vs_python.md](06_r_vs_python.md).

## References

- Cusanovich et al. 2018, *Cell* — scATAC atlas and LSI methodology.
- Stuart et al. 2021, *Nat Methods* — Signac.
- Zhang et al. 2024, *Nat Methods* — SnapATAC2.
