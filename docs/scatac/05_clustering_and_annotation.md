# 05. Clustering and annotation — from math to biology

## Goal

Group cells in the LSI embedding into clusters, then assign a biological cell type name to each cluster.

## Part A: clustering — mathematical groups

### Working in 30D, not 2D

Cluster decisions happen in the LSI 30-dimensional space, not the UMAP 2D projection. UMAP is visualization only. The 30D space captures real biological variance; UMAP loses information for the sake of a 2D picture.

### KNN graph construction

For each cell, identify the K nearest neighbors (default K = 20) in LSI space:

```r
pbmc <- FindNeighbors(pbmc, reduction = "lsi", dims = 2:30)
```

Output: each cell becomes a node; edges connect to its 20 nearest neighbors. Signac also computes a **Shared Nearest Neighbor (SNN)** graph: edge weights reflect how many neighbors two cells share.

Why SNN over raw KNN: raw KNN's "top 20 neighbors" is a hard cutoff. SNN smooths this — cells with overlapping neighborhoods get stronger edges. Cluster boundaries become more organic.

### Leiden algorithm

Community detection on the SNN graph, maximizing modularity:

- **Modularity** measures the quality of a graph partition: dense connectivity within groups, sparse between.
- Modularity > 0.4 is good, > 0.6 is excellent.

```r
pbmc <- FindClusters(pbmc, algorithm = 3, resolution = 0.8)
# algorithm = 3 → Leiden (recommended)
# algorithm = 1 → Louvain (classic, older)
```

Leiden (Traag et al. 2019) fixes a Louvain bug — well-connected communities are guaranteed. It is the current default.

### The resolution parameter

The single most critical knob. Smaller resolution gives fewer clusters; larger gives more:

| Resolution | Cluster count | Biological interpretation |
|-----------|---------------|--------------------------|
| 0.1 | 3–5 | Broadest: T, B, myeloid |
| 0.3 | 6–10 | CD4, CD8, B, NK, mono, DC |
| 0.5 | 10–15 | plus subtypes |
| **0.8** | **15–20** | **typical for PBMC 10k** |
| 1.5 | 25–35 | fine states |
| 3.0 | 50+ | over-clustered, noise dominant |

This pipeline uses `resolution = 0.8`, producing 18 clusters.

Resolution selection guide:
- Broad cell type identification: 0.3–0.5
- Subtype detection: 0.8–1.2
- State discovery: 1.5–2.5
- Verification: cluster-specific marker genes should be identifiable

### Cluster stability

How real is a cluster? Check:
1. **Multi-resolution consistency.** The `clustree` package (Zappia 2018) visualizes how clusters split as resolution rises. Stable clusters persist across multiple resolutions.
2. **Marker sharpness.** A cluster with no clearly enriched markers may be an over-split of an actual cluster.
3. **Cross-framework consistency.** If R Signac and Python SnapATAC2 recover similar cluster counts, the biology is robust.

In this repository R found 18 clusters and Python found 13–15. Same biology, different granularity. Both are valid.

## Part B: UMAP — visualization

UMAP (Uniform Manifold Approximation and Projection) reduces the 30D embedding to 2D for plotting.

```r
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 2:30, seed.use = 42)
DimPlot(pbmc, label = TRUE) + NoLegend()
```

Important: UMAP is stochastic. Each seed produces a different layout. Cluster identity is preserved; visual appearance differs. Do not cross-compare UMAP shapes across runs — cluster centroid signatures are the reproducible units.

Distance interpretation: UMAP preserves neighborhood, not distance. Two clusters far apart on the UMAP are not necessarily biologically distant. To assess biological distance, examine gene expression comparison directly.

## Part C: cell type annotation

Give each cluster a biological label. Three approaches, in order of sophistication.

### Approach 1: gene activity plus canonical markers

The pipeline we implemented.

**Steps:**

1. Compute gene activity — a proxy for expression from peaks near each gene:
   ```r
   gene.activities <- GeneActivity(pbmc)
   ```
   The function sums peak accessibility over each gene's body plus promoter (upstream 2 kb).

2. Normalize and log-transform:
   ```r
   pbmc[["RNA"]] <- CreateAssayObject(counts = gene.activities)
   pbmc <- NormalizeData(pbmc, assay = "RNA")
   ```

3. Overlay canonical markers on UMAP:
   ```r
   markers <- c("MS4A1", "CD3E", "CD8A", "CD4", "NKG7", "CD14", "LYZ", "FCGR3A")
   FeaturePlot(pbmc, features = markers, ncol = 4, max.cutoff = "q95")
   ```

4. Assign labels manually based on marker enrichment.

**Caveats of gene activity as proxy:**
- Longer genes are over-estimated (more peaks fall within them).
- Distal enhancers (>2 kb) are excluded.
- Chromatin priming can look like expression when no transcription is happening.

For immune markers (CD3E, MS4A1, CD14), with strong promoter accessibility, the proxy works well. For fine subtypes, it degrades.

### Canonical PBMC markers

| Marker | Cell type |
|--------|-----------|
| **MS4A1** (CD20) | B cell |
| **CD3E, CD3D** | Pan-T cell |
| **CD4** | CD4 T (also on monocytes) |
| **CD8A, CD8B** | CD8 T |
| **NKG7, KLRD1** | NK cell |
| **CD14** | Classical monocyte |
| **FCGR3A** (CD16) | CD16 non-classical monocyte and NK |
| **LYZ** | Monocyte and DC |
| **FCER1A** | Dendritic cell |
| **PPBP** | Platelet |
| **CD34** | HSPC |

### Approach 2: reference-based label transfer (recommended)

Use Azimuth (R) or SingleR (R) or CellTypist (Python) to transfer labels from a curated reference atlas.

**Azimuth PBMC (R):**
```r
library(Azimuth)
pbmc <- RunAzimuth(pbmc, reference = "pbmcref")
DimPlot(pbmc, group.by = "predicted.celltype.l2", label = TRUE)
```

Reference: Hao 2021 PBMC atlas (165K cells, L1/L2/L3 hierarchical labels). Result fields include `predicted.celltype.l1/l2/l3` and confidence scores.

**Advantages:**
- Reproducible: same reference produces same labels
- Fine subtype resolution (30+ cell types at L2)
- Confidence score per cell

**Disadvantages:**
- Requires reference installation and download (~500 MB one-time)
- Reference biology may not match query biology (disease, treatment, different tissue)

For scATAC-only workflows, Azimuth operates on gene activity, so proxy limitations apply.

### Approach 3: ChromVAR TF motif activity (scATAC-native)

Skip the RNA proxy entirely by using TF motif accessibility as the identity signal.

```r
library(chromVAR)
pbmc <- AddMotifs(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = ...)
pbmc <- RunChromVAR(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)
DefaultAssay(pbmc) <- "chromvar"
FeaturePlot(pbmc, features = c("GATA3", "TBX21", "EBF1"))
```

Lineage-defining TFs give identity:

| TF | Cell type indicator |
|----|--------------------|
| SPI1 (PU.1) | Myeloid |
| EBF1 | B cell |
| TBX21 (T-bet) | Th1, NK |
| GATA3 | Th2, NK |
| FOXP3 | Treg |
| RUNX1 | Hematopoietic |

ChromVAR is peak-based and gene-independent, matching the ATAC-first philosophy. It is future work for this repository.

## Multi-resolution refinement

Real workflows often use two passes:
1. **Coarse resolution (0.3–0.5)** to separate major cell types.
2. **Per-cluster subclustering** — extract one cluster (e.g. all T cells), subset, and rerun LSI + UMAP + Leiden. Fine subtypes (naive CD4, memory CD4, effector CD8) emerge.

This iterative refinement is standard for immune datasets. This repository's pipeline uses a single pass for simplicity.

## After annotation

With cell types identified, downstream analyses become possible:
- **Differential accessibility:** which peaks differ between cell types? Cluster-marker peaks generate biological hypotheses.
- **Regulatory network inference:** motif enrichment in cluster-marker peaks identifies cell-type-specific TFs.
- **Trajectory analysis:** naive-to-effector transitions can be mapped in the LSI embedding (Monocle3, scFates).
- **GWAS overlay:** non-coding disease variants can be mapped to specific cell types via peak overlap (Corces 2020 pattern).

## What this pipeline delivers

Final output:
- ~9,850 cells (R Signac) or ~10,149 cells (Python SnapATAC2)
- ~18 clusters (Signac) or ~13–15 clusters (SnapATAC2)
- Canonical marker overlay identifies:
  - T cell populations (CD4, CD8, naive, memory): ~50%
  - B cells: ~10%
  - Monocytes (CD14 classical + CD16 non-classical): ~25%
  - NK cells: ~5–10%
  - Rare types (DC, pDC, HSPC): <5%

## Common pitfalls

1. **Over-splitting from high resolution.** At resolution 1.5+, 50+ clusters emerge, most doublets or noise.
2. **Marker ambiguity.** CD4 is expressed in both CD4 T cells and monocytes. Combined markers (CD3E+ and CD4+ = CD4 T; CD3E- and CD4+ = monocyte) resolve.
3. **Over-interpreting gene activity.** High GAPDH activity is not a cell type signal; it means the gene body is broadly accessible.
4. **UMAP shape as biology.** Cluster distances in UMAP are visualization artifacts. Verify by inspecting the embedding directly.
5. **Reference bias.** An Azimuth PBMC reference may misclassify cells from patients with cancer or immune disease. Reference-query domain match matters.

## References

- Traag et al. 2019 — Leiden algorithm
- McInnes et al. 2018 — UMAP
- Stuart et al. 2021, *Nat Methods* — Signac annotation vignette
- Hao et al. 2021, *Cell* — Azimuth PBMC reference
- Schep et al. 2017, *Nat Methods* — ChromVAR
