# 06. Cell annotation and downstream — Multiome's payoff

## Annotation — Multiome's easy path

In scATAC-only pipelines annotation was difficult: gene activity proxy, long-gene bias, marker gene promoter accessibility not always tracking expression. Multiome simplifies this completely: the RNA modality provides canonical marker gene expression directly.

## Approach 1: automated reference-based annotation (recommended)

This is the preferred method — reproducible, evidence-based, and provides confidence scores.

### R: Azimuth PBMC reference

Azimuth (Hao 2021 PBMC atlas, 165K cells) provides hierarchical L1/L2/L3 labels. Multiome-friendly: WNN and Azimuth were designed together.

**Install (once):**
```r
renv::install("satijalab/azimuth")
```

**Use:**
```r
library(Azimuth)

# Azimuth runs on RNA modality (SCT is normalized)
DefaultAssay(pbmc) <- "SCT"

# Label transfer from PBMC reference
pbmc <- RunAzimuth(pbmc, reference = "pbmcref")

# Result fields:
# - pbmc$predicted.celltype.l1 (broad: T, B, Mono, NK, DC, ...)
# - pbmc$predicted.celltype.l2 (subtypes: CD4 Naive, CD8 TEM, ...)
# - pbmc$predicted.celltype.l3 (fine subtypes)
# - pbmc$predicted.celltype.l2.score (confidence 0-1)

# Assign cell_type from Azimuth L2
pbmc$cell_type <- pbmc$predicted.celltype.l2

# Visualize
DimPlot(pbmc, reduction = "wnn.umap",
        group.by = "cell_type", label = TRUE, repel = TRUE) +
        NoLegend()
```

On first call, Azimuth downloads the PBMC reference (~500 MB, one-time). Subsequent calls use the cache.

### Python: CellTypist Immune_All_Low

CellTypist (Domínguez 2022 *Science*) provides pretrained ML models for immune tissues. The `Immune_All_Low.pkl` model has 32 immune cell type labels.

**Install (once):**
```bash
pip install celltypist
```

**Use:**
```python
import celltypist
from celltypist import models

# Download PBMC-tuned model (~10 MB, cached after first call)
models.download_models(model='Immune_All_Low.pkl', force_update=False)
pbmc_model = models.Model.load(model='Immune_All_Low.pkl')

# Annotate RNA modality. CellTypist needs log1p-normalized data (not scaled).
# If Block 5 saved rna.raw before sc.pp.scale, CellTypist reads from .raw.X.
if rna.raw is None:
    rna_for_annot = rna.copy()
    rna_for_annot.X = rna_for_annot.layers['counts'].copy()
    sc.pp.normalize_total(rna_for_annot, target_sum=1e4)
    sc.pp.log1p(rna_for_annot)
    input_ad = rna_for_annot
else:
    input_ad = rna

predictions = celltypist.annotate(
    input_ad,
    model=pbmc_model,
    majority_voting=True,           # cluster-consensus
    over_clustering=rna.obs['leiden'],   # use our Leiden clusters
)

# Attach labels
rna.obs = rna.obs.join(predictions.predicted_labels)
mdata.obs['cell_type'] = rna.obs['majority_voting'].values
atac.obs['cell_type'] = mdata.obs['cell_type'].values

# Visualize
mu.pl.umap(mdata, color='cell_type', legend_loc='on data',
           legend_fontsize=8)
```

`majority_voting=True` plus `over_clustering=rna.obs['leiden']`: within each cluster, the majority label becomes the cluster's label. This enforces cluster-level consistency.

### Comparison

| Dimension | Azimuth (R) | CellTypist (Python) |
|-----------|-------------|---------------------|
| Reference | 165K PBMC atlas (Hao 2021) | Immune_All_Low pretrained (356K) |
| L1/L2/L3 hierarchy | Yes | Typically single L2-equivalent |
| Confidence score | `.score` field | `probability_matrix` |
| Install time | ~5 min (GitHub + reference) | 30 sec (pip install) |
| Reference size | ~500 MB (one-time) | ~10 MB per model |
| Multiome awareness | Designed with WNN | RNA-only (still effective) |
| Rare cell types | pDC, HSPC, MAIT included | Similar |

Both are scientific standards and work well for PBMC. Azimuth integrates more naturally with WNN in R; CellTypist has a lighter Python install.

## Approach 2: cluster-averaged marker inspection (optional verification)

Before or after automated annotation, sanity-check clusters against canonical markers.

**R:**
```r
avg_expr <- AverageExpression(
  pbmc,
  features = markers,
  group.by = "seurat_clusters",
  assays = "SCT"
)$SCT

round(avg_expr, 2)
```

Example output:
```
        0     1     2     3     4     5    ...
MS4A1  0.02  0.03  0.05  0.02  0.04  3.85  ...
CD3E   0.05  4.20  3.85  2.98  0.15  0.08  ...
CD14   4.02  0.10  0.15  0.08  0.20  0.15  ...
NKG7   0.30  1.20  0.85  4.12  4.85  0.10  ...
```

Each column is a cluster; rows are markers. Highest marker per cluster hints at its cell type. Cluster 0 has highest CD14 → monocyte. Cluster 5 has highest MS4A1 → B cell.

**Python equivalent:**
```python
import pandas as pd

cluster_means = pd.DataFrame(
    index=markers,
    columns=sorted(rna.obs["leiden"].unique()),
)
for cluster in cluster_means.columns:
    cells = rna.obs["leiden"] == cluster
    for gene in markers:
        cluster_means.at[gene, cluster] = rna[cells, gene].X.mean()

print(cluster_means.round(2))
```

## Downstream 1: peak-to-gene linkage

This is Multiome's main analytical strength — connecting regulatory elements to the genes they control.

```r
DefaultAssay(pbmc) <- "ATAC"

# Compute region stats (required for LinkPeaks)
pbmc <- RegionStats(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)

# Link peaks to specific genes
pbmc <- LinkPeaks(
  pbmc,
  peak.assay       = "ATAC",
  expression.assay = "SCT",   # RNA (SCT normalized)
  genes.use        = c("MS4A1", "CD8A", "CD14", "NKG7")
)

# View links
Links(pbmc)
```

**Logic:** For each gene, correlate peak accessibility with gene expression across a 500 kb window. High correlation plus short distance suggests a regulatory link.

**Multiome advantage:** paired values in the same cells. scATAC-only relied on reference scRNA integration, which introduces integration error. Multiome computes the link directly.

**Visualization:**
```r
CoveragePlot(pbmc, region = "CD8A",
             features = "CD8A",
             expression.assay = "SCT",
             extend.upstream = 5000,
             extend.downstream = 5000)
```

Output: CD8A locus, flanking peaks, linkage arcs, and per-cluster expression side by side. Regulatory landscape and expression aligned.

## Downstream 2: differential accessibility per cluster

Cluster-marker peaks generate biological hypotheses: which enhancer is specifically open in this cell type?

```r
DefaultAssay(pbmc) <- "ATAC"

# CD8 T effector cluster vs other T cells
da_peaks <- FindMarkers(
  pbmc,
  ident.1 = "CD8 T effector",
  ident.2 = c("CD4 T naive", "CD4 T memory"),
  test.use = "LR",
  latent.vars = "nCount_ATAC",
  min.pct = 0.05
)

head(da_peaks, 20)
```

**Test choice:** `test.use = "LR"` (logistic regression) plus `latent.vars = "nCount_ATAC"` addresses the depth confound. The modern alternative is pseudobulk plus DESeq2 (Squair 2021 argument).

Next step: motif enrichment on differential peaks (ChromVAR) — which TFs are unique to CD8 T effectors?

## Downstream 3: motif enrichment and TF activity

Per-cell TF motif accessibility with ChromVAR:

```r
library(chromVAR)
library(JASPAR2020)

pfm <- getMatrixSet(JASPAR2020, opts = list(species = 9606))
pbmc <- AddMotifs(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
pbmc <- RunChromVAR(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)

DefaultAssay(pbmc) <- "chromvar"
FeaturePlot(pbmc, features = c("EBF1", "TBX21", "GATA3"),
            reduction = "wnn.umap")
```

**Extra value in Multiome:** TF motif activity (ATAC) plus TF mRNA expression (RNA) side by side. High correlation → active TF. Discordant → post-translational regulation or priming state.

## Downstream 4: trajectory analysis

For cell state transitions:
- R: Monocle3, slingshot
- Python: scFates, CellRank, scVelo

In Multiome trajectory analysis, RNA velocity plus chromatin dynamics can be examined together. Example question: in a naive-to-effector T cell transition, which peaks open first and which genes are expressed after?

## Beyond this pipeline — roadmap

The Multiome tutorial covers the basic pipeline. Advanced sprint topics:

1. **Detailed cell type annotation** (Azimuth L1/L2/L3 hierarchy comparison)
2. **Systematic peak-to-gene linkage** (all genes, network inference)
3. **Motif enrichment in cluster-marker peaks** (regulatory network)
4. **Paired differential accessibility and expression analysis**
5. **Trajectory inference** (developmental cell state transitions)
6. **Multi-sample integration** (multiple Multiome samples)
7. **Batch correction** (Harmony, scVI)

Each is a deep dive on its own. The basic pipeline is the substrate.

## R vs Python tool ecosystem (annotation and downstream)

| Task | R (Seurat + Signac) | Python (muon + scanpy) |
|------|--------------------|------------------------|
| Manual annotation | `RenameIdents` | `obs.map(dict)` |
| Reference-based | Azimuth (mature) | scArches (deep learning) |
| Peak-to-gene | `LinkPeaks` (built-in) | Manual Pearson corr or scenicplus |
| Motif enrichment | ChromVAR (built-in) | pycisTopic + pycisTarget |
| Differential accessibility | `FindMarkers(test.use="LR")` | scanpy's rank_genes_groups |
| Trajectory | Monocle3, slingshot | scFates, CellRank, scVelo |

The R Multiome ecosystem is more mature and integrated. Python is more modern and ML-heavy, GPU-friendlier. Choice depends on research question and team expertise.

## Common pitfalls

1. **DefaultAssay mismatch.** RNA markers are in `SCT` assay; ATAC peaks in `ATAC` assay. Wrong assay makes FeaturePlot blank.
2. **Cluster label sync (Python).** `mdata.obs['leiden']` exists; `rna.obs['leiden']` does not by default. Manual sync before visualization.
3. **`LinkPeaks` slow.** Running on all genes is expensive. Prefer `genes.use = c(...)` with a targeted list.
4. **Azimuth reference mismatch.** A healthy PBMC reference misclassifies disease samples. Verify reference-query biology match.
5. **Motif enrichment scale.** ChromVAR + JASPAR works for small datasets but does not scale to 100K+ cells. Consider SnapATAC2 + peakVI for atlas-scale.

## References

- Ma et al. 2020, *Cell* — SHARE-seq (paired multiome pioneer)
- Cao et al. 2018, *Science* — sci-CAR (first paired multi-omic method)
- Hao et al. 2021, *Cell* — Seurat v4 (Azimuth PBMC + WNN)
- Domínguez et al. 2022, *Science* — CellTypist reference
- Schep et al. 2017, *Nat Methods* — ChromVAR TF motif activity
- Corces et al. 2020, *Nat Genet* — Multiome + GWAS variant overlay
- Bravo González-Blas et al. 2023, *Nat Methods* — SCENIC+ (Python multi-modal GRN)
