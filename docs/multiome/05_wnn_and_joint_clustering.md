# 05. WNN (Weighted Nearest Neighbor) and joint clustering

## Why WNN — the failure of simple concatenation

The simplest multi-modal integration is to concatenate the two modality embeddings:

```
Joint embedding = [PC1, PC2, ..., PC30, LSI2, LSI3, ..., LSI50]
                   ← RNA ~30 dim →      ← ATAC ~48 dim →
```

Then run KNN plus Leiden on the joined space.

**Why this is insufficient:**
1. **Scale mismatch.** RNA PCA and ATAC LSI have different value ranges; one modality dominates.
2. **Signal-to-noise varies per cell.** Cell A may be crisp in RNA but noisy in ATAC. Equal weighting is wrong.
3. **Modality bias.** Different feature counts (RNA 30 dim vs ATAC 48 dim) means the higher-dim modality carries more weight.

**WNN's solution:** adaptive per-cell modality weighting.

## The WNN algorithm (Hao et al. 2021 *Cell*)

**Central idea:** For each cell, determine which modality better predicts its identity, then use that modality's weight for that cell.

**Steps:**

**1. Per-modality within-modality KNN:**
- KNN (default K=20) within the RNA embedding for each cell
- KNN within the ATAC embedding for each cell

**2. Cross-modality prediction accuracy:**
- Are cell A's RNA neighbors consistent in ATAC embedding?
- Are cell A's ATAC neighbors consistent in RNA embedding?
- These consistencies form modality-specific "prediction scores" for cell A.

**3. Per-cell modality weight:**
```
w_RNA[cell A] = softmax(prediction_score_RNA[cell A])
w_ATAC[cell A] = 1 - w_RNA[cell A]
```

Examples:
- A cell with clear RNA-side clustering: w_RNA high (0.7), w_ATAC low (0.3)
- A cell in a chromatin priming state: w_ATAC high (0.6), w_RNA low (0.4)

**4. Weighted joint neighbor graph:**
```
joint_distance[A, B] = w_RNA[A] * distance_RNA[A, B] +
                       w_ATAC[A] * distance_ATAC[A, B]
```

Recompute KNN on joint distance → joint neighbor graph.

**5. Downstream:** Leiden clustering plus UMAP on the joint graph.

## Signac API (R)

```r
pbmc <- FindMultiModalNeighbors(
  pbmc,
  reduction.list = list("pca", "lsi"),
  dims.list      = list(1:30, 2:50)
)
```

**Arguments:**
- `reduction.list` — modality reduction names, in order
- `dims.list` — which dimensions of each reduction to use
- Optional: `k.nn` (default 20), `snn.graph.name` (default `"wsnn"`)

**Results attached to SeuratObject:**
- `pbmc[["weighted.nn"]]` — joint KNN graph
- `pbmc[["wsnn"]]` — weighted SNN graph (used for clustering)
- `pbmc@meta.data$RNA.weight` — per-cell RNA modality weight
- `pbmc@meta.data$ATAC.weight` — per-cell ATAC modality weight

## muon API (Python)

```python
mu.pp.neighbors(mdata, key_added='wnn')
```

**Results attached to MuData:**
- `mdata.uns['wnn']` — neighbor graph metadata
- `mdata.obsp['wnn_connectivities']` — sparse connectivity matrix
- `mdata.obsp['wnn_distances']` — sparse distance matrix

**Note:** muon's `mu.pp.neighbors` is not the exact mathematical equivalent of Seurat's WNN; it is a similar cross-modal weighted neighbor approach. Results are close but not identical.

**Multi-modal partition (WNN alternative):**
```python
mu.tl.leiden(mdata, resolution=0.8)   # no neighbors_key
```
This is a different algorithm — partition optimization on per-modality graphs directly, not on a joint graph. Result differs from WNN.

## Joint clustering

**R (post-WNN):**
```r
pbmc <- FindClusters(pbmc, graph.name = "wsnn",
                      algorithm = 3,   # Leiden
                      resolution = 0.8,
                      random.seed = 42)
```

`graph.name = "wsnn"` is critical. The default `pbmc_snn` would use the RNA-only graph, wasting the WNN work.

**Python (on the WNN joint graph):**
```python
sc.tl.leiden(
    mdata,
    obsp='wnn_connectivities',   # joint connectivity matrix
    resolution=0.8,
    key_added='leiden',
    random_state=42,
)
```

The `obsp='wnn_connectivities'` argument directs scanpy to the joint matrix. `mu.tl.leiden(neighbors_key='wnn')` does not work because muon's leiden expects per-modality graphs.

## Joint UMAP

**R:**
```r
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn",
                reduction.name = "wnn.umap",
                reduction.key = "wnnUMAP_")
DimPlot(pbmc, reduction = "wnn.umap", label = TRUE)
```

`nn.name = "weighted.nn"` feeds the WNN's joint KNN into UMAP.

**Python:**
```python
mu.tl.umap(mdata, neighbors_key='wnn')
mu.pl.umap(mdata, color='leiden')
```

## The three-way UMAP comparison

Visualize the benefit of WNN by plotting three UMAPs side by side:

**R:**
```r
p1 <- DimPlot(pbmc, reduction = "umap.rna", label = TRUE) + NoLegend() + ggtitle("RNA UMAP")
p2 <- DimPlot(pbmc, reduction = "umap.atac", label = TRUE) + NoLegend() + ggtitle("ATAC UMAP")
p3 <- DimPlot(pbmc, reduction = "wnn.umap", label = TRUE) + NoLegend() + ggtitle("WNN joint UMAP")
p1 + p2 + p3 + plot_layout(ncol = 3)
```

**Python:**
```python
fig, axes = plt.subplots(1, 3, figsize=(18, 6))
sc.pl.umap(rna, color='leiden', ax=axes[0], show=False, title='RNA UMAP')
sc.pl.umap(atac, color='leiden', ax=axes[1], show=False, title='ATAC UMAP')
mu.pl.umap(mdata, color='leiden', ax=axes[2], show=False, title='WNN joint UMAP')
```

**Critical for Python:** the joint leiden lives in `mdata.obs['leiden']`. Per-modality objects (`rna.obs`, `atac.obs`) do not have it. Sync before plotting:
```python
rna.obs['leiden'] = mdata.obs['leiden'].values
atac.obs['leiden'] = mdata.obs['leiden'].values
```

**Interpretation:**
- **RNA-only:** major cell types separate clearly; fine subtypes may not
- **ATAC-only:** different cluster boundaries; chromatin similarity view. Some cluster splits absent from RNA
- **WNN joint:** groups confirmed by both modalities; the most robust cluster structure. Rare states (priming, transitions) become visible.

If new clusters appear in WNN that neither RNA nor ATAC show independently, the Multiome value is proven.

## Resolution parameter for WNN

Resolution selection is more critical in WNN than in single-modality analysis. Two modalities add information; cluster count naturally grows:

- **0.3–0.5:** major cell types (5–8 clusters)
- **0.8:** standard for Multiome PBMC 10k, 15–20 clusters
- **1.5:** fine subtype detection, 25–35 clusters

Our workflow uses `resolution = 0.8`, yielding 19 clusters on WNN joint.

## Modality weight interpretation (advanced)

**R post-WNN:**
```r
pbmc$RNA.weight   # per-cell RNA weight
pbmc$ATAC.weight  # per-cell ATAC weight

# Cluster-level average
pbmc@meta.data %>% group_by(seurat_clusters) %>%
                    summarise(mean_rna_w = mean(RNA.weight))
```

**Biological interpretation:**
- High-RNA-weight clusters: transcriptional identity dominant (active cell types)
- High-ATAC-weight clusters: regulatory landscape dominant, transcription not yet reflecting the state (priming, transition)

This clue points to which clusters occupy which biological states. Optional advanced analysis.

## multiVI — the deep learning alternative

WNN is classical. The modern deep learning alternative is multiVI (Ashuach 2023).

**multiVI approach:**
- Separate VAE encoders per modality
- Joint latent space (shared z)
- Native batch correction
- Better rare cell type detection
- GPU recommended

**Setup:**
```bash
pip install scvi-tools
```

```python
import scvi
scvi.model.MULTIVI.setup_mudata(mdata, batch_key='sample')
model = scvi.model.MULTIVI(mdata)
model.train()
mdata.obsm['X_multivi'] = model.get_latent_representation()
# KNN + Leiden + UMAP on this latent
```

**Trade-off:**
- WNN: simple, fast, no GPU, mature
- multiVI: complex, requires ML expertise, GPU recommended, better at atlas scale

This pipeline uses WNN as sufficient. multiVI is future sprint territory.

## Common pitfalls

1. **`dims.list` order.** If RNA is the first element of `reduction.list`, use RNA dims in the first element of `dims.list`. Swapping means wrong dims per modality.
2. **`graph.name` choice.** `FindClusters(graph.name="pbmc_snn")` uses the RNA-only graph. Post-WNN, use `"wsnn"`.
3. **muon `neighbors_key` KeyError.** muon's leiden expects per-modality; for joint WNN clustering use `sc.tl.leiden(obsp='wnn_connectivities')`.
4. **Sync leiden across modalities.** Python joint clustering populates only `mdata.obs['leiden']`; per-modality obs need manual sync for plotting.
5. **UMAP shape stochastic.** Same seed reproduces layout across runs, but changes between runs of different seeds. Cluster identity is preserved; orientation may change.

## References

- Hao et al. 2021, *Cell* — WNN original paper
- Ashuach et al. 2023, *Nat Methods* — multiVI deep integration
- Bredikhin et al. 2022, *Genome Biol* — muon multi-modal framework
- Seurat WNN vignette: https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis
- Signac Multiome vignette: https://stuartlab.org/signac/articles/pbmc_multiomic
