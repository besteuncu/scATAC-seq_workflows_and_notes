# 04. Modality processing — RNA and ATAC handled separately

## Why separate processing

The two modalities have fundamentally different statistical natures:
- **RNA:** count magnitude is informative (NB distribution), features count ~20K (genes)
- **ATAC:** sparse binary-like (~95% zero), features count 100K+ (peaks)

Applying the same normalization and dimensionality reduction to both would be wrong:
- TF-IDF on RNA makes no sense (magnitude information is present; normalize scales, not signals)
- PCA on ATAC fails (sparse binary matrices break the centering assumption)

Each modality goes through its own classical pipeline, and WNN combines them later.

## RNA processing — SCTransform + PCA

### SCTransform (R)

**What it does:** Poisson-based variance stabilization. Combines depth normalization, variance stabilization, and top variable feature selection in a single function.

```r
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc, verbose = FALSE)
```

**Mathematical background** (Hafemeister and Satija 2019, *Genome Biol*):
1. Regress out cell sequencing depth via Pearson residual regression
2. Model overdispersion per gene
3. Compute standardized residuals
4. Select top variable features (default 3000)

**Classic alternative:**
```r
pbmc <- NormalizeData(pbmc)          # log(count/depth * 10000 + 1)
pbmc <- FindVariableFeatures(pbmc)   # dispersion-based, top 2000
pbmc <- ScaleData(pbmc)              # mean 0, variance 1
```

**Why SCTransform is preferred:** modern approach, better depth-confound handling, and downstream PCA is more "biological". The classic pipeline is older but functional.

**Critical detail:** SCTransform creates a new assay called `SCT`. Downstream PCA runs on the `SCT` assay (SCTransform switches DefaultAssay automatically).

### highly_variable_genes (Python)

The scanpy `seurat_v3` flavor does similar work:

```python
# HVG on RAW counts (order-critical!)
sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor='seurat_v3', layer='counts')

# Then normalize + log + scale
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)
sc.pp.scale(rna, max_value=10)
```

**Why `layer='counts'`:** `seurat_v3` runs loess variance regression on raw counts. If normalize_total is applied first, non-integer values break loess. `layer='counts'` explicitly reads the raw layer regardless of what `.X` currently contains.

**Alternative:** `flavor='cell_ranger'` or default (`'seurat'`) — simpler but not as robust as `seurat_v3`.

### PCA (both stacks)

Standard and unchanged.

**R:**
```r
pbmc <- RunPCA(pbmc)   # default 50 PCs
ElbowPlot(pbmc, ndims = 50)   # elbow typically 20-30
```

**Python:**
```python
sc.tl.pca(rna, n_comps=50)
sc.pl.pca_variance_ratio(rna, n_pcs=50)
```

**PC interpretation:** The first 5–10 PCs typically capture major cell type separations. Beyond that: fine subtypes or noise.

Example PBMC 10k Multiome PC1–PC5:
- PC1: T cell (RPS/RPL ribosomal) vs Monocyte (VCAN, LYZ)
- PC2: NK/CTL (GNLY, NKG7) vs B (MS4A1, PAX5)
- PC3: Naive T (LEF1, CCR7) vs Effector CTL (GZMH, KLRD1)
- PC4: pDC (TCF4, LILRA4, IRF8) discrimination
- PC5: CD16+ non-classical monocyte

Directly biology-driven PCs indicate high data quality.

## ATAC processing — TF-IDF + LSI

### TF-IDF + LSI (R)

Same as scATAC-only pipeline (see scATAC docs). Detail in [../scatac/04_normalization_and_dim_reduction.md](../scatac/04_normalization_and_dim_reduction.md).

```r
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")
pbmc <- RunSVD(pbmc)
```

**Multiome-specific note:** WNN downstream uses `dims = 2:50` (dropping LSI-1 depth confound). This is more dims than scATAC-only's typical `2:30`. Extra dims add information WNN can use.

### TF-IDF + LSI (Python)

```python
# muon's atac module
ac.pp.tfidf(atac, scale_factor=1e4)
ac.tl.lsi(atac)

# Manual drop of LSI-1 (muon does not auto-handle)
atac.obsm['X_lsi'] = atac.obsm['X_lsi'][:, 1:]
atac.varm["LSI"] = atac.varm["LSI"][:, 1:]
atac.uns["lsi"]["stdev"] = atac.uns["lsi"]["stdev"][1:]
```

**Why manual drop:** muon's non-spectral LSI implementation does not handle the depth confound automatically the way SnapATAC2's spectral does. Slicing `[:, 1:]` manually skips it.

## Modality-only UMAP (optional comparison)

Producing per-modality UMAPs before WNN is useful for comparison:

**R:**
```r
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:30,
                reduction.name = "umap.rna", reduction.key = "rnaUMAP_")
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 2:50,
                reduction.name = "umap.atac", reduction.key = "atacUMAP_")
```

**Python:**
```python
sc.pp.neighbors(rna, n_pcs=30)
sc.tl.umap(rna)   # rna.obsm['X_umap']

sc.pp.neighbors(atac, use_rep="X_lsi")
sc.tl.umap(atac)   # atac.obsm['X_umap']
```

**Interpretation:**
- RNA-only UMAP: major cell types separate cleanly; fine subtypes sometimes blur
- ATAC-only UMAP: chromatin similarity view; different cluster boundaries
- Neither is "ground truth"; WNN joint UMAP provides the richest view

## Modality comparison — biological interpretation

Multiome's real value is that the two modalities capture different information per cell.

**RNA excels at:**
- Actively expressed genes (high transcription)
- Cell type identity markers (CD3, CD19, CD14)
- Immediate-early genes (activation state)

**ATAC excels at:**
- Regulatory landscape (enhancer/promoter accessibility)
- TF motif availability (chromatin priming)
- Lineage-defining epigenetic memory (developmental history)

**Example biological scenario:**
- Naive T cell CD8A: ATAC open, RNA low → priming, ready for TCR activation
- Effector T cell CD8A: ATAC open, RNA high → active transcription
- B cell CD8A: ATAC closed, RNA absent → silent, wrong lineage
- Memory T cell CD8A: ATAC half-open, RNA medium → intermediate state

WNN separates these four states at per-cell resolution.

## Dimensionality choice

**RNA PCA dims:**
- Standard: 1:30 (top 30 PCs)
- Alternative: 1:15 (aggressive), 1:50 (permissive)
- Guide: ElbowPlot inflection + biological requirement (more dims to detect rare cell types)

**ATAC LSI dims:**
- Standard: 2:30 (drop LSI-1)
- Multiome: 2:50 more common (extra dims add info)
- Guide: DepthCor plot (does LSI-2 still correlate with depth?) + biological need

**WNN dims:**
```r
FindMultiModalNeighbors(pbmc,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:30, 2:50))
```

Separate dim choice per modality. RNA is the first element (dims 1:30), ATAC the second (dims 2:50).

**Critical:** RNA first, ATAC second. `reduction.list` order matches `dims.list` order. Swapping them selects wrong dimensions per modality.

## Common pitfalls

1. **Not managing `DefaultAssay`.** RNA needs `DefaultAssay(pbmc) <- "RNA"` (or `"SCT"` after SCTransform); ATAC needs `<- "ATAC"`. Missing this makes FeaturePlot return empty plots.
2. **`seurat_v3` HVG order.** `sc.pp.highly_variable_genes(flavor='seurat_v3')` before `normalize_total`. Wrong order raises "expects raw counts" warning plus a potential loess failure.
3. **Forgetting to drop LSI-1.** Python muon does not auto-handle; manual `[:, 1:]` slice required.
4. **PCA versus LSI confusion.** RNA needs PCA, ATAC needs LSI. Confusing them destroys the analysis.
5. **Reduction dim count.** Both modalities' reductions must exist before WNN (`pca` and `lsi` names). Missing one raises "reduction not found" in `FindMultiModalNeighbors`.

## References

- Hafemeister and Satija 2019, *Genome Biol* — SCTransform paper
- Cusanovich et al. 2018, *Cell* — scATAC LSI methodology
- Stuart et al. 2019, *Cell* — Seurat integration (SCTransform + PCA)
- Bredikhin et al. 2022, *Genome Biol* — muon multi-modal processing docs
