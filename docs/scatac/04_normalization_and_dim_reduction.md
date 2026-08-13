# 04. Normalization and dimensionality reduction — TF-IDF and LSI

## Goal

Turn a sparse binary-like cell × peak matrix (roughly 10K cells × 100K peaks, 95%+ zeros) into a meaningful 30-dimensional embedding. The embedding is the substrate for downstream clustering and UMAP.

## Why classic scRNA methods fail here

The standard scRNA-seq pipeline is:
1. `NormalizeData` (log-normalize)
2. `FindVariableFeatures`
3. `ScaleData`
4. `RunPCA`

This does not translate to scATAC-seq for three reasons:
- **Count magnitude carries no information.** Values are 0/1/2; log-scaling adds no signal.
- **Extreme sparsity.** 95%+ zeros makes PCA's centering assumption catastrophic (dense representation blows up).
- **Feature count 10× larger.** ~100K peaks versus ~20K genes makes computation expensive and overfitting more likely.

A different normalization plus dimensionality reduction pair is required. TF-IDF plus LSI is the standard.

## TF-IDF, borrowed from NLP

**TF-IDF** = Term Frequency × Inverse Document Frequency.

Developed in the 1980s in information retrieval for sparse binary document × term matrices. The problem it solves is structurally identical to what scATAC-seq faces:

| NLP | scATAC |
|-----|--------|
| Document | Cell |
| Term (word) | Peak |
| Document length | Cell sequencing depth |
| Common words ("the", "and") | Ubiquitous peaks (promoters) |
| Informative rare words | Cell-type-specific enhancers |

### TF — Term Frequency

Formula:
```
TF[i, j] = count[i, j] / total_count[i]
```

The count of peak `j` in cell `i`, divided by the total count in cell `i`. This is **depth normalization**: a deep-sequenced cell's high counts are relativized.

### IDF — Inverse Document Frequency

Formula:
```
IDF[j] = log(N_cells / N_cells_with_peak[j])
```

Peak `j`'s inverse frequency across cells. Ubiquitous peaks get low IDF (small weight). Rare peaks get high IDF (large weight).

Example with 10,000 cells:
- A promoter peak accessible in 9,500 cells: IDF = log(10000/9500) ≈ 0.05 → very small weight
- A CD8 T cell enhancer peak accessible in 500 cells: IDF = log(10000/500) ≈ 3.0 → large weight

Biologically meaningful: informative peaks (cell-type discriminators) get up-weighted; universally accessible peaks get muted.

### Combined form

Signac uses:
```
TF-IDF[i, j] = log(1 + TF[i, j] × IDF[j] × 10000)
```

The `log(1 + x × 10000)` form stabilizes variance (mirrors scRNA's `NormalizeData`). Values fall in a reasonable numeric range for SVD.

Function in Signac: `RunTFIDF(pbmc)`.
Function in SnapATAC2: implicit inside `snap.pp.select_features`.

## Feature selection

Rather than running SVD on all 100K peaks, first select the top variable features. Two benefits:
- Faster computation (10K × 30K instead of 10K × 100K).
- Less noise (low-variance peaks are not discriminatory).

Signac: `FindTopFeatures(pbmc, min.cutoff = "q5")` retains the top 95% by count.
SnapATAC2: `snap.pp.select_features(data, n_features=25000)` retains the top 25K.

Setting `min.cutoff = "q0"` keeps all peaks (conservative). Setting `"q75"` retains only the top 25% by variance (aggressive, may drop rare-cell-type peaks).

## LSI — Latent Semantic Indexing

LSI is simply SVD applied to the TF-IDF matrix.

### What SVD does

Decompose an M × N matrix:
```
X = U × Σ × V^T
```

Where `U` is the left singular vectors, `Σ` is the diagonal of singular values, and `V^T` is the right singular vectors. Choosing 30 leading vectors approximates the original 100K-dimension space with only 30 dimensions.

Biologically, the first 30 singular vectors are the 30 largest variance directions in the data. Biological structure lives in these directions.

### PCA versus LSI

Same core idea, different framing:
- **PCA**: SVD on centered matrix (`X - mean(X)`).
- **LSI**: SVD on TF-IDF matrix, no centering.

For sparse binary matrices, centering densifies the storage catastrophically. LSI's uncentered form is what makes it feasible.

The name "LSI" is preserved from NLP terminology. Mathematically it is SVD.

## The LSI-1 depth confound

The most important detail in this pipeline.

**The problem:** The largest singular vector (LSI-1) typically correlates strongly with sequencing depth. Its correlation value is 0.7 to 0.9. LSI-1 separates cells by "high depth versus low depth", not by biology.

**Why:** TF-IDF corrects depth partially but not fully. The number of peaks a cell has non-zero values in (feature richness) still scales with depth, and this residual variance dominates the first singular vector.

**Visual verification:**
```r
DepthCor(pbmc)   # Signac
```

Expected output: the LSI-1 bar exceeds 0.7 in absolute value; LSI-2, LSI-3 stay below 0.3.

**Fix:**
- Signac: pass `dims = 2:30` to every downstream call (UMAP, FindNeighbors).
- SnapATAC2: the spectral algorithm handles depth internally, so no manual drop is needed. Optionally slice `data.obsm['X_lsi'][:, 1:]`.

## Choosing the number of dimensions

Convention: 30 dimensions (with LSI-1 dropped, effectively 29).

Alternatives:
- `dims = 2:30` — safe default
- `dims = 2:50` — captures more variance, but includes noise
- `dims = 2:15` — too few, biological structure lost

Elbow plot for guidance:
```r
ElbowPlot(pbmc, ndims = 50, reduction = "lsi")
```

The elbow typically appears at 15–25 dimensions. Beyond that is noise. In practice, 30 dimensions works for most datasets.

## Sanity checks after dim reduction

1. **DepthCor plot** — LSI-1 correlation > 0.7 (as expected); LSI-2+ correlations well below 0.4.
2. **LSI embedding shape** — `dim(pbmc@reductions$lsi@cell.embeddings)` should be `N_cells × 30`.
3. **UMAP on `dims = 2:30`** — should show cluster structure, not a single blob.
4. **Depth overlay UMAP** — `FeaturePlot(pbmc, features = "nCount_peaks")` on UMAP should NOT show a strong gradient. A gradient means LSI-1 still leaks into downstream dims.

## SnapATAC2 spectral embedding

SnapATAC2's `snap.tl.spectral` is a modernized variant of LSI:
- Depth confound handled automatically (no manual dim-1 drop)
- Graph-based spectral decomposition
- Results closely track LSI-based methods, roughly 95% concordance downstream

Pipeline:
```python
snap.pp.select_features(data, n_features=25000)
snap.tl.spectral(data)
```

Result stored in `data.obsm['X_spectral']`.

## Common pitfalls

1. **Forgot to drop LSI-1.** Passing `dims = 1:30` instead of `dims = 2:30` is the most common error.
2. **Skipped depth overlay check.** Even after dropping LSI-1, other dimensions may correlate with depth. Visual inspection is required.
3. **Aggressive `min.cutoff`.** Setting `"q75"` on a small dataset drops rare-cell-type peaks and their signals with them.
4. **No seed for UMAP.** Without `seed.use = 42`, layouts change run to run.
5. **PCA instead of SVD.** `RunPCA` is inappropriate for scATAC; use `RunSVD`.

## Why this stage is the "money block"

Everything downstream operates on the LSI embedding:
- Clustering uses it
- UMAP uses it
- Trajectory analysis uses it
- Batch correction uses it

If the embedding is compromised, downstream analysis is compromised. This stage is the engineered feature representation that scATAC needs to work at all. The "aha" moment of scATAC — cells grouping by type — happens here.

## References

- Deerwester et al. 1990 — original LSI paper (NLP context)
- Cusanovich et al. 2018, *Cell* — first LSI on scATAC
- Stuart et al. 2021, *Nat Methods* — Signac TF-IDF and LSI implementation
- Zhang et al. 2024, *Nat Methods* — SnapATAC2 spectral embedding
