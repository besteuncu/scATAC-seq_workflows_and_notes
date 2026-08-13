# 06. R (Signac) vs Python (SnapATAC2) — side-by-side comparison

## Motivation for parallel implementation

This repository processes the same 10x PBMC 10k ATACv2 dataset through two independent stacks:
- **R:** Seurat + Signac (Stuart lab, Nat Methods 2021)
- **Python:** scanpy + SnapATAC2 (Zhang, Nat Methods 2024)

Reasons to run both:
1. **Reproducibility validation.** If both frameworks converge on the same biological result, the analysis is robust to framework choice.
2. **Portfolio breadth.** A modern bioinformatician benefits from familiarity with both R and Python single-cell ecosystems.
3. **Framework fitness for purpose.** Each tool has strengths in specific contexts.
4. **Interview preparation.** "Why R and not Python?" becomes answerable with concrete experience.

## Framework philosophy contrast

**Signac (R):**
- Peak-first: trusts Cellranger's called peaks and builds on them
- Memory-first: object lives in RAM (Seurat convention)
- Feature-rich: peak-to-gene linkage, ChromVAR, footprinting built in
- Vignette-heavy: Stuart lab publishes vignettes for most use cases
- Community: large (Seurat ecosystem, millions of users)

**SnapATAC2 (Python):**
- Fragment-first: begins from raw fragments, calls own peaks or tiles the genome
- Disk-backed: object stored in .h5ad on disk, read on demand (atlas-scale friendly)
- Rust backend: substantial speed advantage on large datasets
- scverse-integrated: anndata, scanpy, scvi-tools, muon ecosystem
- Community: smaller but growing rapidly

## Step-by-step comparison

### 1. Data loading

**Signac:**
```r
counts <- Read10X_h5(counts_path)      # peak matrix
metadata <- read.csv(metadata_path)     # metadata
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  fragments = frag_path,                # link, not load
  genome = "hg38",
  min.cells = 10
)
pbmc <- CreateSeuratObject(counts = chrom_assay, assay = "peaks")
```

**SnapATAC2:**
```python
data = snap.pp.import_fragments(
  fragment_file=str(frag_path),          # read fragments directly
  chrom_sizes=snap.genome.hg38,
  file="pbmc.h5ad",                      # backed mode
  min_num_fragments=200,
)
```

**Differences:**
- Signac: peak matrix plus fragments (two files)
- SnapATAC2: fragments only (peaks called later)
- Signac: in-memory; SnapATAC2: disk-backed

### 2. QC metrics

**Signac:**
```r
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc, fast = FALSE)
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments /
                            pbmc$passed_filters * 100
```

**SnapATAC2:**
```python
snap.metrics.tsse(data, snap.genome.hg38)
snap.metrics.frag_size_distr(data)
```

**Differences:**
- Signac provides four metrics ready (TSS, nucleosome, FRiP, blacklist)
- SnapATAC2 has two metrics built in (TSS, fragment size); FRiP is computed later after peak calling
- **Scale mismatch:** Signac TSS is 2–10; SnapATAC2 tsse is 5–30. Same concept, different formulas. Do not copy thresholds.

### 3. Filtering

**Signac:**
```r
pbmc <- subset(
  x = pbmc,
  subset = nCount_peaks > 3000 &
           nCount_peaks < 100000 &
           pct_reads_in_peaks > 40 &
           nucleosome_signal < 4 &
           TSS.enrichment > 2
)
```

**SnapATAC2:**
```python
snap.pp.filter_cells(
  data,
  min_counts=3000,
  min_tsse=2,
  max_counts=100000,
)
```

**Differences:**
- Signac uses 5 thresholds (explicit, verbose)
- SnapATAC2 uses 3 thresholds (minimal, opinionated)
- SnapATAC2 deprioritizes nucleosome and FRiP filters

### 4. Normalization

**Signac:**
```r
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q5")
```

**SnapATAC2:**
```python
snap.pp.add_tile_matrix(data, bin_size=500)
snap.pp.select_features(data, n_features=25000)
```

**Differences:**
- Signac: TF-IDF on peak matrix (~90K peaks)
- SnapATAC2: spectral embedding on tile matrix (~6M 500 bp bins); different feature space

### 5. Dimensionality reduction

**Signac:**
```r
pbmc <- RunSVD(pbmc)
DepthCor(pbmc)                          # inspect LSI-1 depth confound
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 2:30)   # drop LSI-1
```

**SnapATAC2:**
```python
snap.tl.spectral(data)                  # handles depth confound internally
```

**Differences:**
- Signac: LSI plus manual `dims = 2:30`
- SnapATAC2: spectral, automatic depth handling

### 6. Clustering

**Signac:**
```r
pbmc <- FindNeighbors(pbmc, reduction = "lsi", dims = 2:30)
pbmc <- FindClusters(pbmc, algorithm = 3, resolution = 0.8)
```

**SnapATAC2:**
```python
snap.pp.knn(data)
snap.tl.leiden(data, resolution=0.8)
```

**Differences:** None substantive. Same Leiden algorithm, same KNN approach; only syntax differs.

### 7. UMAP

**Signac:**
```r
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 2:30, seed.use = 42)
DimPlot(pbmc, label = TRUE)
```

**SnapATAC2:**
```python
snap.tl.umap(data)
snap.pl.umap(data, color="leiden")
```

**Differences:** UMAP shapes differ due to seed and implementation, but neighborhood structure is preserved.

### 8. Gene activity and annotation

**Signac:**
```r
gene.activities <- GeneActivity(pbmc)
pbmc[["RNA"]] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(pbmc, assay = "RNA")
FeaturePlot(pbmc, features = c("MS4A1", "CD3E"))
```

**SnapATAC2 + scanpy:**
```python
gene_matrix = snap.pp.make_gene_matrix(data, gene_anno=snap.genome.hg38, inplace=False)
sc.pp.normalize_total(gene_matrix)
sc.pp.log1p(gene_matrix)
sc.pl.umap(gene_matrix, color=["MS4A1", "CD3E"])
```

**Differences:**
- Signac auto-normalizes via `NormalizeData`
- SnapATAC2 requires explicit `sc.pp.normalize_total + log1p` (easy to forget; results plot as mostly-dark if skipped)

## Results comparison on this repository's data

Same 10x PBMC 10k ATACv2 dataset through both stacks:

| Metric | Signac (R) | SnapATAC2 (Python) |
|--------|-----------|--------------------|
| Cells after filtering | 9,850 | 10,149 |
| Feature space | 165,434 peaks | 6,062,095 tiles (500bp) → 25,000 selected |
| Cluster count (resolution 0.8) | 18 | 13–15 |
| Major cell types identified | T, B, monocyte, NK — separated | T, B, monocyte, NK — separated |
| Canonical marker overlay | Clean separation in cluster-marker plots | Same |
| Runtime (setup + pipeline) | ~15 min | ~25 min (includes I/O overhead) |
| Peak memory | ~4 GB RAM | ~1.5 GB RAM (disk-backed) |

**Cell count difference (3%):** Signac's filter is slightly stricter (nucleosome and FRiP included). SnapATAC2's is looser. Biology is the same.

**Cluster count difference (18 vs 13):** Different feature spaces (peak vs tile) and different algorithms. Major cell types match; sub-clustering differs.

**UMAP shape:** Different visual layouts. But cluster neighborhood structure (which cluster is near which) is similar.

## When to use which

| Use case | Recommendation | Reason |
|----------|---------------|--------|
| Small dataset (<50K cells) | Signac | Vignette-heavy, fast prototyping |
| Atlas scale (1M+ cells) | SnapATAC2 | Rust backend, disk-backed |
| Peak-based specific analysis (footprinting, ChromVAR, LinkPeaks) | Signac | Built-in comprehensive |
| Python ML/DL integration (scVI, deep learning) | SnapATAC2 | scverse ecosystem |
| Multiome (WNN integration) | Signac | Seurat WNN mature |
| Multiome (multiVI VAE) | SnapATAC2 + scvi-tools | Native VAE support |
| CV: R-heavy lab | Signac | Ergonomic fit |
| CV: Python-heavy lab | SnapATAC2 | Ergonomic fit |

## Portfolio narrative

In a CV or interview:

"In scATAC-seq analysis I processed the same 10x PBMC 10k dataset through two independent stacks: Signac + Seurat in R and SnapATAC2 + scanpy in Python. Both frameworks recovered the same major PBMC populations (T subtypes, B, monocytes, NK) with similar cluster counts and marker patterns. The 3% cell-count and 20% cluster-count differences arose from framework-specific filter and algorithm choices; the biological findings were robust. Framework choice depends on ecosystem, scale, and downstream analysis needs."

This conveys three things:
1. Technical breadth (two frameworks)
2. Scientific rigor (reproducibility validation)
3. Judgment (understanding framework trade-offs)

## Migration paths

**Signac to SnapATAC2:**
- Fragments file is reused via `snap.pp.import_fragments`
- Convert Seurat object to SingleCellExperiment with `zellkonverter`, then to AnnData
- Peak matrix is not converted to tile matrix; instead proceed in the SnapATAC2 tile-space

**SnapATAC2 to Signac:**
- Start from the same fragments file
- Generate peak matrix (SnapATAC2's own peak calling or external MACS2)
- Create Signac object with `CreateChromatinAssay + CreateSeuratObject`

## Common framework-specific pitfalls

**Signac:**
1. Forgetting to drop LSI-1 in `dims = 2:30`
2. Wrong `DefaultAssay` set, causing FeaturePlot to return blank
3. Invalid fragment path causing TSSEnrichment to hang
4. Deprecated `NucleosomeSignal` and `TSSEnrichment`; new API is `ATACqc()`

**SnapATAC2:**
1. `data.obs` is a polars DataFrame, not pandas — use `[:].columns` instead of `.columns`
2. `.h5ad` file lock issue in backed mode — cannot open the same file from two processes
3. `snap.genome.hg19` triggers UCSC download; use dict fallback if network is slow
4. Distinction between `filtered_feature_bc_matrix.h5` and `raw_..._matrix.h5`

## References

- Stuart et al. 2021, *Nat Methods* — Signac
- Zhang et al. 2024, *Nat Methods* — SnapATAC2
- Bredikhin et al. 2022, *Genome Biol* — muon (multi-modal Python)
- Wolf et al. 2018, *Genome Biol* — scanpy (Python single-cell base)
- Hao et al. 2021, *Cell* — Seurat v4 (WNN, R base framework)
