# 02. Multiome pipeline overview — 8 blocks

## General flow

The scATAC-only 6-block pipeline is extended with WNN integration to become 8 blocks:

```
Load (RNA + ATAC) → QC per modality → Filter intersection
                                             |
Process RNA (SCT + PCA) ← ← ← ← ← ← ← Process ATAC (TF-IDF + LSI)
       |                                     |
       └────────── WNN integration ──────────┘
                          |
                Joint UMAP + clustering
                          |
                 Annotation (RNA markers)
```

## The 8 blocks: biology, math, meaning

| Block | Biological question | Mathematical operation | What it tells us |
|-------|--------------------|-----------------------|------------------|
| 1. Setup | Is the environment ready? | Library load + path check | Reproducibility |
| 2. Load | How do the two modalities join? | H5 parse → MuData / SeuratObject | Paired cell × RNA + cell × ATAC |
| 3. QC | Which cells are healthy in both modalities? | Per-modality metrics (TSS, nucleosome, mito%) | Quality score per cell |
| 4. Filter | Intersection filter (RNA AND ATAC good) | Threshold AND | Analysis-ready cell subset |
| 5. Process RNA | Which genes are discriminative in RNA? | SCTransform + PCA | RNA embedding (30D) |
| 6. Process ATAC | Which peaks are discriminative in ATAC? | TF-IDF + LSI dims 2:50 | ATAC embedding (48D) |
| 7. WNN integration | How to combine the two modalities? | Weighted Nearest Neighbor | Joint neighbor graph |
| 8. Cluster + Annotate | What cell type per cluster? | Leiden on joint graph + reference | Biological identity |

## Block details

### Block 1: setup

**Purpose:** Load libraries, define paths, verify files.

**Multiome-specific:** R also loads `BSgenome.Hsapiens.UCSC.hg38` (for downstream motif work). Python adds `muon` and `mudata`.

**Duration:** ~5 seconds.

### Block 2: load Multiome data

**Purpose:** Parse `filtered_feature_bc_matrix.h5` into two assays and combine into a multi-modal object.

**Mathematics:** The H5 contains two matrices — "Gene Expression" (~36K genes × ~11K cells) and "Peaks" (~140K peaks × ~11K cells). Signac attaches them as separate assays in a SeuratObject. muon wraps them as two AnnDatas inside a MuData.

**Critical:** Both assays share the same cell barcode set. This symmetry is a guarantee for all downstream steps — `pbmc$RNA` and `pbmc$ATAC` will always contain the same ~11K cells. Fragment file is linked (not loaded) to the ATAC assay for memory efficiency.

**Duration:** ~30 seconds.

### Block 3: QC per modality

**Purpose:** Compute per-cell QC metrics independently for RNA and ATAC.

**RNA metrics:**
- `nCount_RNA` — total mRNA counts
- `nFeature_RNA` — number of unique genes detected
- `percent.mt` — mitochondrial gene fraction (nucleus prep quality)

**ATAC metrics:**
- `nCount_ATAC` — total ATAC fragments (peak-level counts)
- `TSS.enrichment` — Signac scale, cell integrity
- `nucleosome_signal` — mono/free ratio

Each modality has its own failure mode. Mito > 20% in RNA suggests a dying cell; TSS < 2 in ATAC suggests ambient DNA. Thresholds are independent.

**Duration:** RNA is instant, ATAC TSS enrichment takes 3–8 minutes.

### Block 4: filter cells (intersection)

**Purpose:** Retain only cells passing both RNA and ATAC thresholds.

**Threshold example:**
```r
subset = nCount_ATAC > 1000 & nCount_ATAC < 100000 &
         nCount_RNA > 1000 & nCount_RNA < 25000 &
         nucleosome_signal < 2 & TSS.enrichment > 1 &
         percent.mt < 20
```

**Why intersection:** A cell strong in RNA but weak in ATAC will contribute noise to WNN's ATAC-side signal. Strict filtering improves cluster quality.

**Expected outcome:** 85–95% of cells retained (~11K → ~10.5K).

**Python note:** `mu.pp.intersect_obs(mdata)` is required after independent modality filtering; otherwise WNN throws mismatched cell errors.

### Block 5: process RNA (SCTransform + PCA)

**Purpose:** Normalize RNA, select variable features, run PCA → RNA-only cell embedding.

**Mathematics (R):**
- `SCTransform` — Poisson-based variance stabilization; combines depth normalization, variance stabilization, and top variable feature selection
- `RunPCA` — principal component analysis on 50 PCs

**Mathematics (Python):**
- `highly_variable_genes(flavor='seurat_v3')` on raw counts (order-critical)
- `normalize_total + log1p` for depth normalization
- `scale` for mean 0, variance 1
- `pca` for 50 PCs

**Why PCA works for RNA (unlike LSI for ATAC):** RNA counts are less sparse (70–80% zero vs 95%+ in ATAC), and count magnitude is meaningful. PCA's centering does not fail on RNA data.

**Critical:** In Multiome, RNA is the "gold standard annotation" source. WNN's integration quality depends on RNA embedding quality.

**Duration:** ~5 minutes. Elbow plot shows an inflection at PC 20–30.

### Block 6: process ATAC (TF-IDF + LSI)

**Purpose:** Normalize ATAC with sparse-appropriate method, then reduce dimensions.

**Same as scATAC-only pipeline:**
- `RunTFIDF` — NLP-style normalization
- `FindTopFeatures(min.cutoff="q0")` — top variable peaks
- `RunSVD` — LSI (SVD on TF-IDF), 50 dims

**Multiome-specific:** WNN uses `dims = 2:50` (drop LSI-1 depth confound). This is more dims than scATAC-only's typical `2:30`; the extra dims add information that WNN can use.

**Python:** muon's LSI does not automatically handle depth confound. Manually slice:
```python
atac.obsm['X_lsi'] = atac.obsm['X_lsi'][:, 1:]
```

**Duration:** ~2–3 minutes (TF-IDF + SVD).

### Block 7: WNN integration

**Purpose:** Combine the two modality embeddings with adaptive per-cell weighting. WNN determines which modality is more informative for each cell.

**Mathematics:**
1. Per-modality KNN (K=20) within each cell's embedding
2. Per-cell modality prediction accuracy → modality weight per cell
3. Weighted KNN → joint neighbor graph
4. Result stored in `weighted.snn` (R) or `wnn_connectivities` (Python)

**Why weighted:** Modalities have different signal-to-noise per cell. Simple concatenation ignores this; WNN adapts.

**Signac API:**
```r
pbmc <- FindMultiModalNeighbors(
  pbmc,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:30, 2:50)
)
```

**muon API:**
```python
mu.pp.neighbors(mdata, key_added='wnn')
```

**Duration:** 2–4 minutes.

### Block 8: cluster + joint UMAP + annotation

**Purpose:** Leiden clustering on WNN joint graph + UMAP visualization + reference-based automated annotation.

**Cluster (R):**
```r
pbmc <- FindClusters(pbmc, graph.name = "wsnn",
                      algorithm = 3, resolution = 0.8)
```

**Cluster (Python):**
```python
sc.tl.leiden(mdata, obsp='wnn_connectivities',
             resolution=0.8, key_added='leiden')
```

**Note:** `mu.tl.leiden(neighbors_key='wnn')` fails because muon expects per-modality graphs. Use `sc.tl.leiden` with `obsp` argument for joint graph clustering.

**Automated annotation:**
- R: Azimuth PBMC reference (`RunAzimuth(pbmc, reference = "pbmcref")`)
- Python: CellTypist `Immune_All_Low.pkl` model

Both provide labels plus confidence scores, replacing manual marker-based annotation.

**Duration:** 1–2 minutes for clustering; 2–5 minutes for Azimuth; ~30 seconds for CellTypist.

## Expected results (10x PBMC Multiome 10k)

- Cell count after filtering: ~10.5K (from 11.9K)
- Cluster count (resolution 0.8): 15–20
- Major cell type composition:
  - CD4 T (~30%): naive, memory, effector subsets
  - CD8 T (~20%): naive, memory, effector, MAIT/NKT
  - Monocytes (~20%): CD14 classical, CD16 non-classical
  - B cells (~10%): naive, memory
  - NK cells (~10%): CD16, CD56 bright/dim
  - Rare (<5%): DC, pDC, HSPC

## R versus Python side by side

| Block | Seurat + Signac (R) | muon + scanpy (Python) |
|-------|--------------------|------------------------|
| Load | `Read10X_h5` → SeuratObject | `mu.read_10x_h5` → MuData |
| QC RNA | `PercentageFeatureSet` | `sc.pp.calculate_qc_metrics(qc_vars=['mt'])` |
| QC ATAC | `NucleosomeSignal + TSSEnrichment` | `ac.tl.nucleosome_signal + ac.tl.tss_enrichment` |
| Filter | `subset(..., subset = ...)` | `mu.pp.filter_obs + mu.pp.intersect_obs` |
| RNA process | `SCTransform + RunPCA` | `sc.pp.hvg + normalize + scale + sc.tl.pca` |
| ATAC process | `RunTFIDF + FindTopFeatures + RunSVD` | `ac.pp.tfidf + ac.tl.lsi` |
| Integration | `FindMultiModalNeighbors` (WNN) | `mu.pp.neighbors` (WNN equivalent) |
| Cluster | `FindClusters(graph.name="wsnn")` | `sc.tl.leiden(obsp='wnn_connectivities')` |
| UMAP | `RunUMAP(nn.name="weighted.nn")` | `mu.tl.umap(neighbors_key='wnn')` |
| Annotate | Azimuth PBMC reference | CellTypist Immune_All_Low |

## Beyond Block 8

The basic pipeline enables advanced downstream:
- **Peak-to-gene linkage:** `LinkPeaks(pbmc, peak.assay="ATAC", expression.assay="SCT")` — Multiome's main analytical strength
- **Differential accessibility per cluster:** pseudobulk + DESeq2
- **Motif enrichment in cluster-marker peaks:** ChromVAR
- **Cell trajectory (developmental):** Monocle3 or scFates

These are follow-up sprint topics building on the basic pipeline.
