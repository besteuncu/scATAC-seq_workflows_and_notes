# 03. Multiome QC and intersection filtering

## Why Multiome QC differs from scATAC-only

The scATAC-only pipeline used four metrics. Multiome requires independent metric groups for each modality:

**RNA metrics:**

| Metric | Measures | Threshold |
|--------|----------|-----------|
| `nCount_RNA` | Total mRNA count | 1,000 – 25,000 |
| `nFeature_RNA` | Number of unique genes | 200 – 5,000 |
| `percent.mt` | Mitochondrial gene fraction | < 20% |

**ATAC metrics:**

| Metric | Measures | Threshold |
|--------|----------|-----------|
| `nCount_ATAC` | Total ATAC fragments | 1,000 – 100,000 |
| `TSS.enrichment` | TSS-proximal signal density | > 1 (Signac scale) |
| `nucleosome_signal` | Fragment length nucleosome pattern | < 2 |

## Why intersection, not union

**Intersection filter:** A cell is retained only if it passes both RNA and ATAC thresholds. **Union** would retain a cell passing either.

**Multiome requires intersection.** Consider two scenarios:

*Scenario A — RNA good, ATAC bad (TSS = 0.5, nucleosome = 3.5):*
- RNA analysis returns reasonable results
- ATAC pipeline injects noise into this cell's peak values
- WNN integration weights the cell's ATAC contribution near zero, but the cell still contributes "dirty" information to the joint embedding
- **Drop this cell.** It does not represent proper biology in the Multiome context.

*Scenario B — RNA bad (mito 35%), ATAC good:*
- Dead-cell signature in RNA
- ATAC is likely ambient or artifact
- **Drop this cell.**

Intersection filter is a Multiome-defining feature. Union retains ~20% more cells but degrades cluster quality.

## Selecting RNA versus ATAC thresholds

Empirical values for 10x Multiome PBMC 10k:

| Metric | Lower | Upper | Reason |
|--------|-------|-------|--------|
| RNA count | 1000 | 25000 | Lower: shallow droplet; upper: doublet candidate |
| RNA feature | 200 | 5000 | Lower: not a cell; upper: doublet |
| Mito % | — | 20% | Upper: dead or stressed cell |
| ATAC count | 1000 | 100000 | Lower: shallow; upper: doublet |
| TSS enrichment | 1 | — | Lower: ambient DNA |
| Nucleosome signal | — | 2 | Upper: overtagmented |

**Calibrating for your dataset:** Inspect violin plots and histograms. Retain the distribution's "body" and drop the tails. Thresholds are dataset-specific.

## Cellranger-ARC metadata differences

The `per_barcode_metrics.csv` from Cellranger-ARC (Multiome) differs from `singlecell.csv` (scATAC-only). Multiome-specific extras include:
- `is_cell` — Cellranger's own cell call
- `excluded_reason` — filter decision reason
- `gex_umis_count` — RNA UMI total (Multiome context)
- `atac_fragments` — ATAC fragment total
- `atac_tss_fragments` — TSS-region fragment count
- `linked_atac_barcode` / `linked_gex_barcode` — cross-modality barcode consistency (should always match)

These are available for custom QC calculations but are not required by the standard Signac/muon workflow.

## Fragment length histogram

`FragmentHistogram(pbmc, group.by = "nucleosome_group")` should show two peaks for healthy ATAC:
- ~50 bp (nucleosome-free)
- ~200 bp (mono-nucleosomal)

Multiome protocols split libraries between RNA and ATAC, so per-cell ATAC coverage is slightly lower than scATAC-only. The histogram shape should remain bimodal. Loss of bimodality (single peak or uniform) indicates a preparation problem.

## Doublet detection — Multiome's easier path

In scATAC-only, doublet detection was difficult (ArchR required). In Multiome, RNA can be inspected with scRNA doublet tools:

**R side:**
```r
library(scDblFinder)
sce <- as.SingleCellExperiment(DietSeurat(pbmc, assays = "RNA"))
sce <- scDblFinder(sce)
pbmc$doublet_score <- sce$scDblFinder.score
pbmc$is_doublet <- sce$scDblFinder.class == "doublet"

# Add to filter
pbmc <- subset(pbmc, subset = is_doublet == FALSE & ...)
```

**Python side:**
```python
import scrublet
scrub = scrublet.Scrublet(rna.X)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
rna.obs['doublet_score'] = doublet_scores
rna.obs['is_doublet'] = predicted_doublets
```

Typical doublet fraction for PBMC 10k: 5–8%. Removing doublets improves WNN quality substantively. This pipeline does not include doublet filtering by design (simplicity); future work.

## Interpreting the QC violin plot

Standard 6-panel display:
```r
VlnPlot(pbmc, features = c("nCount_RNA", "nFeature_RNA", "percent.mt",
                            "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
        ncol = 3, log = TRUE)
```

Look for:
1. **RNA count and feature distribution:** unimodal is ideal; right-tailed with long upper tail = doublet candidates.
2. **Percent mitochondrial:** left peak (<20%) is the main population; the right tail = dead cells.
3. **ATAC count:** right-skewed; extreme upper values indicate doublets.
4. **TSS enrichment:** peak around 5–10; left tail (<2) is below the filter threshold.
5. **Nucleosome signal:** main peak around 0.5–1; tail above 2 fails the filter.

Red flag: bimodal distribution. Either two populations (batch, dying vs healthy) or the filter's cut-off point falls in a valid region.

## Post-filter verification

`summary(pbmc@meta.data[, ...])` should confirm:
- Minimum values above threshold (filter held)
- Maximum values below threshold (upper limit held)
- Median values matching a high-quality dataset

**Expected for PBMC 10k Multiome after filter:**
- ~10,500 cells (from 11,900)
- RNA count median ~4,000
- RNA feature median ~1,900
- Mito % median ~10 (limit 20)
- ATAC count median ~23,000
- TSS.enrichment median ~5.5
- Nucleosome signal median ~0.67

These values reflect a high-quality Multiome dataset.

## Common pitfalls

1. **Union filter as accidental default.** Use `A AND B`, not `A OR B`. Independent modality evaluation requires intersection.
2. **Forgetting `mu.pp.intersect_obs`.** After per-modality filtering in Python, this call synchronizes cell names. Without it, WNN throws a cell-count mismatch error.
3. **Metadata column names differ from scATAC.** Cellranger-ARC (Multiome) columns differ from Cellranger-ATAC (`atac_fragments` vs `passed_filters`). Verify with `colnames(pbmc@meta.data)`.
4. **Wrong mito prefix.** Human uses `MT-`; mouse uses `mt-`. Wrong prefix returns `percent.mt = 0` for all cells (fake healthy).
5. **Skipping doublet detection.** Multiome doublet fraction is 5–10%. Optional but improves analysis quality.

## References

- 10x Multiome QC recommendations: https://kb.10xgenomics.com/hc/en-us/articles/360061165691
- Signac Multiome QC vignette: https://stuartlab.org/signac/articles/pbmc_multiomic
- muon QC docs: https://muon.readthedocs.io/en/latest/omics/atac.html
