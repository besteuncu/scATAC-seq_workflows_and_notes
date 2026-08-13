# 03. Quality control — metrics, thresholds, decisions

## Purpose

Cellranger-ATAC returns roughly 10,000 candidate cell barcodes after its own cell-calling filter. This is a coarse pass; the filtered output still contains:
- Ambient DNA (empty droplets and DNA fragments)
- Dead or damaged cells (broken nuclei)
- Doublets (two cells sharing a barcode)
- Overtagmented or degraded samples

QC removes these to prevent noise contamination downstream.

## Four key metrics

### 1. TSS enrichment score

**Measures:** Whether a cell's fragments are enriched at transcription start sites.

**Formula:** Tn5 cut density in the TSS ±100 bp window divided by flanking background density (±2 kb). Higher ratio, better cell.

**Biological interpretation:** Healthy nuclei show open chromatin concentrated at TSS regions of active promoters. High TSS enrichment indicates fragments come from biologically meaningful locations. Low TSS enrichment indicates fragments from random genomic regions, characteristic of ambient DNA or broken nuclei.

**Threshold conventions:**
- `> 2` — Signac vignette minimum, widely accepted
- `> 4` — good quality
- `> 7` — ENCODE excellent standard
- `< 2` — likely not a cell; filter out

**Signac function:** `TSSEnrichment(pbmc, fast = FALSE)`
**SnapATAC2 function:** `snap.metrics.tsse(data, snap.genome.hg38)`

Note that the two frameworks compute TSS enrichment on different scales. Signac's `TSS.enrichment` typically ranges 2–10; SnapATAC2's `tsse` ranges 5–30. Same biological property, different formulas. Do not copy thresholds across frameworks without adjustment.

### 2. Nucleosome signal

**Measures:** Whether the fragment length distribution captures the nucleosome pattern.

**Formula:** Ratio of mono-nucleosomal fragments (147–294 bp) to nucleosome-free fragments (<147 bp).

**Biological interpretation:** Healthy Tn5 tagmentation produces a bimodal distribution:
- ~50 bp peak (nucleosome-free, Tn5 cuts within a nucleosome-free stretch)
- ~200 bp peak (mono-nucleosomal, Tn5 cuts on either side of one nucleosome)

The bimodality itself signals clean tagmentation.

**Threshold:**
- `< 4` — acceptable
- `> 4` — overtagmented or degraded; filter out

**Signac function:** `NucleosomeSignal(pbmc)` (deprecated; new API `ATACqc()`).

Fragment length histogram (`FragmentHistogram(pbmc, group.by = "nucleosome_group")`) should show two peaks. A single peak or uniform distribution indicates a preparation problem.

### 3. Fraction of Reads in Peaks (FRiP)

**Measures:** What fraction of a cell's fragments fall inside peaks.

**Formula:** `pct_reads_in_peaks = peak_region_fragments / passed_filters * 100`.

**Biological interpretation:** High FRiP indicates fragments come from biologically meaningful regulatory regions. Low FRiP indicates noise, background, or ambient DNA.

**Threshold:**
- `> 40%` — good quality
- `> 60%` — excellent
- `< 30%` — likely not a cell or ambient-dominated

Source: Cellranger metadata (`singlecell.csv`) or computed post-hoc.

### 4. Total fragment count

**Measures:** Sequencing depth per cell.

**Biological interpretation:** Very low counts imply shallow droplets. Very high counts imply doublets (two cells' fragments summed).

**Threshold:**
- Lower bound: 3,000 (Signac vignette) or 1,000 (permissive)
- Upper bound: 100,000 (doublet upper limit)

ATACv2 chemistry has median ~28,000 fragments per cell; older ATACv1 had ~5,000. Modern chemistry is much deeper.

## Additional metric: blacklist ratio

**Measures:** Fraction of fragments that overlap ENCODE blacklist regions (repetitive and low-mappability).

**Biological interpretation:** High blacklist ratio indicates mapping artifacts.

**Threshold:** `< 0.05`.

Note: with ATACv2 chemistry and Cellranger-ATAC v2, this value is often 0 for all cells, because v2's upstream processing handles blacklist differently. In this case, drop the blacklist filter.

## Interpreting QC plots

Standard 4-panel violin plot:
```r
VlnPlot(pbmc, features = c("nCount_peaks", "TSS.enrichment",
                            "nucleosome_signal", "pct_reads_in_peaks"),
        pt.size = 0.1, ncol = 4)
```

Look for:
1. **Distribution shape** — long-tailed unimodal is ideal; bimodal suggests two populations, one likely of lower quality.
2. **Lower tail** — cells below threshold, will be filtered.
3. **Upper tail (fragment count)** — extreme values are doublet candidates.
4. **Outliers** — sparse points are doublets; clustered points may be a rare biological population.

Fragment length histogram: expect two peaks (nucleosome-free ~50 bp, mono-nucleosomal ~200 bp). A single peak or uniform distribution indicates a preparation problem.

## Applying the filter

Use intersection of multiple thresholds (AND, not OR):

```r
pbmc <- subset(
  x = pbmc,
  subset = nCount_peaks       > 3000 &
           nCount_peaks       < 100000 &
           pct_reads_in_peaks > 40 &
           nucleosome_signal  < 4 &
           TSS.enrichment     > 2
)
```

Typical outcome: 85–95% of Cellranger-approved cells pass, leaving 8,500–9,500 cells. ATACv2 datasets often lose fewer cells due to higher baseline quality.

The SnapATAC2 equivalent is more compact:
```python
snap.pp.filter_cells(
    data,
    min_counts=3000,
    min_tsse=2,
    max_counts=100000,
)
```

SnapATAC2 uses only three metrics by default. Nucleosome and FRiP filters are omitted by design — SnapATAC2 relies on downstream cluster-aware peak calling to clean things up.

## Post-filter verification

```r
summary(pbmc@meta.data[, c("nCount_peaks", "TSS.enrichment",
                            "nucleosome_signal", "pct_reads_in_peaks")])
```

Confirm minimum values exceed thresholds. Distribution should look cleaner overall.

## Doublet detection

The threshold filter partially catches doublets via high fragment count, but doublet detection is a separate concern.

Best practice: use ArchR's `addDoubletScores()` (LSI-based, scATAC-specific). Signac does not implement this directly; the workaround is to briefly pass through ArchR just for doublet scoring, then return to Signac.

Alternatives (scRNA-focused but adaptable): `scDblFinder` (R), `Scrublet` (Python).

This repository's pipeline does not include doublet detection for simplicity. ArchR integration is future work.

## Decision log

| Metric | Threshold | Reason |
|--------|-----------|--------|
| `nCount_peaks > 3000` | lower | remove shallow droplets |
| `nCount_peaks < 100000` | upper | remove doublet candidates |
| `TSS.enrichment > 2` | lower | cell integrity |
| `nucleosome_signal < 4` | upper | remove overtagmented cells |
| `pct_reads_in_peaks > 40` | lower | remove ambient DNA-dominated |
| `blacklist_ratio < 0.05` | upper | remove mapping artifacts (skip for ATACv2) |

## Common pitfalls

1. **Too strict filter.** If less than 70% of cells pass, thresholds are too tight for this dataset. Relax them.
2. **Too permissive filter.** If clusters look messy or many contain fewer than 20 cells, doublets or noise are getting through. Tighten thresholds.
3. **Framework scale mismatch.** Signac TSS versus SnapATAC2 TSSE use the same name but different formulas. Do not copy thresholds across frameworks.
4. **Independent nucleosome and TSS.** Good nucleosome signal does not guarantee good TSS enrichment, and vice versa. Both must clear their thresholds.
5. **Cellranger metadata changes.** ATACv2's `singlecell.csv` has different column names and values than ATACv1. Some columns like `blacklist_region_fragments` may be all zeros. Always inspect the metadata columns.

## References

- ENCODE ATAC-seq standards: https://www.encodeproject.org/atac-seq/
- Signac QC vignette: https://stuartlab.org/signac/articles/pbmc_vignette
- SnapATAC2 QC tutorial: https://kzhang.org/SnapATAC2/tutorials/pbmc.html
