# 01. Biology and assays — what scATAC-seq measures

## Where scATAC sits in the regulatory cascade

Gene regulation is a layered process, and each modality reads a different layer:

```
DNA sequence → Chromatin state → TF binding → Transcription → Protein
    |               |                |              |
   WGS         ATAC-seq         ChIP-seq        RNA-seq
```

| Assay | Layer measured | Question answered |
|-------|----------------|-------------------|
| WGS | Nucleotide sequence | What is the code? |
| **scATAC-seq** | **Open chromatin (Tn5 accessible)** | **Which code is potentially readable in each cell?** |
| ChIP-seq / CUT&RUN | Occupancy of one TF or histone mark | Where is this specific protein bound? |
| scRNA-seq | mRNA levels per cell | What was transcribed? |

scATAC uniquely reads **regulatory potential** (which enhancers are accessible for any TF to bind), while ChIP measures **realization** (where one specific TF actually sat), and RNA reads the downstream product.

## The biology of Tn5 tagmentation

**ATAC** = Assay for Transposase-Accessible Chromatin.

Chromatin exists in two functional states:
- **Closed (heterochromatin):** nucleosomes tightly packed, DNA inaccessible, TF binding blocked, transcription silent.
- **Open (euchromatin):** nucleosome-free regions, DNA available for TF binding, transcription possible.

Tn5 transposase enzyme:
1. Can only physically access **open chromatin** (blocked by tightly packed histones).
2. Cuts accessible DNA and inserts two sequencing adapters.
3. Each cut event produces one **fragment** that is later sequenced.

The result: a genome-wide map of open chromatin per cell.

## Fragment file format

The `fragments.tsv.gz` output stores one Tn5 cut event per line:

```
chr1    30500    30720    AAACGAAAGAAACGAA-1    2
```

Columns:
- `chr1` — chromosome
- `30500` — fragment start position
- `30720` — fragment end position (220 bp long)
- `AAACGAAAGAAACGAA-1` — cell barcode (10x)
- `2` — read count (usually 1, occasionally 2 with PCR duplicates)

A typical cell yields 10,000–100,000 fragments. Each fragment is direct evidence that a genomic region was accessible in that cell.

## Peaks: biologically meaningful regions

Individual fragments are not interpreted in isolation. **Peaks** are genomic regions where fragments concentrate, indicating regulatory elements (enhancers, promoters, TF binding sites).

**Peak calling** aggregates fragments across cells and detects local density enrichments against background. The standard tool is **MACS2** (Model-based Analysis of ChIP-Seq).

Typical output: a `peaks.bed` file with 50,000–200,000 peaks for human PBMC datasets. Each peak spans a few hundred base pairs.

Terminology:
- **Peak** = biological object (a regulatory element)
- **Fragment** = physical measurement (a Tn5 cut event)

## The cell × peak matrix

Once peaks are called, per-cell fragment counts inside each peak form a matrix:

```
              peak_1   peak_2   peak_3   ...
cell_AAA-1      2        0        1
cell_AAB-1      1        1        0
cell_AAC-1      0        0        2
```

This is the direct analog of the cell × gene expression matrix in scRNA-seq. Now:
- **Feature** = peak (not gene)
- **Cell** = 10x barcode
- **Value** = fragment count (typically 0/1/2, since a diploid cell has at most two accessible alleles per site)

All downstream analysis (normalization, dimensionality reduction, clustering, annotation) operates on this sparse binary-like matrix.

## Why single-cell resolution matters

Bulk ATAC-seq returns average signal per sample. Two biologically distinct scenarios collapse to the same measurement:

- **Scenario A:** A peak is half-open in every cell.
- **Scenario B:** A peak is fully open in half the cells and closed in the other half.

Bulk cannot distinguish these. scATAC-seq recovers the per-cell distribution, revealing **cell-type-specific regulatory landscape**. A T cell's active enhancers differ from a B cell's, and single-cell resolution captures this directly.

## scATAC vs scRNA — comparing the two layers

| Dimension | scRNA-seq | scATAC-seq |
|-----------|-----------|------------|
| Unit measured | UMI count (mRNA molecule) | Fragment (Tn5 cut event) |
| Feature | Gene (~20K) | Peak (~50–200K) |
| Reads | Expression | Chromatin accessibility |
| Physical meaning | Realized expression | Regulatory potential |
| Position in cascade | Downstream (result) | Upstream (cause) |
| Typical value range | 0–50+ per gene | 0/1/2 (diploid limit) |
| Statistical distribution | Negative binomial | Bernoulli / binary fits better |
| Sparsity | 60–80% zeros | 95–98% zeros |

The critical biological distinction:
- scRNA count: "this gene is actively transcribed right now"
- scATAC fragment: "this region is accessible enough for TF binding"

Combining both modalities per cell (the Multiome case) exposes four states:
- ATAC open + RNA present → gene actively expressed
- **ATAC open + RNA absent → priming** (chromatin ready, transcription not yet started). A hallmark of developmental cell states.
- ATAC closed + RNA absent → silent, canonical
- ATAC closed + RNA present → rare, likely measurement error or declining state

Priming is what scATAC uniquely detects when paired with RNA. It is why differentiation, regeneration, and cancer plasticity studies increasingly rely on scATAC or Multiome.

## References

- Buenrostro et al. 2013, *Nat Methods* — ATAC-seq original method.
- Buenrostro et al. 2015, *Nature* — first single-cell ATAC-seq.
- Cusanovich et al. 2018, *Cell* — single-cell atlas of chromatin accessibility.
- Stuart et al. 2021, *Nat Methods* — Signac framework.
- Zhang et al. 2024, *Nat Methods* — SnapATAC2 framework.
