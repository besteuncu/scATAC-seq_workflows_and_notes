# 01_qc.R
# Signac PBMC 10k ATACv2 pipeline — QC stage
# Dataset: 10x Genomics ATACv2 PBMC 10k (hg38, Cellranger-ATAC v2, March 2022)

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(EnsDb.Hsapiens.v86)
  library(GenomeInfoDb)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

# --- Paths (shared with SnapATAC2 workflow at hub-level D:/scATAC-seq/data/) ---
data_dir <- "D:/scATAC-seq/data/pbmc_10k_atacv2"
sample   <- "10k_pbmc_ATACv2_nextgem_Chromium_Controller"

counts_path   <- file.path(data_dir, paste0(sample, "_filtered_peak_bc_matrix.h5"))
metadata_path <- file.path(data_dir, paste0(sample, "_singlecell.csv"))
frag_path     <- file.path(data_dir, paste0(sample, "_fragments.tsv.gz"))

# --- Output directories (workflow-local) ---
dir.create("data/rds", showWarnings = FALSE, recursive = TRUE)
dir.create("results",  showWarnings = FALSE, recursive = TRUE)

# --- Sanity check ---
stopifnot(file.exists(counts_path))
stopifnot(file.exists(metadata_path))
stopifnot(file.exists(frag_path))

message("Block 1 complete. Paths verified, output dirs ready.")


# --- Block 2: Load data + build ChromatinAssay ---

message("Loading peak matrix (.h5) ...")
counts <- Read10X_h5(counts_path)

message("Loading per-cell metadata (.csv) ...")
metadata <- read.csv(metadata_path, header = TRUE, row.names = 1)

message("Building ChromatinAssay (fragments linked) ...")
chrom_assay <- CreateChromatinAssay(
  counts       = counts,
  sep          = c(":", "-"),
  fragments    = frag_path,
  genome       = "hg38",
  min.cells    = 10,
  min.features = 200
)

message("Wrapping in Seurat object ...")
pbmc <- CreateSeuratObject(
  counts    = chrom_assay,
  assay     = "peaks",
  meta.data = metadata
)

message(sprintf("Loaded %d cells x %d peaks (pre-filter).", ncol(pbmc), nrow(pbmc)))
pbmc


# --- Block 3: Add annotation + compute QC metrics ---

message("Loading Ensembl v86 (hg38) gene annotation ...")
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"   # Ensembl "1" -> UCSC "chr1"
Annotation(pbmc) <- annotation

message("Computing NucleosomeSignal ...")
pbmc <- NucleosomeSignal(pbmc)

message("Computing TSSEnrichment (fast=FALSE for higher-resolution TSS profile) ...")
pbmc <- TSSEnrichment(pbmc, fast = FALSE)

message("Computing pct_reads_in_peaks + blacklist_ratio ...")
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio    <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

# Categorical grouping for QC plots (used in Block 4)
pbmc$high.tss         <- ifelse(pbmc$TSS.enrichment > 2, "High", "Low")
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, "High (>4)", "Low (<=4)")

message("QC metrics computed. Summary of key metrics:")
summary(pbmc@meta.data[, c("nCount_peaks", "TSS.enrichment", "nucleosome_signal",
                           "pct_reads_in_peaks", "blacklist_ratio")])



# --- Block 4: QC visualization + cell filtering ---

message("Building QC violin panel ...")

p_qc <- VlnPlot(
  object   = pbmc,
  features = c("nCount_peaks", "TSS.enrichment", "nucleosome_signal",
               "pct_reads_in_peaks"),
  pt.size  = 0.1,
  ncol     = 4
)
ggsave("results/qc_violin.png", p_qc, width = 16, height = 4, dpi = 150)

message("Building fragment length histogram ...")
p_frag <- FragmentHistogram(pbmc, group.by = "nucleosome_group")
ggsave("results/qc_fragment_histogram.png", p_frag, width = 8, height = 4, dpi = 150)

message("Cell counts before filtering:")
n_before <- ncol(pbmc)
cat(sprintf("  Total: %d\n", n_before))

message("Applying QC filter (Signac vignette thresholds, ATACv2-adjusted) ...")
pbmc <- subset(
  x = pbmc,
  subset = nCount_peaks > 3000 &
    nCount_peaks < 100000 &
    pct_reads_in_peaks > 40 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

n_after <- ncol(pbmc)
cat(sprintf("Cell counts after filtering: %d (kept %.1f%%)\n",
            n_after, 100 * n_after / n_before))

message("Post-filter summary:")
summary(pbmc@meta.data[, c("nCount_peaks", "TSS.enrichment",
                           "nucleosome_signal", "pct_reads_in_peaks")])


# --- Block 5: Feature selection + LSI + Clustering + UMAP ---

message("TF-IDF normalization ...")
pbmc <- RunTFIDF(pbmc)

message("Selecting top variable features (q5 = top 95% by count) ...")
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q5")

message("Running SVD (LSI) ...")
pbmc <- RunSVD(pbmc)

# Visualize LSI-1 depth correlation
message("DepthCor plot ...")
p_depth <- DepthCor(pbmc) + ggtitle("LSI component vs sequencing depth")
ggsave("results/lsi_depthcor.png", p_depth, width = 6, height = 4, dpi = 150)

lsi1_cor <- cor(pbmc@reductions$lsi@cell.embeddings[, 1], pbmc$nCount_peaks)
message(sprintf("LSI-1 correlation with nCount_peaks: %.3f", lsi1_cor))
message("If |cor| > 0.7, dropping LSI-1 (dims 2:30) is mandatory.")

message("UMAP on LSI dims 2:30 ...")
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 2:30, seed.use = 42)

message("KNN graph + Leiden clustering (resolution 0.8) ...")
pbmc <- FindNeighbors(pbmc, reduction = "lsi", dims = 2:30)
pbmc <- FindClusters(pbmc, algorithm = 3, resolution = 0.8, random.seed = 42)

n_clusters <- length(unique(pbmc$seurat_clusters))
message(sprintf("Found %d clusters", n_clusters))

message("UMAP plot with cluster labels ...")
p_umap <- DimPlot(pbmc, label = TRUE, pt.size = 0.4) + NoLegend() +
  ggtitle(sprintf("Leiden clusters on LSI (n=%d, res=0.8)", n_clusters))
ggsave("results/umap_clusters.png", p_umap, width = 6, height = 5, dpi = 150)

# Save clustered object for reuse in later blocks
saveRDS(pbmc, "data/rds/03_clustered.rds")
message("Saved clustered object to data/rds/03_clustered.rds")


# --- Block 6: Cell type annotation via gene activity ---

message("Computing gene activity (peak sum over gene body + promoter) ...")
gene.activities <- GeneActivity(pbmc)

# Add as a new RNA assay
pbmc[["RNA"]] <- CreateAssayObject(counts = gene.activities)

message("Normalizing gene activity ...")
pbmc <- NormalizeData(
  object              = pbmc,
  assay               = "RNA",
  normalization.method = "LogNormalize",
  scale.factor        = median(pbmc$nCount_RNA)
)

# Switch default assay to RNA for marker plotting
DefaultAssay(pbmc) <- "RNA"

# Canonical PBMC marker genes
marker_genes <- c(
  "MS4A1",   # B cell (CD20)
  "CD3E",    # T cell (pan)
  "CD8A",    # CD8 T
  "CD4",     # CD4 T (also expressed on monocytes)
  "NKG7",    # NK
  "CD14",    # Classical monocyte
  "LYZ",     # Monocyte / DC
  "FCGR3A"   # CD16 monocyte / NK
)

# Keep only markers found in the assay
marker_genes <- intersect(marker_genes, rownames(pbmc))
message(sprintf("Plotting %d markers", length(marker_genes)))

p_markers <- FeaturePlot(
  object     = pbmc,
  features   = marker_genes,
  pt.size    = 0.3,
  max.cutoff = "q95",
  ncol       = 4
) & NoAxes()
ggsave("results/annotation_markers.png", p_markers, width = 16, height = 8, dpi = 150)

DefaultAssay(pbmc) <- "peaks"

# Save annotated object
saveRDS(pbmc, "data/rds/04_annotated.rds")
message("Saved annotated object to data/rds/04_annotated.rds")