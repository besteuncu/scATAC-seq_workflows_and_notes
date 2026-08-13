# multiome_pipeline.R - Block 1: Setup

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

# Paths (hub-level data)
data_dir <- "D:/scATAC-seq/data/pbmc_multiome_10k"
sample   <- "pbmc_granulocyte_sorted_10k"

h5_path       <- file.path(data_dir, paste0(sample, "_filtered_feature_bc_matrix.h5"))
frag_path     <- file.path(data_dir, paste0(sample, "_atac_fragments.tsv.gz"))
metrics_path  <- file.path(data_dir, paste0(sample, "_per_barcode_metrics.csv"))

stopifnot(file.exists(h5_path))
stopifnot(file.exists(frag_path))
stopifnot(file.exists(metrics_path))

message("Block 1 complete. Paths verified.")

# --- Block 2: Load both modalities ---

message("Loading Multiome H5 (both RNA + ATAC counts) ...")
inputdata.10x <- Read10X_h5(h5_path)

# Split into two matrices
rna_counts  <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# Metadata
metadata <- read.csv(metrics_path, header = TRUE, row.names = 1)

# Build Seurat with RNA first
pbmc <- CreateSeuratObject(counts = rna_counts, meta.data = metadata)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Add ATAC as ChromatinAssay
# Filter peaks to standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use    <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts   <- atac_counts[as.vector(grange.use), ]

# Get gene annotation
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

chrom_assay <- CreateChromatinAssay(
  counts     = atac_counts,
  sep        = c(":", "-"),
  fragments  = frag_path,
  min.cells  = 10,
  annotation = annotation
)
pbmc[["ATAC"]] <- chrom_assay

message(sprintf("Loaded Multiome: %d cells, %d genes, %d peaks",
                ncol(pbmc), nrow(pbmc[["RNA"]]), nrow(pbmc[["ATAC"]])))
pbmc

# --- Block 3: QC per modality ---

# --- RNA side ---
message("Computing RNA QC metrics ...")
# already done: nFeature_RNA, nCount_RNA, percent.mt

# --- ATAC side ---
DefaultAssay(pbmc) <- "ATAC"

message("Computing ATAC QC: NucleosomeSignal + TSSEnrichment ...")
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc, fast = FALSE)

# Metrics preview
message("Per-cell QC summary:")
summary(pbmc@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt",
                           "nCount_ATAC", "TSS.enrichment", "nucleosome_signal")])

# QC violin panel (all metrics)
DefaultAssay(pbmc) <- "RNA"
p_qc <- VlnPlot(
  pbmc,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt",
               "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 3,
  pt.size = 0.1,
  log = TRUE
)
ggsave("qc_violin_multiome.png", p_qc, width = 14, height = 8, dpi = 150)


# --- Block 4: Cell filtering (intersection) ---

n_before <- ncol(pbmc)

pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC       < 100000 &
    nCount_ATAC       > 1000 &
    nCount_RNA        < 25000 &
    nCount_RNA        > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment    > 1 &
    percent.mt        < 20
)

n_after <- ncol(pbmc)
cat(sprintf("Filtered: %d -> %d cells (kept %.1f%%)\n",
            n_before, n_after, 100 * n_after / n_before))

# --- Block 5: Process RNA (SCTransform + PCA) ---

DefaultAssay(pbmc) <- "RNA"

message("Running SCTransform normalization + PCA ...")
pbmc <- SCTransform(pbmc, verbose = FALSE)
pbmc <- RunPCA(pbmc)

# Elbow plot to pick dims
p_elbow <- ElbowPlot(pbmc, ndims = 50)
ggsave("rna_elbow.png", p_elbow, width = 6, height = 4, dpi = 150)

# UMAP on RNA only (for comparison later with joint UMAP)
pbmc <- RunUMAP(pbmc, dims = 1:30, reduction.name = "umap.rna",
                reduction.key = "rnaUMAP_")

# --- Block 6: Process ATAC (TF-IDF + LSI) ---

DefaultAssay(pbmc) <- "ATAC"

message("Running TF-IDF + FindTopFeatures + SVD (LSI) ...")
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")
pbmc <- RunSVD(pbmc)

# UMAP on ATAC only (dims 2:30)
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 2:50,
                reduction.name = "umap.atac",
                reduction.key = "atacUMAP_")

# --- Block 7: WNN integration ---

message("Running WNN (Weighted Nearest Neighbor) ...")
pbmc <- FindMultiModalNeighbors(
  pbmc,
  reduction.list = list("pca", "lsi"),
  dims.list      = list(1:30, 2:50)   # RNA PCA dims 1:30, ATAC LSI dims 2:50
)

# Joint UMAP
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn",
                reduction.name = "wnn.umap",
                reduction.key = "wnnUMAP_")

# Joint clustering
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3,
                     resolution = 0.8, random.seed = 42)

n_clusters <- length(unique(pbmc$seurat_clusters))
message(sprintf("Joint WNN found %d clusters", n_clusters))

# --- Block 8: Joint UMAP + automated annotation (Azimuth PBMC) ---

# Install once: renv::install("satijalab/azimuth")
library(Azimuth)

# Run Azimuth on RNA modality (uses SCT normalized)
DefaultAssay(pbmc) <- "SCT"
pbmc <- RunAzimuth(pbmc, reference = "pbmcref")

# Result fields added to meta.data:
#   predicted.celltype.l1  (broad: T, B, Mono, NK, DC, ...)
#   predicted.celltype.l2  (subtypes: CD4 Naive, CD8 TEM, ...)
#   predicted.celltype.l3  (fine subtypes)
#   predicted.celltype.l2.score  (confidence 0-1)

# Assign cell_type from Azimuth L2 (medium granularity)
pbmc$cell_type <- pbmc$predicted.celltype.l2

# Three UMAPs comparison
p1 <- DimPlot(pbmc, reduction = "umap.rna",
              group.by = "cell_type", label = TRUE, repel = TRUE) +
  NoLegend() + ggtitle("RNA UMAP")
p2 <- DimPlot(pbmc, reduction = "umap.atac",
              group.by = "cell_type", label = TRUE, repel = TRUE) +
  NoLegend() + ggtitle("ATAC UMAP")
p3 <- DimPlot(pbmc, reduction = "wnn.umap",
              group.by = "cell_type", label = TRUE, repel = TRUE) +
  NoLegend() + ggtitle("WNN joint UMAP (Azimuth annotated)")

p_all <- p1 + p2 + p3 + plot_layout(ncol = 3)
ggsave("umap_three_way_annotated.png", p_all, width = 18, height = 6, dpi = 150)

# Confidence overlay
p_conf <- FeaturePlot(pbmc, features = "predicted.celltype.l2.score",
                      reduction = "wnn.umap", min.cutoff = 0.5) +
  ggtitle("Azimuth L2 confidence")
ggsave("azimuth_confidence.png", p_conf, width = 6, height = 5, dpi = 150)

# Save annotated object
saveRDS(pbmc, "multiome_annotated.rds")
