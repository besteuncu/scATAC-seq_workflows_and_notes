# 08_link_peaks.R
# Multiome peak-gene linkage via Signac LinkPeaks
# Prerequisite: multiome_pipeline.R sonuna kadar çalıştırılmış (multiome_annotated.rds var)

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(BSgenome.Hsapiens.UCSC.hg38)   # GC içeriği için genom sekansı
  library(EnsDb.Hsapiens.v86)            # gen koordinatları
  library(ggplot2)
})

set.seed(42)

# ---------- Objeyi yükle ----------
pbmc <- readRDS("multiome_annotated.rds")

# ---------- Adım 1: RegionStats (GC + accessibility metrics) ----------
DefaultAssay(pbmc) <- "ATAC"

message("RegionStats hesaplanıyor (1-2 dk) ...")
pbmc <- RegionStats(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)
# Her peak için GC%, sequence length, dinucleotide freq. LinkPeaks background
# eşleştirmesi bu metriklere göre yapılıyor.

# ---------- Adım 2: LinkPeaks ----------
message("LinkPeaks çalışıyor (15-30 dk) ...")
pbmc <- LinkPeaks(
  object           = pbmc,
  peak.assay       = "ATAC",
  expression.assay = "SCT",
  genes.use        = NULL,     # NULL = variable genler; belirli set için: c("CD8A","MS4A1",...)
  distance         = 5e5,      # ±500 kb pencere
  min.cells        = 10,       # peak en az 10 hücrede accessible
  method           = "pearson",
  n_sample         = 200,      # her peak için 200 background peak (GC + access matched)
  pvalue_cutoff    = 0.05,
  score_cutoff     = 0.05
)
# İç akış: (peak, gen) çifti ±500kb → Pearson r → GC/access matched 200 background
# → null dağılım → z-score → p-value → filtre.

# ---------- Adım 3: Sonuçları incele ----------
links_gr <- Links(pbmc)
cat(sprintf("\nToplam link sayısı: %d\n", length(links_gr)))

links_df <- as.data.frame(links_gr)
links_df <- links_df[order(-abs(links_df$zscore)), ]

cat("\n--- En güçlü 20 link ---\n")
print(head(links_df, 20))

write.csv(links_df, "peak_gene_links.csv", row.names = FALSE)
cat("\npeak_gene_links.csv kaydedildi.\n")

# Özet istatistik
cat(sprintf("\nLink uzaklık istatistiği (bp):\n"))
print(summary(links_df$width))
# Kısa link'ler promoter-yakın (güvenilir); uzun link'ler distal enhancer.

cat(sprintf("\nHer gene başına link sayısı — top 20:\n"))
print(head(sort(table(links_df$gene), decreasing = TRUE), 20))

# ---------- Adım 4: Görselleştirme (coverage plot + link arkları) ----------
DefaultAssay(pbmc) <- "ATAC"

marker_genes <- c("CD8A", "MS4A1", "CD4", "GNLY", "LYZ", "IL7R")

for (gene in marker_genes) {
  message(sprintf("Coverage plot: %s", gene))
  p <- tryCatch({
    CoveragePlot(
      object            = pbmc,
      region            = gene,
      features          = gene,
      expression.assay  = "SCT",
      extend.upstream   = 5e4,
      extend.downstream = 5e4,
      group.by          = "cell_type"
    )
  }, error = function(e) {
    message(sprintf("  Skip %s: %s", gene, e$message))
    NULL
  })
  if (!is.null(p)) {
    ggsave(paste0("coverage_", gene, ".png"), p,
           width = 12, height = 8, dpi = 150)
  }
}
# Panel dizini: (üst) gen anotasyonu + link arkları, (orta) peak coverage cell type başına,
# (alt) gen ekspresyon violin cell type başına.

# ---------- Kaydet ----------
saveRDS(pbmc, "multiome_annotated_linked.rds")
message("08_link_peaks.R tamamlandı.")
