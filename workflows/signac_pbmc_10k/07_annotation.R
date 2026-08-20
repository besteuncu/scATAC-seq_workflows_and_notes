# 07_annotation.R
# scATAC-only automated annotation via Gene Activity + Azimuth (+ optional internal Multiome ref)
# Prerequisite: scatac_pipeline.R sonuna kadar çalışmış olmalı (pbmc objesi ATAC assay + LSI + umap.atac ile hazır)
# Ya da: pbmc <- readRDS(".../pbmc_scatac.rds")

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(EnsDb.Hsapiens.v86)
  library(Azimuth)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

# ---------- Adım A: Gene Activity ----------
DefaultAssay(pbmc) <- "ATAC"
# ATAC assay default: GeneActivity() peak matrisini + fragment file'ı bu assay'den alır.

message("GeneActivity hesaplanıyor (10-15 dk beklenir) ...")
gene.activities <- GeneActivity(pbmc)
# Her gen için TSS-2000..gene_end aralığında fragment sayısı; sparse Matrix (genes x cells) döner.

pbmc[["GeneActivity"]] <- CreateAssayObject(counts = gene.activities)
# Yeni bir assay olarak ekle; RNA-benzeri işleyeceğiz.

pbmc <- NormalizeData(
  object = pbmc,
  assay  = "GeneActivity",
  normalization.method = "LogNormalize",
  scale.factor = median(pbmc$nCount_GeneActivity)
)
# LogNormalize + dataset-adaptive scale factor; Azimuth referansı log-normalized RNA bekler.


# ---------- Adım B: Azimuth (external PBMC referans) ----------
DefaultAssay(pbmc) <- "GeneActivity"
# Azimuth query default assay'i olarak gene activity; RNA gibi işleyecek.

message("RunAzimuth çalışıyor (2-3 dk; ilk seferde referans indirilir ~5 dk) ...")
pbmc <- RunAzimuth(pbmc, reference = "pbmcref", assay = "GeneActivity")
# CCA anchor + label transfer. Sonuç: predicted.celltype.l1/l2/l3 + .score

pbmc$cell_type_atac        <- pbmc$predicted.celltype.l2
pbmc$annotation_confidence <- pbmc$predicted.celltype.l2.score
# L2 (orta granülite) prod etiketi; confidence 0-1


# ---------- Adım C: Görselleştirme + değerlendirme ----------
p_umap <- DimPlot(pbmc, reduction = "umap.atac",
                  group.by = "cell_type_atac",
                  label = TRUE, repel = TRUE) +
  NoLegend() + ggtitle("scATAC-only Azimuth annotation (L2)")

p_conf <- FeaturePlot(pbmc, reduction = "umap.atac",
                      features = "annotation_confidence",
                      min.cutoff = 0.3, max.cutoff = 1.0) +
  ggtitle("Azimuth confidence per cell")

ggsave("annotation_atac_azimuth.png", p_umap + p_conf,
       width = 14, height = 6, dpi = 150)

cat("\n--- Hücre tipi dağılımı ---\n")
print(table(pbmc$cell_type_atac))

cat(sprintf("\nConfidence > 0.5 olan hücre oranı: %.1f%%\n",
            100 * mean(pbmc$annotation_confidence > 0.5)))
# Beklenti: >%75 sağlıklı; <%50 → gene activity gürültülü, alternatif dene.

saveRDS(pbmc, "pbmc_scatac_annotated.rds")


# ---------- Adım D (opsiyonel): Internal Multiome referans ile label transfer ----------
# Multiome pipeline'ı çalıştırdıysan bunu da yap ve karşılaştır.
# Yoksa bu bloğu atla.

multi_path <- "D:/scATAC-seq/scATAC-seq_workflows_and_notes/workflows/multiome_pbmc_10k_signac/multiome_annotated.rds"

if (file.exists(multi_path)) {
  message("Internal Multiome referans yükleniyor ...")
  multi <- readRDS(multi_path)
  DefaultAssay(multi) <- "SCT"
  # Multiome RNA (SCT normalized) referans olacak.

  DefaultAssay(pbmc) <- "GeneActivity"

  transfer.anchors <- FindTransferAnchors(
    reference       = multi,
    query           = pbmc,
    reduction       = "cca",
    reference.assay = "SCT",
    query.assay     = "GeneActivity",
    features        = intersect(VariableFeatures(multi),
                                rownames(pbmc[["GeneActivity"]]))
  )
  # CCA ile ortak low-dim uzay; MNN anchor. features = ortak variable genes.

  predictions <- TransferData(
    anchorset        = transfer.anchors,
    refdata          = multi$cell_type,
    weight.reduction = pbmc[["lsi"]],
    dims             = 2:30
  )
  # weight.reduction = LSI (scATAC'ın native uzayı); dim 1 teknik korelasyon, atla.

  pbmc <- AddMetaData(pbmc, metadata = predictions)
  pbmc$cell_type_internal <- pbmc$predicted.id

  # Karşılaştırma
  cat("\n--- External (Azimuth) vs Internal (Multiome) confusion matrix ---\n")
  print(table(External = pbmc$cell_type_atac,
              Internal = pbmc$cell_type_internal))

  p_int <- DimPlot(pbmc, reduction = "umap.atac",
                   group.by = "cell_type_internal",
                   label = TRUE, repel = TRUE) +
    NoLegend() + ggtitle("Internal Multiome referans annotation")

  ggsave("annotation_atac_internal.png",
         p_umap + p_int, width = 14, height = 6, dpi = 150)

  saveRDS(pbmc, "pbmc_scatac_annotated.rds")
} else {
  message("Multiome referans dosyası bulunamadı; Adım D atlandı.")
}

message("07_annotation.R tamamlandı.")
