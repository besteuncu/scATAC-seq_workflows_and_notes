# install_packages.R
# One-time install of all R dependencies for the Signac PBMC 10k pipeline.
# Run after `renv::init(bare = TRUE)` in the workflow directory.

# BiocManager: R'ın Bioconductor paketlerine köprü. CRAN'da yok, önce kurulur.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# ---------- CRAN packages (general R ecosystem) ----------
cran_pkgs <- c(
  "Seurat",       # single-cell analysis framework (base for Signac)
  "Signac",       # scATAC-seq extension of Seurat
  "ggplot2",      # plotting (Seurat depends on it, explicit for clarity)
  "patchwork",    # combine multiple ggplots into one figure
  "hdf5r",        # read .h5 files (cellranger peak matrix output)
  "here"          # project-root-relative path management
)

# ---------- Bioconductor packages (genomics + motif analysis) ----------
bioc_pkgs <- c(
  # Genome annotation and sequence
  "EnsDb.Hsapiens.v86",           # Ensembl gene annotation for hg38
  "BSgenome.Hsapiens.UCSC.hg38",  # hg38 genome sequence (~800 MB)
  "GenomeInfoDb",                 # chromosome naming conventions
  "GenomicRanges",                # interval algebra on genome
  
  # Motif analysis (ChromVAR + JASPAR)
  "chromVAR",       # per-cell TF motif activity deviations
  "motifmatchr",    # scan peaks for TF binding motifs
  "JASPAR2020",     # TF motif position weight matrix database
  "TFBSTools"       # motif utility functions
)

# ---------- Install ----------
install.packages(cran_pkgs)
BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)

cat("\nInstallation complete. Now run:\n")
cat("  renv::snapshot()\n")
cat("to lock package versions to renv.lock.\n")