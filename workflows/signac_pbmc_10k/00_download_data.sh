#!/usr/bin/env bash
# Download 10x Genomics ATACv2 PBMC 10k dataset (hg38, Cellranger-ATAC v2, 2022).

set -euo pipefail

# Resolve hub-level data dir (three levels up from workflow folder)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="$SCRIPT_DIR/../../../data/pbmc_10k_atacv2"
mkdir -p "$DATA_DIR"

BASE="https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller"

FILES=(
  "10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz"
  "10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz.tbi"
  "10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5"
  "10k_pbmc_ATACv2_nextgem_Chromium_Controller_singlecell.csv"
)

for f in "${FILES[@]}"; do
  target="$DATA_DIR/$f"
  if [[ -f "$target" ]]; then
    echo "[skip] $f already exists"
  else
    echo "[download] $f"
    curl -L -o "$target" "$BASE/$f"
  fi
done

echo
echo "Download complete. Files in $DATA_DIR:"
ls -lh "$DATA_DIR"