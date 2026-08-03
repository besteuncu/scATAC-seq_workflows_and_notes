#!/usr/bin/env bash
# Download the 10x Genomics PBMC 10k scATAC v1.0.1 dataset.
# Usage: bash 00_download_data.sh
# Expected total size: about 2 GB.

set -euo pipefail

DATA_DIR="data/atac_v1_pbmc_10k"
mkdir -p "$DATA_DIR"

BASE="https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k"

FILES=(
  "atac_v1_pbmc_10k_fragments.tsv.gz"
  "atac_v1_pbmc_10k_fragments.tsv.gz.tbi"
  "atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
  "atac_v1_pbmc_10k_singlecell.csv"
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