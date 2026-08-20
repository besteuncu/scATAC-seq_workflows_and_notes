"""
07_annotation.py
scATAC-only automated annotation via Gene Activity + CellTypist (+ optional internal Multiome ref)

Prerequisite: scatac_pipeline.ipynb sonuna kadar çalışmış, işlenmiş AnnData diskte.
Ya da: adata = snap.read("pbmc_scatac_processed.h5ad")

WSL Ubuntu, scatac-py conda env içinden çalıştırılacak.
"""

# ---------- Setup ----------
import snapatac2 as snap
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import celltypist
from celltypist import models
import matplotlib.pyplot as plt

sc.settings.set_figure_params(dpi=100, facecolor='white')

# adata yüklü değilse:
# adata = snap.read("/mnt/d/scATAC-seq/data/pbmc_10k_atacv2/pbmc_scatac_processed.h5ad")


# ---------- Adım A: Gene activity matrisi ----------
print("Gene matrix hesaplanıyor (5-10 dk) ...")
gene_matrix = snap.pp.make_gene_matrix(
    adata,
    gene_anno=snap.genome.hg38,
    upstream=2000,
    downstream=0,
    use_x=False
)
# exp-decay ağırlıklı (decay_length=5kb default); upstream=2kb promoter yakalar.
# use_x=False → raw fragment counts kullan. Sonuç: AnnData (cells x genes).

# Kalite: az sayıda hücrede sinyali olan gen'leri at (gürültü)
sc.pp.filter_genes(gene_matrix, min_cells=5)

# CP10K + log1p (CellTypist referansıyla uyumlu normalizasyon)
sc.pp.normalize_total(gene_matrix, target_sum=1e4)
sc.pp.log1p(gene_matrix)


# ---------- Adım B: CellTypist annotation (external referans) ----------
print("CellTypist Immune_All_Low modeli yükleniyor ...")
try:
    model = models.Model.load("Immune_All_Low.pkl")
except Exception:
    models.download_models(force_update=False, model=["Immune_All_Low.pkl"])
    model = models.Model.load("Immune_All_Low.pkl")
# İlk seferde ~10MB indir. Model: 32 immune cell type, logistic regression.

# CellTypist majority_voting için gene_matrix'te clustering olmalı
if 'leiden' not in gene_matrix.obs.columns:
    sc.pp.neighbors(gene_matrix, use_rep=None, n_neighbors=15)
    sc.tl.leiden(gene_matrix, resolution=0.8, random_state=42)
    # neighbors için PCA gerekli — pipeline'da yoksa hesapla:
    # sc.pp.scale(gene_matrix, max_value=10); sc.tl.pca(gene_matrix, n_comps=50)

print("CellTypist annotate çalışıyor (<1 dk) ...")
predictions = celltypist.annotate(
    gene_matrix,
    model=model,
    majority_voting=True,
    over_clustering="leiden"
)
# majority_voting: hücre-başı tahmin + local neighborhood çoğunluk oyu (noise azaltır).

# Ana AnnData'ya (ATAC objesi) taşı
adata.obs['cell_type_celltypist'] = predictions.predicted_labels.majority_voting.values
adata.obs['celltypist_confidence'] = predictions.probability_matrix.max(axis=1)


# ---------- Adım C: Görselleştirme + değerlendirme ----------
sc.pl.umap(
    adata,
    color=['cell_type_celltypist', 'celltypist_confidence'],
    legend_loc='on data',
    legend_fontsize=8,
    save='_scatac_celltypist.png'
)

print("\n--- Hücre tipi dağılımı ---")
print(adata.obs['cell_type_celltypist'].value_counts())

conf_frac = (adata.obs['celltypist_confidence'] > 0.5).mean()
print(f"\nConfidence > 0.5 olan hücre oranı: {100*conf_frac:.1f}%")
# Beklenti: >%75 sağlıklı; <%50 → alternatif metod dene.

adata.write('pbmc_scatac_annotated.h5ad')


# ---------- Adım D (opsiyonel): Internal Multiome referans (ingest / label transfer) ----------
# Multiome muon dataset'i işlenmiş ve annotated ise, scanpy.tl.ingest ile
# etiketleri scATAC'a transfer et. Bu, dataset-spesifik context'i korur.

import os
multi_path = "/mnt/d/scATAC-seq/scATAC-seq_workflows_and_notes/workflows/multiome_pbmc_10k_muon/multiome_annotated.h5mu"

if os.path.exists(multi_path):
    import muon as mu
    print("Internal Multiome referans yükleniyor ...")
    mdata = mu.read(multi_path)
    ref_rna = mdata.mod['rna']
    # Multiome RNA modality; celltype_celltypist ya da leiden var beklenir.

    # Ortak gen seti üzerinden ingest
    common_genes = list(set(ref_rna.var_names) & set(gene_matrix.var_names))
    print(f"Ortak gen sayısı: {len(common_genes)}")

    ref_sub   = ref_rna[:, common_genes].copy()
    query_sub = gene_matrix[:, common_genes].copy()

    # Referans PCA + neighbors + UMAP hazır olmalı; ingest bunu bekler
    if 'X_pca' not in ref_sub.obsm:
        sc.pp.scale(ref_sub, max_value=10)
        sc.tl.pca(ref_sub, n_comps=50)
        sc.pp.neighbors(ref_sub, n_neighbors=15)
        sc.tl.umap(ref_sub)

    # Query'yi referansa yerleştir
    sc.tl.ingest(query_sub, ref_sub, obs='cell_type', embedding_method='umap')
    # cell_type kolonu ref_sub.obs'ta olmalı; yoksa 'leiden' ya da hangi label varsa

    adata.obs['cell_type_internal'] = query_sub.obs['cell_type'].values

    # Karşılaştırma
    print("\n--- External (CellTypist) vs Internal (Multiome) confusion ---")
    print(pd.crosstab(adata.obs['cell_type_celltypist'],
                      adata.obs['cell_type_internal']))

    sc.pl.umap(
        adata,
        color=['cell_type_celltypist', 'cell_type_internal'],
        legend_loc='on data',
        legend_fontsize=8,
        save='_scatac_annotation_comparison.png'
    )

    adata.write('pbmc_scatac_annotated.h5ad')
else:
    print("Multiome referans dosyası bulunamadı; Adım D atlandı.")

print("07_annotation.py tamamlandı.")
