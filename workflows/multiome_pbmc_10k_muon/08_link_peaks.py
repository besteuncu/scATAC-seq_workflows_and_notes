"""
08_link_peaks.py
Multiome peak-gene linkage — Signac LinkPeaks'in manuel Python implementasyonu.

Prerequisite: multiome_pipeline.ipynb sonuna kadar çalışmış (multiome_annotated.h5mu var)

Yaklaşım (Signac LinkPeaks'i takip):
  1. Her gen için ±500kb içindeki aday peak'leri bul
  2. Her (peak, gen) çifti için Pearson korelasyon hesapla
  3. GC + accessibility matched 200 background peak → null dağılım
  4. z-score, p-value, filtre

Bu implementasyon eğitim amaçlı — Signac'in Rcpp optimizasyonu kadar hızlı değil.
Büyük dataset'te SCENIC+ ya da Signac tercih edilir.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import muon as mu
import anndata as ad
from scipy.stats import pearsonr, norm
from scipy.sparse import issparse
from tqdm import tqdm
import pyranges as pr

# İsteğe bağlı: GC hesabı için genom fasta gerekir. pyfaidx ile:
# from pyfaidx import Fasta
# genome = Fasta("/path/to/hg38.fa")

# ---------- Config ----------
MULTIOME_PATH = "multiome_annotated.h5mu"
DISTANCE      = 500_000       # ±500 kb pencere
MIN_CELLS     = 10            # peak en az 10 hücrede
N_BACKGROUND  = 200           # background peak sayısı
PVAL_CUTOFF   = 0.05
SCORE_CUTOFF  = 0.05
GENES_USE     = None          # None = HVG'leri kullan; test için: ["CD8A","MS4A1","CD4","GNLY","LYZ","IL7R"]


# ---------- 1) Multiome dataset yükle ----------
print("Multiome dataset yükleniyor ...")
mdata = mu.read(MULTIOME_PATH)
rna  = mdata.mod['rna']
atac = mdata.mod['atac']

# Ortak hücreler
common_cells = list(set(rna.obs_names) & set(atac.obs_names))
rna  = rna[common_cells].copy()
atac = atac[common_cells].copy()
print(f"Ortak hücre sayısı: {len(common_cells)}")


# ---------- 2) Gen koordinatları (Ensembl / GTF) ----------
# Gerçek dünyada bir GTF ya da EnsDb'nin Python muadili (pyensembl) gerekir.
# Basitleştirmek için: RNA AnnData'da var['chrom','start','end'] varsa kullan;
# yoksa gtfparse ile GTF oku.
if not all(c in rna.var.columns for c in ['chrom', 'start', 'end']):
    print("UYARI: rna.var'da chrom/start/end yok.")
    print("Devam etmek için pyensembl ile gen koordinatları eklenmeli:")
    print("  from pyensembl import EnsemblRelease")
    print("  ens = EnsemblRelease(98)  # hg38, ilk sefer indirir")
    print("  ens.download(); ens.index()")
    print("  def get_coords(sym):")
    print("      try:")
    print("          g = ens.genes_by_name(sym)[0]")
    print("          return pd.Series([g.contig, g.start, g.end])")
    print("      except: return pd.Series([None,None,None])")
    print("  rna.var[['chrom','start','end']] = rna.var_names.to_series().apply(get_coords)")
    raise SystemExit("Gen koordinatı eksik, script durduruldu.")

gene_coords = rna.var[['chrom', 'start', 'end']].copy()
gene_coords['gene']  = rna.var_names
gene_coords['chrom'] = gene_coords['chrom'].astype(str)
gene_coords['tss']   = gene_coords['start']    # + strand kabul; - strand için end kullanılabilir


# ---------- 3) Peak koordinatları ----------
# atac.var_names formatı: "chr1-1000-2000" veya "chr1:1000-2000"
def parse_peak(name):
    for sep in [':', '-']:
        if sep in name:
            parts = name.replace(':', '-').split('-')
            if len(parts) >= 3:
                return parts[0], int(parts[1]), int(parts[2])
    return None, None, None

peak_coords = pd.DataFrame(
    [parse_peak(p) for p in atac.var_names],
    columns=['chrom', 'start', 'end'],
    index=atac.var_names
)
peak_coords['peak']   = atac.var_names
peak_coords['gc']     = np.nan  # aşağıda hesaplanacak (isteğe bağlı)
peak_coords['access'] = np.asarray((atac.X > 0).sum(axis=0)).flatten()  # kaç hücrede accessible


# ---------- 4) Peak filtre (min.cells) ----------
peak_ok = peak_coords['access'] >= MIN_CELLS
peaks_use = peak_coords[peak_ok].copy()
atac_use = atac[:, peak_ok.values].copy()
print(f"Peak filtresi: {peak_coords.shape[0]} → {peaks_use.shape[0]}")


# ---------- 5) Gen seti ----------
if GENES_USE is None:
    sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor='seurat_v3', layer='counts')
    genes_test = rna.var_names[rna.var['highly_variable']].tolist()
else:
    genes_test = [g for g in GENES_USE if g in rna.var_names]
genes_test = [g for g in genes_test if g in gene_coords.index]
print(f"Test edilecek gen sayısı: {len(genes_test)}")


# ---------- 6) Yardımcı fonksiyonlar ----------
def get_expr_vec(gene):
    """rna assay'inden gen ekspresyon vektörü (dense)"""
    x = rna[:, gene].X
    return np.asarray(x.todense()).flatten() if issparse(x) else np.asarray(x).flatten()

def get_access_vec(peak):
    """atac assay'inden peak accessibility vektörü (dense)"""
    x = atac_use[:, peak].X
    return np.asarray(x.todense()).flatten() if issparse(x) else np.asarray(x).flatten()

def find_candidate_peaks(gene, distance=DISTANCE):
    """gen TSS'inden ±distance içindeki peak'leri döndür"""
    gc = gene_coords.loc[gene]
    mask = (
        (peaks_use['chrom'] == gc['chrom']) &
        (peaks_use['end']   > gc['tss'] - distance) &
        (peaks_use['start'] < gc['tss'] + distance)
    )
    return peaks_use[mask].index.tolist()

def match_background(peak, n=N_BACKGROUND, use_gc=False):
    """
    GC + accessibility matched background peak seç.
    use_gc=False → sadece accessibility matching (GC olmadan basitleştirilmiş).
    """
    target_access = peaks_use.loc[peak, 'access']
    # +/- %20 access bandında olan peak'ler
    lo, hi = target_access * 0.8, target_access * 1.2
    pool = peaks_use[(peaks_use['access'] >= lo) & (peaks_use['access'] <= hi)].index
    pool = pool[pool != peak]
    if len(pool) < n:
        return pool.tolist()
    return np.random.choice(pool, size=n, replace=False).tolist()


# ---------- 7) Ana döngü — link hesabı ----------
print("Peak-gene link hesabı başlıyor ...")
results = []

for gene in tqdm(genes_test, desc="Genes"):
    candidates = find_candidate_peaks(gene)
    if len(candidates) == 0:
        continue

    expr = get_expr_vec(gene)
    if expr.std() == 0:
        continue

    for peak in candidates:
        access = get_access_vec(peak)
        if access.std() == 0:
            continue

        # Gerçek korelasyon
        r_actual, _ = pearsonr(access, expr)

        # Background null dağılım
        bg_peaks = match_background(peak)
        bg_rs = []
        for bp in bg_peaks:
            bg_access = get_access_vec(bp)
            if bg_access.std() == 0:
                continue
            r, _ = pearsonr(bg_access, expr)
            bg_rs.append(r)

        if len(bg_rs) < 20:
            continue

        bg_mean = np.mean(bg_rs)
        bg_sd   = np.std(bg_rs)
        if bg_sd == 0:
            continue

        z = (r_actual - bg_mean) / bg_sd
        pval = 2 * norm.sf(abs(z))

        results.append({
            'peak': peak,
            'gene': gene,
            'chrom': peaks_use.loc[peak, 'chrom'],
            'start': peaks_use.loc[peak, 'start'],
            'end':   peaks_use.loc[peak, 'end'],
            'score': r_actual,
            'zscore': z,
            'pvalue': pval
        })

links_df = pd.DataFrame(results)


# ---------- 8) Multiple testing correction + filtre ----------
if len(links_df) > 0:
    from statsmodels.stats.multitest import multipletests
    _, qvals, _, _ = multipletests(links_df['pvalue'], method='fdr_bh')
    links_df['qvalue'] = qvals

    links_filt = links_df[
        (links_df['pvalue'] < PVAL_CUTOFF) &
        (links_df['score'].abs() > SCORE_CUTOFF)
    ].copy()
    links_filt = links_filt.sort_values('zscore', key=abs, ascending=False)

    print(f"\nToplam test edilen peak-gene çifti: {len(links_df)}")
    print(f"Filtreden geçen anlamlı link: {len(links_filt)}")
    print("\n--- En güçlü 20 link ---")
    print(links_filt.head(20))

    links_filt.to_csv("peak_gene_links.csv", index=False)
    print("\npeak_gene_links.csv kaydedildi.")

    # Özet
    print("\nHer gen başına link sayısı — top 20:")
    print(links_filt['gene'].value_counts().head(20))
else:
    print("Hiç anlamlı link bulunamadı — parametreleri gevşet ya da gen setini genişlet.")


# ---------- 9) Basit görselleştirme ----------
# Not: muon'da Signac'in CoveragePlot muadili yok — pyGenomeTracks veya IGV gerekir.
# Buradan sonrası isteğe bağlı; sadece link istatistiklerini plot ediyoruz.
if len(links_df) > 0:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    axes[0].hist(links_df['score'], bins=50)
    axes[0].axvline(0, color='k', ls='--')
    axes[0].set_xlabel('Pearson r'); axes[0].set_title('Correlation distribution')

    axes[1].hist(links_df['zscore'], bins=50)
    axes[1].axvline(0, color='k', ls='--')
    axes[1].set_xlabel('z-score'); axes[1].set_title('Background-corrected z')

    # Uzaklık dağılımı (TSS'e)
    dist = []
    for _, row in links_df.iterrows():
        tss = gene_coords.loc[row['gene'], 'tss']
        d = (row['start'] + row['end']) / 2 - tss
        dist.append(d)
    axes[2].hist(dist, bins=50)
    axes[2].axvline(0, color='k', ls='--')
    axes[2].set_xlabel('Peak center − TSS (bp)')
    axes[2].set_title('Peak-TSS distance')

    plt.tight_layout()
    plt.savefig("link_stats.png", dpi=150)
    print("link_stats.png kaydedildi.")

print("08_link_peaks.py tamamlandı.")
