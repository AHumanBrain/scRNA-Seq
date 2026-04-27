# Create conda environment
conda create -n scrna python=3.10 -y
conda activate scrna

# Install Scanpy and dependencies
pip install scanpy matplotlib seaborn pandas numpy leidenalg

# Verify
python -c "import scanpy; print(scanpy.__version__)"

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

sc.settings.verbosity = 3
sc.settings.figdir = './figures/'

# Download data
import urllib.request, tarfile, os
url = 'https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz'
urllib.request.urlretrieve(url, 'pbmc3k.tar.gz')
with tarfile.open('pbmc3k.tar.gz') as tar:
    tar.extractall('.')

# Load into AnnData object
adata = sc.read_10x_mtx(
    'filtered_gene_bc_matrices/hg19/',
    var_names='gene_symbols',
    cache=True
)
adata.var_names_make_unique()

print(adata)
# AnnData object with n_obs x n_vars = 2700 x 32738

# Calculate QC metrics
# Mitochondrial genes start with "MT-"
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                            log1p=False, inplace=True)

# Visualize QC metrics
sc.pl.violin(adata,
             ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

# Filter cells:
# - At least 200 genes detected
# - No more than 2500 genes (doublet filter)
# - Less than 5% mitochondrial reads (dead cell filter)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

print(f"Cells after QC: {adata.n_obs}")
print(f"Genes after QC: {adata.n_vars}")

# Normalize: scale each cell to 10,000 total counts, then log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Save raw normalized counts for later (needed for DE analysis)
adata.raw = adata

# Find highly variable genes (top 2000)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3,
                             min_disp=0.5)
print(f"Highly variable genes: {adata.var.highly_variable.sum()}")

# Plot variable genes
sc.pl.highly_variable_genes(adata)

# Keep only highly variable genes for downstream analysis
adata = adata[:, adata.var.highly_variable]

# Scale data (zero mean, unit variance) — needed for PCA
sc.pp.scale(adata, max_value=10)

# Run PCA
sc.tl.pca(adata, svd_solver='arpack')

# Plot variance explained by each PC
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)
# Look for the "elbow" — where adding more PCs stops helping

# Visualize cells in PCA space
sc.pl.pca(adata, color='CST3')  # CST3 = monocyte marker

# Build neighborhood graph using top 10 PCs
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Cluster using Leiden algorithm (better than Louvain)
sc.tl.leiden(adata, resolution=0.5)
print(f"Number of clusters: {adata.obs['leiden'].nunique()}")

# Run UMAP for visualization
sc.tl.umap(adata)

# Plot UMAP colored by cluster
sc.pl.umap(adata, color=['leiden'], legend_loc='on data',
           title='PBMC clusters', frameon=False)

# Find marker genes for each cluster
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# Get top markers as a dataframe
markers = sc.get.rank_genes_groups_df(adata, group=None)
print(markers.head(20))

# Dot plot of top markers per cluster
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5)

# Known PBMC marker genes:
# CD3D, CD3E, CD3G  → T cells (all T cells)
# CD4               → CD4+ T cells (helper T cells)
# CD8A, CD8B        → CD8+ T cells (cytotoxic T cells)
# MS4A1, CD79A      → B cells
# GNLY, NKG7        → NK cells
# CD14, LYZ         → CD14+ Monocytes
# FCGR3A, MS4A7     → CD16+ Monocytes
# FCER1A, CST3      → Dendritic cells
# PPBP              → Megakaryocytes (platelets)

# Visualize known markers on UMAP
marker_genes = ['CD3D', 'CD4', 'CD8A', 'MS4A1', 'GNLY',
                'CD14', 'LYZ', 'FCGR3A', 'FCER1A', 'PPBP']
sc.pl.umap(adata, color=marker_genes, ncols=5)

# Assign cell type labels to clusters
# (adjust based on YOUR cluster numbers)
cell_type_map = {
    '0': 'CD4+ T cells',
    '1': 'CD14+ Monocytes',
    '2': 'CD4+ T cells',
    '3': 'B cells',
    '4': 'CD8+ T cells',
    '5': 'NK cells',
    '6': 'CD14+ Monocytes',
    '7': 'Dendritic cells',
    '8': 'CD16+ Monocytes',
}
adata.obs['cell_type'] = adata.obs['leiden'].map(cell_type_map)

# Final annotated UMAP
sc.pl.umap(adata, color='cell_type', legend_loc='on data',
           title='PBMC Cell Types', frameon=False,
           palette='tab10')

# Python: Compare CD4 T cells vs CD8 T cells
sc.tl.rank_genes_groups(adata, 'cell_type', groups=['CD4+ T cells'],
                         reference='CD8+ T cells', method='wilcoxon')
sc.pl.rank_genes_groups(adata, groups=['CD4+ T cells'], n_genes=20)

# Python
adata.write('pbmc3k_analyzed.h5ad')
# Load later: adata = sc.read_h5ad('pbmc3k_analyzed.h5ad')

# Export cell type assignments
adata.obs[['cell_type']].to_csv('cell_type_assignments.csv')