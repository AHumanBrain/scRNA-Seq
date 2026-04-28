#!/usr/bin/env python3
"""
scRNA-Seq Integrated Analysis Pipeline
=======================================
Accepts one or more 10x MEX matrix directories as input.
- Single input: runs standard Scanpy analysis (no Harmony)
- Multiple inputs: merges datasets and applies Harmony batch correction

Features: Scrublet doublet detection, multi-resolution Leiden clustering,
CellTypist automated cell annotation, and Harmony batch integration.
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import sys
import os
import glob
import scrublet as scr
import celltypist
from celltypist import models

sc.settings.verbosity = 3
sc.settings.figdir = './figures/'
os.makedirs('./figures/', exist_ok=True)

# =========================================================================
# INPUT HANDLING — Auto-detect MEX directories from arguments
# =========================================================================
if len(sys.argv) < 2:
    print("Error: At least one input directory required.")
    print("Usage: analyze.py <matrix_dir_1> [matrix_dir_2] [matrix_dir_3] ...")
    sys.exit(1)

input_dirs = sys.argv[1:]
print(f"Received {len(input_dirs)} input director(ies): {input_dirs}")


def find_mex_dir(base_path):
    """Recursively search for a MEX matrix directory (contains matrix.mtx*)."""
    # Check if the directory itself contains the matrix
    if glob.glob(os.path.join(base_path, 'matrix.mtx*')):
        return base_path
    # Search subdirectories (handles 10x nested structures like filtered_gene_bc_matrices/hg19/)
    for root, dirs, files in os.walk(base_path):
        for f in files:
            if f.startswith('matrix.mtx'):
                return root
    return None


def detect_var_names(mex_dir):
    """Detect whether the MEX directory uses 'genes.tsv' or 'features.tsv'."""
    if os.path.exists(os.path.join(mex_dir, 'features.tsv.gz')) or \
       os.path.exists(os.path.join(mex_dir, 'features.tsv')):
        return 'gene_symbols'
    return 'gene_symbols'


# Load all datasets
datasets = []
dataset_labels = []

for input_dir in input_dirs:
    mex_dir = find_mex_dir(input_dir)
    if mex_dir is None:
        print(f"Warning: No MEX matrix found in {input_dir}, skipping.")
        continue

    label = os.path.basename(os.path.dirname(input_dir)) or os.path.basename(input_dir)
    # Clean up label
    label = label.replace('filtered_gene_bc_matrices', '').replace('filtered_feature_bc_matrix', '').strip('/')
    if not label or label == '.':
        label = f"sample_{len(datasets)+1}"

    print(f"Loading dataset '{label}' from {mex_dir}...")
    adata = sc.read_10x_mtx(mex_dir, var_names='gene_symbols', cache=False)
    adata.var_names_make_unique()
    adata.obs['dataset'] = label
    datasets.append(adata)
    dataset_labels.append(label)
    print(f"  → {adata.n_obs} cells, {adata.n_vars} genes")

if len(datasets) == 0:
    print("Error: No valid datasets found. Exiting.")
    sys.exit(1)

multi_sample = len(datasets) > 1

# =========================================================================
# MERGE DATASETS (if multiple inputs)
# =========================================================================
if multi_sample:
    # Find common genes across all datasets
    common_genes = datasets[0].var_names
    for ds in datasets[1:]:
        common_genes = common_genes.intersection(ds.var_names)
    print(f"Common genes across {len(datasets)} datasets: {len(common_genes)}")

    # Subset to common genes
    for i in range(len(datasets)):
        datasets[i] = datasets[i][:, common_genes].copy()

    adata = ad.concat(datasets, label='dataset_key', keys=dataset_labels)
    print(f"Concatenated dataset shape: {adata.shape}")
else:
    adata = datasets[0]
    adata.obs['dataset_key'] = dataset_labels[0]
    print(f"Single dataset mode: {adata.shape}")

# =========================================================================
# DOUBLET DETECTION (Scrublet)
# =========================================================================
print("Running Scrublet doublet detection...")
sc.external.pp.scrublet(adata, batch_key='dataset_key')
sc.pl.scrublet_score_distribution(adata, save='_scrublet.png')
adata = adata[~adata.obs.predicted_doublet, :]
print(f"Cells after removing doublets: {adata.n_obs}")

# =========================================================================
# QUALITY CONTROL
# =========================================================================
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, save='_qc_integrated.png')

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
print(f"Cells after QC: {adata.n_obs}")

# =========================================================================
# NORMALIZATION
# =========================================================================
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Save raw (unscaled) for Celltypist
adata.raw = adata

# =========================================================================
# FEATURE SELECTION
# =========================================================================
if multi_sample:
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key='dataset_key')
else:
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata, save='_integrated.png')
adata = adata[:, adata.var.highly_variable]

# =========================================================================
# DIMENSIONALITY REDUCTION
# =========================================================================
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

# =========================================================================
# BATCH CORRECTION (Harmony — only for multi-sample)
# =========================================================================
if multi_sample:
    print("Running Harmony Integration...")
    sc.external.pp.harmony_integrate(adata, 'dataset_key')
    use_rep = 'X_pca_harmony'
else:
    print("Single sample — skipping Harmony integration.")
    use_rep = 'X_pca'

# =========================================================================
# NEIGHBORHOOD & CLUSTERING
# =========================================================================
print(f"Computing neighbors using {use_rep}...")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40, use_rep=use_rep)
sc.tl.umap(adata)

if multi_sample:
    sc.pl.umap(adata, color='dataset_key', save='_integration_batch.png')

# Multi-resolution Leiden clustering
resolutions = [0.3, 0.5, 0.8]
for res in resolutions:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')
sc.pl.umap(adata, color=[f'leiden_{r}' for r in resolutions], save='_leiden_resolutions.png')

# =========================================================================
# AUTOMATED ANNOTATION (CellTypist)
# =========================================================================
print("Running Celltypist for automated labeling...")
models.download_models(force_update=False, model=['Immune_All_Low.pkl'])
model = models.Model.load(model='Immune_All_Low.pkl')
predictions = celltypist.annotate(adata.raw.to_adata(), model=model, majority_voting=True)
adata.obs['celltypist_prediction'] = predictions.predicted_labels['predicted_labels']
adata.obs['celltypist_majority_voting'] = predictions.predicted_labels['majority_voting']

sc.pl.umap(adata, color='celltypist_majority_voting', legend_loc='on data',
           title='Celltypist Auto Annotation', frameon=False, save='_celltypist.png')

# =========================================================================
# OUTPUTS
# =========================================================================
adata.write('pbmc_integrated.h5ad')

export_cols = ['dataset_key', 'leiden_0.5', 'celltypist_majority_voting']
if 'predicted_doublet' in adata.obs.columns:
    export_cols.append('predicted_doublet')
adata.obs[export_cols].to_csv('cell_type_assignments.csv')

print("Analysis complete.")

# =========================================================================
# EXPERIMENTAL: RNA Velocity (scVelo) - COMMENTED OUT
# =========================================================================
# RNA Velocity is highly resource-intensive and requires `spliced` and
# `unspliced` count matrices which are NOT generated by standard CellRanger
# 10x output (they require re-mapping with velocyto or kb-python).
# Due to the 12GB RAM hardware limits, this block remains commented as it
# would cause Out-Of-Memory exhaustion or fail due to missing .loom data.
#
# import scvelo as scv
# print("Running RNA Velocity...")
# scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
# scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
# scv.tl.velocity(adata)
# scv.tl.velocity_graph(adata)
# scv.pl.velocity_embedding_stream(adata, basis='umap', save='_scvelo.png')
