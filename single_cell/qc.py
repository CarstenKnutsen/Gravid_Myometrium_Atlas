'''Goal:Do QC on anndata object from 'create_count_table.py'
Date:230926
Author: Carsten Knutsen
conda_env:uterus_sc
'''


import pandas as pd
import os
import scanpy as sc
from anndata import AnnData
import seaborn as sns
import matplotlib.pylab as plt
import numpy as np
import scanpy.external as sce

figures = '/home/carsten/alvira_bioinformatics/uterus/data/figures/qc'
data = '/home/carsten/alvira_bioinformatics/uterus/data/single_cell_files/scanpy_files'
os.makedirs(figures, exist_ok = True)
os.makedirs(data, exist_ok = True)
sc.set_figure_params(dpi = 300, format = 'png')
sc.settings.figdir = figures

if __name__ == '__main__':
    adata = sc.read(f'{data}/uterus_all_cells_soupx.gz.h5ad')
    print(adata)
    sc.pp.calculate_qc_metrics(adata, qc_vars = ['mt'],expr_type='umis', percent_top=None, log1p=False, inplace=True)
    adata.obs['log10_total_umis'] = np.log10(adata.obs['total_umis'])
    adata.obs['log10_n_genes_by_umis'] = np.log10(adata.obs['n_genes_by_umis'])
    sc.pl.scatter(adata, x='log10_total_umis', y='log10_n_genes_by_umis', show=False, save='genes_by_counts_pretrim_log')
    sc.pl.highest_expr_genes(adata, n_top=20, show=False, save=f'_pretrim')
    sc.pl.violin(adata, ['log10_total_umis'],groupby='Patient',rotation=90,show=False, save = 'counts_pretrim')
    sc.pl.violin(adata, ['log10_n_genes_by_umis'],groupby='Patient', rotation=90,show=False,save = 'genes_pretrim')
    sc.pl.violin(adata, ['pct_umis_mt'],groupby='Patient', rotation=90,show=False,save = 'mt_pretrim')
    sc.pp.filter_cells(adata, min_counts=1000)
    sc.pp.filter_cells(adata, min_genes=500)
    print(adata)
    print(adata.var['seqname'].cat.categories)
    adata = adata[:,~adata.var['seqname'].isin(['chrM'])] #some mito genes not annotated with chromosome thus NA
    print(adata)
    adata = adata[:, ~adata.var['seqname'].isin(['chrM', 'NA'])] #some mito genes not annotated with chromosome thus NA    sc.pp.filter_genes(adata, min_cells=10)
    print(adata)
    sc.pl.scatter(adata, x='log10_total_umis', y='log10_n_genes_by_umis', show = False, save = 'genes_by_counts_posttrim_log')
    sc.pl.highest_expr_genes(adata, n_top=20, show = False, save = f'_posttrim')
    sc.pl.violin(adata, ['log10_total_umis'],rotation=90,show=False,groupby='Patient', save='counts_posttrim')
    sc.pl.violin(adata, ['log10_n_genes_by_umis'],rotation=90,show=False, groupby='Patient',save='genes_posttrim')
    # scrublet
    sce.pp.scrublet(adata, batch_key='Patient')
    sce.pl.scrublet_score_distribution(adata, show=False, save='scrublet_scores')
    print(adata)
    adata = adata[adata.obs['doublet_score'] < 0.1] # lowest recommended trim, results in around 10% bouvlet rate, in line with 10xs estimates. Lose umap blobs
    print(adata)
    adata.write(f'{data}/uterus_processed.gz.h5ad', compression='gzip')
