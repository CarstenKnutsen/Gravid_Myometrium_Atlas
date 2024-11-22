'''Goal:Create count table for all cells after SoupX cleanup of uterus samples
Date:230926
Author: Carsten Knutsen
conda_env:uterus_sc
'''


import os
import scanpy as sc
import numpy as np
import pandas as pd
import string
import anndata
from gtfparse import read_gtf
from anndata import AnnData
from collections import defaultdict

# gtf_fn = '/media/carsten/hdd/genomes/human/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
# data = '/home/carsten/alvira_bioinformatics/uterus/data/single_cell_files/cellranger_output'
# output = '/home/carsten/alvira_bioinformatics/uterus/data/single_cell_files/scanpy_files'
# os.makedirs(output, exist_ok=True)
# if __name__ == '__main__':
#     runs = os.listdir(data)
#     adatas = []
#     gtf = read_gtf(gtf_fn)
#     print(gtf)
#     print(gtf['gene_id'].str.split('.').str[0])
#     gtf['gene_id'] = gtf['gene_id'].str.split('.').str[0]
#     gene_name_dict = pd.Series(gtf.gene_name.values, index=gtf.gene_id).to_dict()
#     for x in gene_name_dict.keys():
#         if gene_name_dict[x] == '':
#             gene_name_dict[x] = x
#     for run in runs:
#         print(run)
#         folder = f'{data}/{run}'
#         adata = sc.read_10x_h5(f'{folder}/outs/filtered_feature_bc_matrix.h5')
#         adata.var_names_make_unique()
#         adata.obs_names = run + '_' + adata.obs_names
#         adata.obs['Patient'] = run.split('_')[-2]
#         adata.obs['Group'] = run.split('_')[0]
#         adata.obs['Contractility'] =[run.split('_')[1] if run.startswith('TN') else None][0]
#         adata.obs['Term'] = ['Term' if run.startswith('T') else 'Preterm'][0]
#         adata.obs['Labor'] = ['Nonlaboring' if run[1]=='N' else 'Laboring'][0]
#         adatas.append(adata.copy())
#     adata = anndata.concat(adatas)
#     for column in ['gene_id', 'gene_type', 'seqname', 'transcript_name', 'protein_id']:
#         temp_dict = pd.Series(gtf[column].values, index=gtf['gene_name']).to_dict()
#         temp_dict.update(pd.Series(gtf[column].values, index=gtf['gene_id']).to_dict())
#         temp_dict = defaultdict(lambda: None, temp_dict)
#         adata.var[column] = [temp_dict[x] for x in adata.var.index]
#     adata.var['mt'] = [True if x == 'chrM' else False for x in adata.var['seqname']]
#     adata.write(f'{output}/uterus_all_cells_soupx.gz.h5ad', compression='gzip')

















gtf_fn = '/media/carsten/hdd/genomes/human/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
data = '/home/carsten/alvira_bioinformatics/uterus/data/single_cell_files/soupx'
output = '/home/carsten/alvira_bioinformatics/uterus/data/single_cell_files/scanpy_files'
os.makedirs(output, exist_ok=True)
def read_adata(folder):
    adata = sc.read_mtx(f'{folder}/matrix.mtx').T
    features = pd.read_csv(f'{folder}/genes.tsv',
                          sep = '\t',
                          header = None)
    bc = pd.read_csv(f'{folder}/barcodes.tsv',
                          sep = '\t',
                    header = None)
    features.rename(columns={0:'gene_id',
                            1: 'gene_symbol',
                            2: 'category'},
                   inplace = True)

    adata.var = features
    adata.obs_names = bc[0]
    adata.var_names = adata.var['gene_id'].values
    return adata
if __name__ == '__main__':
    runs = os.listdir(data)
    adatas = []
    gtf = read_gtf(gtf_fn)
    print(gtf.columns)
    gtf['gene_id'] = gtf['gene_id'].str.split('.').str[0]
    gene_name_dict = pd.Series(gtf.gene_name.values, index=gtf.gene_id).to_dict()
    for x in gene_name_dict.keys():
        if gene_name_dict[x] == '':
            gene_name_dict[x] = x
    for run in runs:
        print(run)
        folder = f'{data}/{run}'
        adata = read_adata(folder)
        adata.var_names = [gene_name_dict[x] for x in adata.var_names]
        adata.var_names_make_unique()
        adata.obs_names = run + '_' + adata.obs_names
        adata.obs['Patient'] = run.split('_')[-2]
        adata.obs['Group'] = run.split('_')[0]
        adata.obs['Contractility'] =[run.split('_')[1] if run.startswith('TN') else 'ND'][0]
        adata.obs['GroupContract'] = adata.obs['Group'].astype('str') + '-' + adata.obs['Contractility'].astype('str')
        adata.obs['Term'] = ['Term' if run.startswith('T') else 'Preterm'][0]
        adata.obs['Labor'] = ['Nonlaboring' if run[1]=='N' else 'Laboring'][0]
        adatas.append(adata.copy())
    adata = anndata.concat(adatas)
    for column in ['gene_id', 'gene_type', 'seqname', 'transcript_name', 'protein_id']:
        temp_dict = pd.Series(gtf[column].values, index=gtf['gene_name']).to_dict()
        temp_dict.update(pd.Series(gtf[column].values, index=gtf['gene_id']).to_dict())
        temp_dict = defaultdict(lambda: None, temp_dict)
        adata.var[column] = [temp_dict[x] for x in adata.var.index]
    adata.var['mt'] = [True if x == 'chrM' else False for x in adata.var['seqname']]
    adata.write(f'{output}/uterus_all_cells_soupx.gz.h5ad', compression='gzip')