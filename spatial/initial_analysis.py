'''Goal:Perform initial qc and clustering '
Date:240823
Author: Carsten Knutsen
conda_env:squidpy
'''

import scanpy as sc
import scanpy.external as sce
import scanpy
import squidpy as sq
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import anndata
from pathlib import Path
import os
sc.logging.print_header()
print(f"squidpy=={sq.__version__}")
output = '/home/carsten/alvira_bioinformatics/uterus/data/vizgen/initial_analysis'
os.makedirs(output, exist_ok=True)
sc.set_figure_params(dpi=300, dpi_save=300,format="png")
sc._settings.ScanpyConfig(max_memory=31,n_jobs=6)
sc.settings.figdir = output
sample_dictionary = {
    # 'Sample_01':{'path':'/media/carsten/hdd/sequencing/uterus_spatial/SampleOutputs/202406191206_MyometriumSample1LY06-19-24_VMSC11602-merfishoutput/region_0/resegment_cellpose/',
    #                               'Treatement':'TNL'},
    #                  'Sample_02':{'path':'/media/carsten/hdd/sequencing/uterus_spatial/SampleOutputs/202406191335_MyometriumSample2-LY6-19-24_VMSC12602-merfishoutput/region_0/resegment_cellpose/',
    #                               'Treatement':'PTL'},
                     'Sample_01':{'path':'/media/carsten/hdd/sequencing/uterus_spatial/SampleOutputs/MyometriumSample3_DHResegmented/resegment_cellpose/',
                                  'Treatement':'TNL'},
                     'Sample_02':{'path':'/media/carsten/hdd/sequencing/uterus_spatial/SampleOutputs/MyometriumSample4_DHResegmented/resegment_cellpose/',
                                  'Treatement':'TNL'},
                     'Sample_03':{'path':'/media/carsten/hdd/sequencing/uterus_spatial/SampleOutputs/MyometriumSample5_DHResegmented/resegment_cellpose/',
                                  'Treatement':'TNL'},
                    }
sc_fn = adata_sc ='/home/carsten/alvira_bioinformatics/uterus/data/single_cell_files/scanpy_files/uterus_processed_celltyped.gz.h5ad'
gene_fn = '/home/carsten/alvira_bioinformatics/uterus/data/pilot/240325_merscope_planning/gene_list_trim.csv'
if __name__ == "__main__":
    adata_sc = sc.read(
        sc_fn)
    gene_df = pd.read_csv(
        gene_fn,
        index_col=0)
    gene_df_ct_markers = gene_df.loc[gene_df['Reason'].str.contains('ct_marker')].index.to_list()
    ct_marker_dict = {}
    for ct in adata_sc.obs['Cell Subtype'].cat.categories:
        if ct == 'Proliferative myeloid':
            continue
        ct_marker_dict[ct] = gene_df.loc[gene_df['Reason'].str.contains(ct)].index.to_list()
    lin_marker_dict = {}
    for ct in adata_sc.obs['Lineage'].cat.categories:
        lin_marker_dict[ct] = gene_df.loc[gene_df['Reason'].str.contains(ct)].index.to_list()
    # for sample, dt in sample_dictionary.items():
    #     spatial_data_path = dt['path']
    #     os.makedirs(f'{output}/{sample}', exist_ok=True)
    #     sc.settings.figdir = f'{output}/{sample}'
    #     adata = sq.read.vizgen(path=spatial_data_path, counts_file='cell_by_gene.csv', meta_file='cell_metadata.csv')
    #     sc.pp.calculate_qc_metrics(adata, percent_top=(1, 2, 5, 10, 20, 50, 100,), inplace=True)
    #     adata.uns['unassinged_transcripts'] = adata.obsm["blank_genes"].to_numpy().sum() / adata.var[
    #         "total_counts"].sum() * 100
    #     sc.pp.filter_cells(adata, min_counts=10)
    #     sc.pp.filter_cells(adata, min_genes=2)
    #     sc.pp.filter_genes(adata, min_counts=1)
    #     sc.pl.highest_expr_genes(adata, save=f'{sample}.png', show=False)
    #     adata.layers["counts"] = adata.X.copy()
    #     sc.pp.normalize_total(adata, inplace=True)
    #     sc.pp.log1p(adata, base=10)
    #     sc.pp.pca(adata)
    #     sc.pp.neighbors(adata)
    #     sc.tl.umap(adata)
    #     sc.tl.leiden(adata)
    #     sc.pl.umap(adata, color=['leiden'], save=f'leiden_{sample}.png', show=False)
    #     sc.tl.rank_genes_groups(adata, groupby='leiden')
    #     sc.pl.rank_genes_groups_dotplot(adata, save=f'leiden_markers_{sample}.png', show=False)
    #     sc.pl.umap(adata, color=adata.var.sort_values('total_counts', ascending=False).head(20).index,
    #                save=f'top20_genes_detected_{sample}.png', show=False)
    #     sc.pl.umap(adata, color=['n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts',
    #                              'pct_counts_in_top_1_genes', 'pct_counts_in_top_2_genes', 'pct_counts_in_top_5_genes',
    #                              'pct_counts_in_top_10_genes', 'pct_counts_in_top_50_genes',
    #                              'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes'], save=f'counts_qc_{sample}.png',
    #                show=False)
    #     sc.pl.violin(adata, ['pct_counts_in_top_1_genes', 'pct_counts_in_top_2_genes', 'pct_counts_in_top_5_genes',
    #                          'pct_counts_in_top_10_genes', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_50_genes',
    #                          'pct_counts_in_top_100_genes'], multi_panel=True, save=f'counts_pct_qc_{sample}.png', show=False)
    #
    #     adata.layers["log1p"] = adata.X.copy()
    #     ct_df = pd.DataFrame(index=adata.obs_names, columns=[f'{ct} score' for ct in ct_marker_dict.keys()])
    #     for ct, genes in ct_marker_dict.items():
    #         name = f'{ct} score'
    #         sc.tl.score_genes(adata, genes, score_name=name)
    #         ct_df[name] = adata.obs[name].copy()
    #     adata.obs['Cell Subtype_ind'] = ct_df.idxmax(axis=1)
    #     sc.pl.umap(adata, color='Cell Subtype_ind', save=f'Cell Subtype_ind_{sample}.png', show=False)
    #     sc.pl.umap(adata, color=[f'{ct} score' for ct in ct_marker_dict.keys()], save=f'celltype_marker_score_{sample}.png',
    #                show=False)
    #     for ct, genes in lin_marker_dict.items():
    #         name = f'{ct} score'
    #         sc.tl.score_genes(adata, genes, score_name=name)
    #     sc.pl.umap(adata, color=[f'{ct} score' for ct in lin_marker_dict.keys()], save=f'lineage_marker_score{sample}.png',
    #                show=False)
    #     os.makedirs(f'{output}/{sample}/gene_umaps', exist_ok=True)
    #     sc.settings.figdir = f'{output}/{sample}/gene_umaps'
    #     for gene in adata.var_names:
    #         sc.pl.umap(adata, color=gene, save=f'{gene}_{sample}', show=False)
    #     dt['adata'] = adata
    #     adata.write(f'{output}/{sample}/{sample}.gz.h5ad',compression='gzip')
    os.makedirs(f'{output}/combined', exist_ok=True)
    sc.settings.figdir = f'{output}/combined'
    # adatas = []
    # for sample, dt in sample_dictionary.items():
    #     spatial_data_path = dt['path']
    #     adata = sq.read.vizgen(path=spatial_data_path, counts_file='cell_by_gene.csv', meta_file='cell_metadata.csv')
    #     adata.obs['Sample'] = sample
    #     adatas.append(adata)
    # adata = anndata.concat(adatas)
    # sc.pp.calculate_qc_metrics(adata, percent_top=(1, 2, 5, 10, 20, 50, 100,), inplace=True)
    # adata.uns['unassinged_transcripts'] = adata.obsm["blank_genes"].to_numpy().sum() / adata.var[
    #     "total_counts"].sum() * 100
    # sc.pp.filter_cells(adata, min_counts=10)
    # sc.pp.filter_cells(adata, min_genes=2)
    # sc.pp.filter_genes(adata, min_counts=1)
    # sc.pl.highest_expr_genes(adata, save=True, show=False)
    # adata.layers["counts"] = adata.X.copy()
    # sc.pp.normalize_total(adata, inplace=True)
    # sc.pp.log1p(adata)
    # sce.pp.magic(adata, knn=10, random_state=7)
    # sc.pp.pca(adata)
    # adata.write(f'{output}/all_runs.gz.h5ad',compression='gzip')
    # sc.pl.pca_variance_ratio(adata,show=False, save=True)
    # sce.pp.harmony_integrate(adata,key='Sample')
    # sc.pp.neighbors(adata, use_rep='X_pca_harmony')
    # adata.write(f'{output}/all_runs.gz.h5ad',compression='gzip')
    # adata.obs['Sample number'] = adata.obs['Sample'].str.split('_').str[1]
    # print(adata)
    # sc.tl.umap(adata)
    # sc.tl.leiden(adata, flavor="igraph",n_iterations=2)
    # adata.write(f'{output}/all_runs.gz.h5ad', compression='gzip')
    adata = sc.read(f'{output}/all_runs.gz.h5ad')
    sc.tl.leiden(adata, resolution=1.5)
    leiden_ct_dictionary = {'0': 'Uterine smooth muscle', '1': 'Uterine smooth muscle', '2': 'Uterine smooth muscle',
                            '3': 'Uterine smooth muscle', '4': 'Uterine smooth muscle', '5': 'Uterine smooth muscle',
                            '6': 'Uterine smooth muscle', '7': 'Fibroblast', '8': 'Uterine smooth muscle',
                            '9': 'Fibroblast', '10': 'Uterine smooth muscle', '11': 'Myeloid',
                            '12': 'Vascular smooth muscle',
                            '13': 'Myeloid', '14': 'Capillary', '15': 'Uterine smooth muscle', '16': 'Fibroblast',
                            '17': 'Capillary', '18': 'Lymphatic EC', '19': 'Uterine smooth muscle', '20': 'Myeloid',
                            '21': 'Uterine smooth muscle', '22': 'Fibroblast', '23': 'Uterine smooth muscle',
                            '24': 'Trophoblast',
                            '25': 'Myeloid', '26': 'Uterine smooth muscle', '27': 'Lymphoid', '28': 'Capillary',
                            '29': 'Uterine smooth muscle', '30': 'Macrovascular', '31': 'Endometrial'}
    adata.obs["Cell Identity"] = (
        adata.obs[f"leiden"].map(leiden_ct_dictionary).astype("str")
    )

    sc.pl.umap(adata, color=['Sample', 'leiden','Cell Identity'], save='leiden.png', show=False)
    sc.tl.rank_genes_groups(adata, groupby='leiden')
    # sc.pl.rank_genes_groups_dotplot(adata, dendrogram=False,save='leiden_markers.png', show=False)
    sc.pl.umap(adata, color=adata.var.sort_values('total_counts', ascending=False).head(20).index,
               save='top20_genes_detected.png', show=False)
    sc.pl.umap(adata, color=['n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts',
                             'pct_counts_in_top_1_genes', 'pct_counts_in_top_2_genes', 'pct_counts_in_top_5_genes',
                             'pct_counts_in_top_10_genes', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_50_genes',
                             'pct_counts_in_top_100_genes'], save='counts_qc.png', show=False)
    sc.pl.violin(adata, ['pct_counts_in_top_1_genes', 'pct_counts_in_top_2_genes', 'pct_counts_in_top_5_genes',
                         'pct_counts_in_top_10_genes', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_50_genes',
                         'pct_counts_in_top_100_genes'], groupby='Sample number',multi_panel=True, save='pct_counts_qc.png', show=False)
    sc.pl.violin(adata, ['n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', ],groupby='Sample number',
                 multi_panel=True, save='counts_qc.png', show=False)
    ct_df = pd.DataFrame(index=adata.obs_names, columns=[f'{ct}' for ct in ct_marker_dict.keys()])
    for ct, genes in ct_marker_dict.items():
        print(ct)
        genes_short = [x for x in genes if x in adata.var_names]
        print(genes_short)
        name = f'{ct} score'
        sc.tl.score_genes(adata, genes_short, score_name=name)
        ct_df[ct] = adata.obs[name].copy()
    adata.obs['Cell Subtype_ind'] = ct_df.idxmax(axis=1)
    sc.pl.umap(adata, color='Cell Subtype_ind', save='Cell Subtype_ind.png', show=False)
    sc.pl.umap(adata, color=[f'{ct} score' for ct in ct_marker_dict.keys()], save='celltype_marker_score.png',
               show=False)
    for ct, genes in lin_marker_dict.items():
        name = f'{ct} score'
        sc.tl.score_genes(adata, genes, score_name=name)
    sc.pl.umap(adata, color=[f'{ct} score' for ct in lin_marker_dict.keys()], save='lineage_marker_score.png',
               show=False)
    os.makedirs(f'{output}/combined/gene_umaps', exist_ok=True)
    sc.settings.figdir = f'{output}/combined/gene_umaps'
    for gene in adata.var_names:
        sc.pl.umap(adata, color=gene, save=gene, show=False)

    sq.gr.spatial_neighbors(adata)
    sq.gr.nhood_enrichment(adata, cluster_key='Cell Identity')
    adata.write(f'{output}/all_runs.gz.h5ad',compression='gzip')
    df = pd.DataFrame(data=adata.uns['Cell Identity_nhood_enrichment']['zscore'],
                      index=adata.obs['Cell Identity'].cat.categories,
                      columns=adata.obs['Cell Identity'].cat.categories)
    sns.set_style('white')
    fig = sns.clustermap(df, cmap='vlag', vmin=-50, vmax=50)
    plt.tight_layout()
    fig.savefig(f'{output}/neighborhood_identity.png', dpi=300)

    sq.gr.nhood_enrichment(adata, cluster_key='Cell Subtype_ind')
    df = pd.DataFrame(data=adata.uns['Cell Subtype_ind_nhood_enrichment']['zscore'],
                      index=adata.obs['Cell Subtype_ind'].cat.categories,
                      columns=adata.obs['Cell Subtype_ind'].cat.categories)
    sns.set_style('white')
    fig = sns.clustermap(df, cmap='vlag', vmin=-50, vmax=50)
    plt.tight_layout()
    fig.savefig(f'{output}/neighborhood_Cell Subtype_ind.png', dpi=300)





