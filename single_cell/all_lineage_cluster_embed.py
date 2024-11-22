'''Goal:Create tissue level embedding and classify cells by lineage'
Date:230926
Author: Carsten Knutsen
conda_env:uterus_sc
'''

import pandas as pd
import os
import scanpy as sc
import leidenalg as la
import itertools
import scanpy.external as sce

figures = '/home/carsten/alvira_bioinformatics/uterus/data/figures/tissue_embedding/all_lineages'
data = '/home/carsten/alvira_bioinformatics/uterus/data/single_cell_files/scanpy_files'
os.makedirs(figures, exist_ok = True)
sc.set_figure_params(dpi = 300, dpi_save = 300, format = 'png')
sc.settings.figdir = figures
if __name__ == '__main__':
    adata = sc.read(f'{data}/uterus_processed.gz.h5ad')
    print(adata)
    adata.layers['soupx'] = adata.X.copy()
    sc.pp.normalize_total(adata, key_added=None, target_sum=1e4)
    adata.layers['cp10k'] = adata.X.copy()
    sc.pp.log1p(adata)
    adata.layers['log'] = adata.X.copy()
    sc.pp.highly_variable_genes(adata,
                                # n_top_genes = 2000,
                                batch_key='Patient'
                                )
    sc.pp.pca(adata, use_highly_variable=True)
    sce.pp.harmony_integrate(adata, 'Patient',adjusted_basis='X_pca')
    sc.pp.neighbors(adata, use_rep='X_pca')

    sc.tl.leiden(adata, key_added='leiden')
    sc.tl.umap(adata, min_dist=0.3)
    sc.pl.umap(adata,color='leiden',show=False,save = 'leiden_original.png')
    sc.pl.umap(adata,color='doublet_score',show=False,save = 'doublet_original.png')

    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups_dotplot(adata, groupby='leiden',
                                    n_genes=int(150 / len(adata.obs['leiden'].unique())), show=False,
                                    save=f'leiden_original_markers.png')

    adata = adata[~adata.obs['leiden'].isin(['5','18'])] # 5 looks to be low quality, 18 RBC

    sc.pp.highly_variable_genes(adata,
                                # n_top_genes = 2000,
                                batch_key='Patient'
                                )
    sc.pp.pca(adata, use_highly_variable=True)
    sce.pp.harmony_integrate(adata, 'Patient', adjusted_basis='X_pca')
    sc.pp.neighbors(adata, use_rep='X_pca')

    sc.tl.leiden(adata, key_added='leiden')
    sc.tl.umap(adata, min_dist=0.3)
    sc.pl.umap(adata, color='leiden', show=False, save='leiden_original2.png')
    sc.tl.dendrogram(adata,'leiden')
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups_dotplot(adata, groupby='leiden',
                                    n_genes=int(150 / len(adata.obs['leiden'].unique())), show=False,
                                    save=f'leiden_original2_markers.png')

    adata = adata[~adata.obs['leiden'].isin(['13'])]  # 13 looks to be low quality

    sce.pp.harmony_integrate(adata, 'Patient', adjusted_basis='X_pca')
    sc.pp.neighbors(adata, use_rep='X_pca')

    sc.tl.leiden(adata, key_added='leiden')
    sc.tl.umap(adata, min_dist=0.3)
    sc.pl.umap(adata, color='leiden', show=False, save='leiden.png')
    sc.tl.dendrogram(adata, 'leiden')
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups_dotplot(adata, groupby='leiden',
                                    n_genes=int(150 / len(adata.obs['leiden'].unique())), show=False,
                                    save=f'leiden_markers.png')
    with pd.ExcelWriter(
            f"{figures}/leiden_markers.xlsx", engine="xlsxwriter"
    ) as writer:
        for ct in adata.obs["leiden"].cat.categories:
            df = sc.get.rank_genes_groups_df(
                adata, key="rank_genes_groups", group=ct
            )
            df.index = df['names']
            df.to_excel(writer, sheet_name=f"{ct} v rest"[:31])
    sc.pl.pca_overview(adata, color='leiden', show=False, save = True)
    sc.pl.pca_variance_ratio(adata, show = False, save ='variance_ratio')
    sc.pl.pca_loadings(adata, components = ','.join([str(x) for x in range(1,10)]), show = False, save = True)
    ## Assign lineage from gene expression
    genes = ['COL1A1', 'PECAM1', 'PTPRC', 'PAEP', 'leiden']
    genedf = sc.get.obs_df(adata, keys=genes)
    grouped = genedf.groupby("leiden")
    mean = grouped.mean()
    mean_t = mean.T
    mean_t.to_csv(f'{figures}/uterus_processed.csv')
    lineage_dict = {}
    for cluster in mean_t.columns:
        gene = mean_t[cluster].idxmax()
        if gene == 'PECAM1':
            lineage_dict[cluster] = 'Endothelial'
        elif gene == 'PTPRC':
            lineage_dict[cluster] = 'Immune'
        elif gene == 'PAEP':
            lineage_dict[cluster] = 'Endometrial'
        elif gene == 'COL1A1':
            lineage_dict[cluster] = 'Mesenchymal'
    adata.obs['Lineage'] = [lineage_dict[x] for x in adata.obs['leiden']]
    adata.write(f'{data}/uterus_assigned_lineage.gz.h5ad', compression='gzip')
    print(adata)
    umaps = f'{figures}/umaps/genes'
    os.makedirs(umaps,exist_ok=True)
    sc.settings.figdir = umaps
    for ct in adata.obs["leiden"].cat.categories:
        df = sc.get.rank_genes_groups_df(
            adata, key="rank_genes_groups", group=ct
        )
        genes = df['names'].values[:5]
        for gene in genes:
            sc.pl.umap(adata, color=gene, title=f'{gene} \n leiden {ct} marker',alpha=0.5, show=False, save=f'_{ct}_{gene}')
    for gene in ['COL1A1',
                 'EPCAM',
                 'PTPRC',
                 'PECAM1',
                 'CDH5',
        'PAEP',
                 'ACTA2',
                 'RYR3',
                 ]:
        sc.pl.umap(adata, color=gene, alpha=0.5, show=False, save=f'_{gene}')
    metadata = f'{figures}/umaps/metadata'
    os.makedirs(metadata, exist_ok=True)
    sc.settings.figdir = metadata
    metadatas = ['log10_n_genes_by_umis',
                 'Patient',
                 'Group',
                 'GroupContract',
                 'Contractility',
                 'Term',
                 'doublet_score',
                 'Labor',
        'leiden',
                 'pct_umis_mt',
        'Lineage']
    for md in metadatas:
        if md == 'leiden':
            sc.pl.umap(adata, legend_loc='on data', color=md, alpha=0.5, show = False,save=f'_{md}')
        else:
            sc.pl.umap(adata, color=md, title=md,alpha=0.5, show=False, save=f'_{md}')
    pd.DataFrame(index=adata.obs_names, data = adata.obsm['X_umap']).to_csv(f'{figures}/allcells_umap.csv')
    with pd.ExcelWriter(f'{figures}/metadata_counts.xlsx', engine='xlsxwriter') as writer:
        obs_list = ['leiden', 'Patient','GroupContract']
        num_obs = len(obs_list) + 1
        for ind in range(0, num_obs):
            for subset in itertools.combinations(obs_list, ind):
                if len(subset) != 0:
                    subset = list(subset)
                    if len(subset) == 1:
                        key = subset[0]
                        adata.obs[key].value_counts().to_excel(writer, sheet_name=key)
                    else:
                        key = "_".join(subset)
                        adata.obs.groupby(subset[:-1])[subset[-1]].value_counts().to_excel(writer, sheet_name=key[:31])

