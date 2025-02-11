'''Goal:Create lineage level embedding and classify cells by cell Subtype'
Date:230926
Author: Carsten Knutsen
conda_env:uterus_sc
'''

import pandas as pd
import os
import scanpy as sc
import scanpy.external as sce
import seaborn as sns
import numpy as np
import matplotlib.pylab as plt
import leidenalg as la
import northstar
import muon as mu
import itertools
import anndata

figures = '/home/carsten/alvira_bioinformatics/uterus/data/figures/tissue_embedding'
data = '/home/carsten/alvira_bioinformatics/uterus/data/single_cell_files/scanpy_files'
os.makedirs(figures, exist_ok=True)
sc.set_figure_params(dpi=300, dpi_save=300, format="png")
sc.settings.figdir = figures

if __name__ == "__main__":
    # read in processed data
    adata = sc.read(f'{data}/uterus_assigned_lineage.gz.h5ad')
    adata = adata[adata.obs['Patient'] != 'P28']
    print(adata)
    # create empty cell type columns for both northstar and cluster identities
    adata.obs["Cell Subtype"] = pd.Series(index=adata.obs.index, data=None, dtype="str")
    # do this for every Lineage
    for Lineage in adata.obs["Lineage"].unique():
        print(Lineage)
        figures_lin = f"{figures}/{Lineage}"
        print(figures_lin)
        os.makedirs(figures_lin, exist_ok=True)
        sc.settings.figdir = figures_lin
        lin_adata = anndata.AnnData(
            adata[adata.obs["Lineage"] == Lineage].layers['log'].copy(),
            obs=adata[adata.obs["Lineage"] == Lineage].obs,
            var=adata[adata.obs["Lineage"] == Lineage].var,
        )
        del lin_adata.var["highly_variable"]
        print(lin_adata)
        if Lineage == 'Endometrial':

            sc.pp.highly_variable_genes(lin_adata,
                                        )
        else:
            sc.pp.highly_variable_genes(lin_adata,
                                        batch_key="Patient"
                                        )
        sc.pl.highly_variable_genes(lin_adata, save="highly_variable", show=False)
        sc.pp.pca(lin_adata, use_highly_variable=True)
        sce.pp.harmony_integrate(lin_adata, key="Patient", max_iter_harmony=20)
        sc.pp.neighbors(lin_adata, use_rep='X_pca_harmony')
        sc.tl.leiden(
            lin_adata,
            key_added=f"leiden_{Lineage}",
        )
        # # leiden cluster identity based off of cell type markers
        leiden_ct_dict = {
            "Mesenchymal": {
                "0":"Uterine smooth muscle",
                "1": "Vessel fibroblast",
                "2": "Vessel fibroblast",
                "3": "Uterine smooth muscle",
                "4": "Vascular smooth muscle",
                "5": "Vessel fibroblast",
                "6": "Trophoblast",
                "7": "Matrix fibroblast",
                "8": "low-quality",
                "9": "low-quality",
                "10": "Vessel fibroblast",
                "11": "low-quality",
                "12": "low-quality",

            },
            "Endometrial": {
                "0": "Glandular",
                "1": "Lumenal",
                "2": "Glandular",
                "3": "Lumenal",
                "4": "Lumenal",
                "5": "Glandular",
                "6": "Glandular",
                "7": "low-quality",
                "8": "Glandular",
                "9": "Lumenal",
                "10": "Ciliated",
                "11": "Epithelial",
                "12": "low-quality",
                "13": "low-quality",
                "14": "low-quality",
                "15": "Glandular",

            },
            "Endothelial": {
                "0": "Lymphatic EC",
                "1": "Capillary",
                "2": "Lymphatic EC",
                "3": "Capillary",
                "4": "Capillary",
                "5": "Venous EC",
                "6": "Arterial EC",
                "7": "Capillary",
                "8": "low-quality",
                "9": "low-quality",
                "10": "Lymphatic EC",
                "11": "Lymphatic EC",

            },
            "Immune": {
                # "0": "Macrophage",
                # "1": "Macrophage",
                # "2": "T cell_1",
                # "3": "Monocyte_1",
                # "4": "Macrophage",
                # "5": "Macrophage",
                # "6": "Macrophage",
                # "7": "Macrophage",
                # "8": "T cell_2",
                # "9": "Monocyte_2",
                # "10": "Monocyte_3",
                # "11": "T cell_3",
                # "12": "Macrophage",
                # "13": "Mast cell",
                # "14": "B cell",
                # "15": "Dendritic cell",
                "0": "Myeloid",
                "1": "T cell",
                "2": "Myeloid",
                "3": "Myeloid",
                "4": "Myeloid",
                "5": "Myeloid",
                "6": "Myeloid",
                "7": "NK cell",
                "8": "Myeloid",
                "9": "Myeloid",
                "10": "low-quality",
                "11": "Dendritic",
                "12": "Myeloid",
                "13": "Basophil",
                "14": "Proliferative myeloid",
                "15": "Plasmacytoid DC",

            }
        }
        lin_adata.obs["Cell Subtype"] = (
            lin_adata.obs[f"leiden_{Lineage}"].map(leiden_ct_dict[Lineage]).astype("str")
        )
        sc.tl.dendrogram(lin_adata, groupby=f"leiden_{Lineage}")
        sc.tl.rank_genes_groups(lin_adata, f"leiden_{Lineage}", method="wilcoxon")
        print(lin_adata.obs[f"leiden_{Lineage}"].cat.categories)
        sc.pl.rank_genes_groups_dotplot(
            lin_adata,
            groupby=f"leiden_{Lineage}",
            n_genes=int(150 / len(lin_adata.obs[f"leiden_{Lineage}"].unique())),
            show=False,
            save=f"{Lineage}_leiden_markers.png",
        )
        # sc.pl.DotPlot(lin_adata, genes, groupby=f"leiden_{Lineage}",).style(
        #     cmap="Reds"
        # ).add_totals().savefig(
        #     os.path.join(figures_lin, f"dotplot_{Lineage}_leiden.png"), dpi=300
        # )
        # plt.clf()
        lin_adata = lin_adata[~lin_adata.obs['Cell Subtype'].isin(['low-quality'])]
        adata.obs["Cell Subtype"].loc[lin_adata.obs.index] = lin_adata.obs[
            "Cell Subtype"
        ]
        # Visualize cell type connections, generate umap for Lineage, save Lineage umap to main object
        sc.tl.paga(lin_adata, groups="Cell Subtype")
        sc.pl.paga(
            lin_adata,
            color="Cell Subtype",
            show=False,
            save=f"_{Lineage}_Cell Subtype.png",
        )
        sc.tl.umap(lin_adata, min_dist=0.05)
        # Add Lineage umaps and leiden clusters to top level
        adata.obs[f"umap_{Lineage}_1"] = np.nan
        adata.obs[f"umap_{Lineage}_2"] = np.nan
        lin_adata.obs[f"umap_{Lineage}_1"] = [x[0] for x in lin_adata.obsm["X_umap"]]
        lin_adata.obs[f"umap_{Lineage}_2"] = [x[1] for x in lin_adata.obsm["X_umap"]]
        adata.obs[f"umap_{Lineage}_1"].loc[lin_adata.obs.index] = lin_adata.obs[
            f"umap_{Lineage}_1"
        ]
        adata.obs[f"umap_{Lineage}_2"].loc[lin_adata.obs.index] = lin_adata.obs[
            f"umap_{Lineage}_2"
        ]
        adata.obs[f"leiden_{Lineage}"] = np.nan
        adata.obs[f"leiden_{Lineage}"].loc[lin_adata.obs.index] = lin_adata.obs[
            f"leiden_{Lineage}"
        ]
        adata.obsm[f"X_umap_{Lineage}"] = adata.obs[
            [f"umap_{Lineage}_1", f"umap_{Lineage}_2"]
        ].to_numpy()
        del adata.obs[f"umap_{Lineage}_1"]
        del adata.obs[f"umap_{Lineage}_2"]
        sc.pl.pca_overview(lin_adata, color=f"leiden_{Lineage}", show=False, save=True)
        sc.pl.pca_variance_ratio(lin_adata, show=False, save="variance_ratio")
        sc.pl.pca_loadings(
            lin_adata, components="1,2,3,4,5", show=False, save="loadings"
        )
        plt.close()
        plt.clf()
        metadata = f"{figures_lin}/metadata"
        os.makedirs(metadata, exist_ok=True)
        sc.settings.figdir = metadata
        for color in [
            'log10_n_genes_by_umis',
            'Cell Subtype',
            'Patient',
            'Group',
            'Contractility',
            'Term',
            'doublet_score',
            'Labor',
            'leiden',
            f"leiden_{Lineage}",
            'pct_umis_mt',
            'Lineage'
        ]:
            if color in ["leiden",f"leiden_{Lineage}"]:
                sc.pl.umap(
                    lin_adata,
                    color=color,
                    alpha=0.5,
                    legend_loc="on data",
                    show=False,
                    save=f"{color}",
                )
            else:
                sc.pl.umap(
                    lin_adata, color=color, alpha=0.5, show=False, save=f"{color}"
                )
        gene_fn = f"{figures_lin}/genes"
        os.makedirs(gene_fn, exist_ok=True)
        sc.settings.figdir = gene_fn
        with pd.ExcelWriter(
                f"{figures_lin}/{Lineage}_leiden_markers.xlsx", engine="xlsxwriter"
        ) as writer:
            for ct in lin_adata.obs[f"leiden_{Lineage}"].cat.categories:
                df = sc.get.rank_genes_groups_df(
                    lin_adata, key="rank_genes_groups", group=ct
                )
                df.to_excel(writer, sheet_name=f"{ct} v rest"[:31])
                genes = df['names'].values[:5]
                for gene in genes:
                    sc.pl.umap(lin_adata, color=gene, title=f'{gene} \n leiden {ct} marker', alpha=0.5, show=False,
                               save=f'_{ct}_{gene}')
        sc.settings.figdir = figures_lin
        # sc.pl.DotPlot(lin_adata, genes, groupby="Cell Subtype",).add_dendrogram(
        #     show=True
        # ).style(cmap="Reds").add_totals().savefig(
        #     os.path.join(figures_lin, f"dotplot_{Lineage}_Cell Subtype.png"), dpi=300
        # )
        # plt.clf()

        pd.DataFrame(index=lin_adata.obs_names, data=lin_adata.obsm["X_umap"]).to_csv(
            f"{figures_lin}/{Lineage}_umap.csv"
        )

        with pd.ExcelWriter(
            f"{figures_lin}/{Lineage}_metadata_counts.xlsx", engine="xlsxwriter"
        ) as writer:
            obs_list = ["Cell Subtype", "Patient"]
            num_obs = len(obs_list) + 1
            for ind in range(0, num_obs):
                for subset in itertools.combinations(obs_list, ind):
                    if len(subset) != 0:
                        subset = list(subset)
                        if len(subset) == 1:
                            key = subset[0]
                            lin_adata.obs[key].value_counts().to_excel(
                                writer, sheet_name=key
                            )
                        else:
                            key = "_".join(subset)
                            lin_adata.obs.groupby(subset[:-1])[
                                subset[-1]
                            ].value_counts().to_excel(writer, sheet_name=key[:31])
        sc.tl.rank_genes_groups(lin_adata, "Cell Subtype", method="wilcoxon")
        sc.pl.rank_genes_groups_dotplot(
            lin_adata,
            groupby="Cell Subtype",
            n_genes=int(100 / len(lin_adata.obs[f"Cell Subtype"].unique())),
            show=False,
            save=f"_{Lineage}_Cell Subtype_markers.png",
        )


        def obs_key_wise_subsampling(adata, obs_key, N):
            '''
            Subsample each class to same cell numbers (N). Classes are given by obs_key pointing to categorical in adata.obs.
            '''
            counts = adata.obs[obs_key].value_counts()
            # subsample indices per group defined by obs_key
            indices = [np.random.choice(adata.obs_names[adata.obs[obs_key] == group], size=N, replace=True) for group in
                       counts.index]
            selection = np.hstack(np.array(indices))
            return adata[selection].copy()


        subsampled = obs_key_wise_subsampling(lin_adata, 'Cell Subtype', 500)
        sc.tl.rank_genes_groups(
            subsampled,
            "Cell Subtype",
            method="wilcoxon",
            pts=True,
            key_added="rank_genes_groups_Cell Subtype",
        )
        with pd.ExcelWriter(
                f"{figures_lin}/{Lineage}_Cell Subtype_markers.xlsx", engine="xlsxwriter"
        ) as writer:
            for ct in subsampled.obs["Cell Subtype"].cat.categories:
                df = sc.get.rank_genes_groups_df(
                    subsampled, key="rank_genes_groups_Cell Subtype", group=ct
                )
                df.index = df['names']
                df.to_excel(writer, sheet_name=f"{ct} v rest"[:31])

        ##  ct v every other ct
        figures_lin_comp = f"{figures_lin}/cell_type_comparisons"
        os.makedirs(figures_lin_comp, exist_ok=True)
        for ct in lin_adata.obs["Cell Subtype"].cat.categories:
            print(ct)
            with pd.ExcelWriter(
                    f"{figures_lin_comp}/{ct}.xlsx", engine="xlsxwriter"
            ) as writer:
                for ct2 in lin_adata.obs["Cell Subtype"].cat.categories:
                    cts_adata = lin_adata[lin_adata.obs["Cell Subtype"].isin([ct, ct2])]
                    try:
                        sc.tl.rank_genes_groups(
                            cts_adata,
                            "Cell Subtype",
                            groups=[ct, ct2],
                            method="wilcoxon",
                            pts=True,
                            key_added="rank_genes_groups_Cell Subtype",
                        )
                        df = sc.get.rank_genes_groups_df(
                            cts_adata, key="rank_genes_groups_Cell Subtype", group=ct
                        )
                        df.index = df['names']
                        df.to_excel(writer, sheet_name=f"{ct} v {ct2}"[:31])
                    except:
                        print(ct)
                        print(ct2)
                        print('no comp')
    print('Done witn Lineages')
    adata = adata[~adata.obs['Cell Subtype'].isna()]
    adata = adata[adata.obs['Cell Subtype'] != 'low-quality']
    sc.tl.umap(adata, min_dist=0.5)
    print(sorted(adata.obs['Cell Subtype'].unique()))
    sc.settings.figdir = f"{figures}/tissue_embedding"
    color_ls = [
                 'leiden',
                 'Lineage',
                 'log10_n_genes_by_umis',
                 "Patient",
                 'GroupContract',
                 'log10_total_umis',
        'Cell Subtype'
                 ]
    for gene in color_ls:
        if gene == 'leiden':
            sc.pl.umap(adata, legend_loc='on data', color=gene, alpha=0.5, show = False,save=f'_{gene}')
        else:
            sc.pl.umap(adata, color=gene, alpha=0.3, show = False, save=f'_{gene}')
    with pd.ExcelWriter(f'{figures}/tissue_embedding/metadata_counts.xlsx', engine='xlsxwriter') as writer:
        obs_list = ['Lineage', 'GroupContract','Cell Subtype',"Patient"]
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

    adata.write(f'{data}/uterus_processed_celltyped.gz.h5ad', compression='gzip')


