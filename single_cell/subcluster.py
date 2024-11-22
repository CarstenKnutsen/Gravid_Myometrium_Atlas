"""
Goal: Run differential expression for comparison groups across labor, both pseudobulk and wilcoxon
Author:Carsten Knutsen
Date:231115
conda_env:pseudobulk
Notes: Adapted from decoupler tutorial https://decoupler-py.readthedocs.io/en/latest/notebooks/pseudobulk.html
"""


import scanpy as sc
import scanpy.external as sce
import pandas as pd
import matplotlib
matplotlib.use('Agg') # https://stackoverflow.com/questions/27147300/matplotlib-tcl-asyncdelete-async-handler-deleted-by-the-wrong-thread
import matplotlib.pyplot as plt
import os



output = "/home/carsten/alvira_bioinformatics/uterus/data/figures"
os.makedirs(output, exist_ok=True)
data = "/home/carsten/alvira_bioinformatics/uterus/data/single_cell_files/scanpy_files"
os.makedirs(output, exist_ok=True)
sc.set_figure_params(dpi=300, dpi_save=300, format="png")
sc.settings.figdir = output
def subcluster_celltype(adata,celltype,output_fol):
    output_ct = f'{output_fol}/subcluster/{celltype}'
    os.makedirs(output_ct, exist_ok=True)
    sc.set_figure_params(dpi=300, format="png")
    sc.settings.figdir = output_ct
    lin_adata = adata[adata.obs['Cell Subtype']==celltype]
    try:
        sc.pp.highly_variable_genes(lin_adata,
                                    batch_key="Patient"
                                    )
    except:
        sc.pp.highly_variable_genes(lin_adata
                            )
    sc.pp.pca(lin_adata, use_highly_variable=True)
    try:
        sce.pp.harmony_integrate(lin_adata, key="Patient", max_iter_harmony=50)
        sc.pp.neighbors(lin_adata, use_rep='X_pca_harmony')

    except:
        sc.pp.neighbors(lin_adata)

    sc.tl.leiden(
        lin_adata,
        key_added=f"leiden_{celltype}",
    )
    sc.tl.umap(lin_adata,min_dist=0.1)
    sc.tl.rank_genes_groups(lin_adata, f"leiden_{celltype}",pts=True, method="wilcoxon")
    print(lin_adata.obs[f"leiden_{celltype}"].cat.categories)
    sc.pl.rank_genes_groups_dotplot(
        lin_adata,
        groupby=f"leiden_{celltype}",
        n_genes=int(100 / len(lin_adata.obs[f"leiden_{celltype}"].unique())),
        show=False,
        save=f"{celltype}_leiden_markers.png",
    )

    for color in [f'leiden_{celltype}','Patient','GroupContract','Cell Subtype']:
        sc.pl.umap(lin_adata, color = color, show=False,save=color)
    with pd.ExcelWriter(
        f"{output_ct}/{celltype}_leiden_markers.xlsx", engine="xlsxwriter") as writer:
        for ld in sorted(lin_adata.obs[f"leiden_{celltype}"].unique()):
            df = sc.get.rank_genes_groups_df(
                lin_adata, key="rank_genes_groups", group=ld
            )
            df.to_excel(writer, sheet_name=f"{ld} v rest"[:31])
        lin_adata.write(f'{output_ct}/{celltype}_adata.gz.h5ad', compression='gzip')
if __name__ == "__main__":
    adata = sc.read(f"{data}/uterus_processed_celltyped.gz.h5ad")
    adata.uns['log1p']['base'] = None
## Subcluster each cell type
for celltype in adata.obs['Cell Subtype'].unique():
    subcluster_celltype(adata,celltype,output)