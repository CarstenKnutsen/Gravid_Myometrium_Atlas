"""
Goal: Run differential expression for comparison groups across labor, both pseudobulk and wilcoxon
Author:Carsten Knutsen
Date:231115
conda_env:pseudobulk
Notes: Adapted from decoupler tutorial https://decoupler-py.readthedocs.io/en/latest/notebooks/pseudobulk.html
"""


import scanpy as sc
import decoupler as dc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg') # https://stackoverflow.com/questions/27147300/matplotlib-tcl-asyncdelete-async-handler-deleted-by-the-wrong-thread
import matplotlib.pyplot as plt
import os
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


output = "/home/carsten/alvira_bioinformatics/uterus/data/figures/deg"
os.makedirs(output, exist_ok=True)
data = "/home/carsten/alvira_bioinformatics/uterus/data/single_cell_files/scanpy_files"
os.makedirs(output, exist_ok=True)
sc.set_figure_params(dpi=300, dpi_save=300, format="png")
sc.settings.figdir = output
comparison_dictionary = {
    "TNL_GC_v_TNL_BC": ["TNL-BC", "TNL-GC"],
    "TNL_GC_v_PNL_ND": ["PNL-ND", "TNL-GC"],
    "TNL_GC_v_TL_ND": ["TL-ND", "TNL-GC"],
    "TNL_v_TL ": ["TL", "TNL"],
}
if __name__ == "__main__":
    adata = sc.read(f"{data}/uterus_processed_celltyped.gz.h5ad")
    # # ### PSeudobulk
    output_p = f"{output}/pseudobulk"
    os.makedirs(output_p, exist_ok=True)
    pdata = dc.get_pseudobulk(
        adata,
        sample_col="Patient",
        groups_col="Cell Subtype",
        layer="soupx",
        mode="sum",
        min_cells=0,
        min_counts=0,
    )
    dc.plot_psbulk_samples(
        pdata,
        groupby=["Patient", "Cell Subtype"],
        figsize=(11, 3),
        dpi=300,
        save=f"{output_p}/pseudo_bulk_samples.png",
    )
    plt.close()
    pp_adata = pdata.copy()
    sc.pp.normalize_total(pp_adata, target_sum=1e6)
    sc.pp.log1p(pp_adata)
    sc.pp.scale(pp_adata, max_value=10)
    sc.tl.pca(pp_adata, n_comps=10)
    sc.pp.neighbors(pp_adata)
    sc.tl.umap(pp_adata, min_dist=0.1)
    sc.pl.umap(
        pp_adata,
        color=["GroupContract", "Cell Subtype"],
        ncols=2,
        size=300,
        show=False,
        save="_pseudobulk_meta.png",
    )

    comp_dict = {}
    for key in comparison_dictionary.keys():
        if key ==  "TNL_v_TL ":
            comparison = 'Group'
        else:
            comparison = 'GroupContract'
        comp_list = comparison_dictionary[key]
        print(comp_list)
        compare_pseudo = pdata[pdata.obs[comparison].isin(comp_list)]
        print(compare_pseudo)
        ct_dict = {}
        for ct in adata.obs["Cell Subtype"].unique():
            ct_adata = compare_pseudo[compare_pseudo.obs["Cell Subtype"] == ct]
            print(ct_adata)
            genes = dc.filter_by_expr(
                ct_adata, group=comparison, min_count=10, min_total_count=15
            )
            ct_adata = ct_adata[:, genes].copy()
            if len(genes) < 100:
                print(ct)
                print("NOT ENOUGH GENES")
                continue
            print(ct_adata)
            if len(ct_adata.obs_names) < 4:
                print(ct)
                print('Not enough samples')
                continue
            if len(ct_adata.obs[comparison].unique())==1:
                print(ct)
                print(ct_adata.obs[comparison].unique())
                print('Only found in one group')
                continue
            dds = DeseqDataSet(
                adata=ct_adata,
                design_factors=comparison,
                ref_level=[
                    comparison,
                    sorted(ct_adata.obs[comparison].unique())[0],
                ],
                refit_cooks=True,
                n_cpus=8,
            )
            dds.deseq2()
            contrast = [comparison] + sorted(
                ct_adata.obs[comparison].unique()
            )
            coeff = f"{comparison}_{contrast[-1]}_vs_{contrast[-2]}"
            stat_res = DeseqStats(dds, contrast=contrast, n_cpus=8)
            stat_res.summary()
            stat_res.lfc_shrink(coeff=coeff)
            results_df = stat_res.results_df
            ct_dict[ct] = results_df.sort_values("pvalue")
        comp_dict[key] = ct_dict
    for key in comp_dict.keys():
        ct_dict = comp_dict[key]
        with pd.ExcelWriter(
            f"{output_p}/{key}_pseudobulk_comparisons.xlsx", engine="xlsxwriter"
        ) as writer:
            for key2 in sorted(ct_dict.keys()):
                ct_df = ct_dict[key2]
                ct_df.to_excel(writer, sheet_name=f"{key2}"[:31])
    ### Wilcoxon rank sum
    output_w = f"{output}/wilcoxon"
    os.makedirs(output_w, exist_ok=True)

    for key in comparison_dictionary.keys():
        if key == "TNL_v_TL ":
            comparison = 'Group'
        else:
            comparison = 'GroupContract'
        comp_list = comparison_dictionary[key]
        compare_adata = adata[adata.obs[comparison].isin(comp_list)]
        ct_dict = {}
        with pd.ExcelWriter(
            f"{output_w}/{key}_wilcoxon_comparisons.xlsx", engine="xlsxwriter"
        ) as writer:
            for ct in sorted(adata.obs["Cell Subtype"].unique()):
                ct_adata = compare_adata[compare_adata.obs["Cell Subtype"] == ct]
                if ct_adata.obs[comparison].value_counts().min() < 5:
                    print(ct)
                    print(ct_adata.obs[comparison].value_counts())
                    print('Too few cells in group')
                    continue
                print(ct_adata)
                try:
                    sc.tl.rank_genes_groups(
                        ct_adata,
                        groupby=comparison,
                        method="wilcoxon",
                        pts=True,
                        key_added="group_contract",
                    )
                    df = sc.get.rank_genes_groups_df(
                        ct_adata, key="group_contract", group=[comp_list[0]]
                    )
                    df.to_excel(writer, sheet_name=f"{ct}"[:31])
                except:
                    print(ct)
                    print('No comparison made')
