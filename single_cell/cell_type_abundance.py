import muon as mu
import os
import scanpy as sc
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker


figures = '/home/carsten/alvira_bioinformatics/uterus/data/figures/cell_type_abundance'
data = '/home/carsten/alvira_bioinformatics/uterus/data/single_cell_files/scanpy_files'
os.makedirs(figures, exist_ok=True)
sc.set_figure_params(dpi=300, dpi_save=300, format="png")
sc.settings.figdir = figures

if __name__ == "__main__":
    adata = sc.read(f'{data}/uterus_processed_celltyped.gz.h5ad')
    df = adata.obs
    p_gc_dict = pd.Series(df.GroupContract.values, index=df.Patient).to_dict()
    df2 = df.groupby('Patient')['Cell Subtype'].value_counts(normalize=True).mul(100).rename(
        '% cells sampled').reset_index()
    df2['Cell Subtype'] = df2['level_1']
    df2['GroupContract'] = [p_gc_dict[x] for x in df2.Patient]
    df2.to_csv(f'{figures}/cellsubtype_abunances.csv')

    print(df2)
    hue_order = ['TNL-BC','TNL-GC','PNL-ND','TL-ND']
    for Lineage in adata.obs['Lineage'].cat.categories:
        lin_adata = adata[adata.obs['Lineage'] == Lineage].copy()
        df3 = df2.loc[df2['Cell Subtype'].isin(lin_adata.obs['Cell Subtype'].cat.categories.tolist())].copy()
        df3['Cell Subtype'] = df3['Cell Subtype'].astype('string')
        order = lin_adata.obs['Cell Subtype'].cat.categories
        sns.catplot(data=df3,
                    x='Cell Subtype',
                    y='% cells sampled',
                    hue='GroupContract',
                    order=order,
                    hue_order = hue_order,
                    errorbar='se',
                    kind='bar')
        plt.xticks(rotation=90)
        plt.yscale('log')
        if Lineage=='Immune':
            plt.ylim(0.01, 20)
        ax = plt.gca()
        ax.yaxis.set_major_formatter(ticker.ScalarFormatter())  # set regular formatting
        plt.title(f'{Lineage}')
        plt.savefig(f'{figures}/{Lineage}_abunance_by_GroupContract_all.png', bbox_inches='tight')

        df_lin = lin_adata.obs
        df_lin_plot = df_lin.groupby('Patient')['Cell Subtype'].value_counts(normalize=True).mul(100).rename(
        f'% {Lineage} cells sampled').reset_index()
        df_lin_plot['GroupContract'] = [p_gc_dict[x] for x in df_lin_plot.Patient]
        df_lin_plot['Cell Subtype'] = df_lin_plot['level_1']
        sns.catplot(data=df_lin_plot,
                    x='Cell Subtype',
                    y=f'% {Lineage} cells sampled',
                    hue='GroupContract',
                    order=order,
                    hue_order=hue_order,
                    errorbar='se',
                    kind='bar')
        plt.xticks(rotation=90)
        plt.yscale('log')
        ax = plt.gca()
        ax.yaxis.set_major_formatter(ticker.ScalarFormatter())  # set regular formatting
        plt.title(f'{Lineage}')
        plt.savefig(f'{figures}/{Lineage}_abunance_by_GroupContract_{Lineage}.png', bbox_inches='tight')
    df = adata.obs
    df2 = df.groupby('Patient')['Lineage'].value_counts(normalize=True).mul(100).rename(
        '% cells sampled').reset_index()
    df2['Lineage'] = df2['level_1']
    df2['GroupContract'] = [p_gc_dict[x] for x in df2.Patient]
    df2.to_csv(f'{figures}/Lineage_abunances.csv')

    sns.catplot(data=df2,
                x='Lineage',
                y='% cells sampled',
                hue='GroupContract',
                hue_order=hue_order,
                errorbar='se',
                kind='bar')
    plt.xticks(rotation=90)
    plt.savefig(f'{figures}/Lineage_abunance_by_GroupContract.png', bbox_inches='tight')

