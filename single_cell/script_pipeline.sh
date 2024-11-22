#!/bin/bash
echo 'starting'
#conda activate soupxR
#Rscript soupX.R
#conda activate uterus_sc
#python create_count_table.py
#python qc.py
#python all_lineage_cluster_embed.py
#python cell_typing.py
#python cell_type_abundance.py
#conda activate pseudobulk
#python differential_expression.py
conda activate msbio
python pathway_analysis_msbio.py
#python pathway_analysis_msbio_celltype_difference.py

