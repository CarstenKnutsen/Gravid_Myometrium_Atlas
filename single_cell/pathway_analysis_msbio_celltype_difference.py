''' Run Metascape on all  DEGs from uterus snRNA
author: Carsten Knutsen
Date: Octonrt 31st 2023
conda environment: "msbio"
MSBio v3.5.230101
'''

import subprocess as sp
import os
import shutil
import pandas as pd
import numpy as np
import time
import json

fns = ['/home/carsten/alvira_bioinformatics/uterus/data/figures/tissue_embedding/Mesenchymal/cell_type_comparisons/Matrix fibroblast.xlsx',
       '/home/carsten/alvira_bioinformatics/uterus/data/figures/tissue_embedding/Mesenchymal/cell_type_comparisons/Uterine smooth muscle.xlsx',
       '/home/carsten/alvira_bioinformatics/uterus/data/figures/tissue_embedding/Endothelial/cell_type_comparisons/Capillary.xlsx'
       ]
output_fol = '/home/carsten/alvira_bioinformatics/uterus/data/figures/deg/cell_type_comparisons_pathway'

def deg_to_metascape(deg_fn,
                         output_fol,
                         msbio_fol = '/home/carsten/MSBio/msbio_v3.5.20230101/',
                         species='9606', # Takes taxonomy id 9606 is human
                     top_n_genes=100
):
    '''Takes a .xlsx file with DEGs and runs metascape on up and downregulated genes from each sheet

    !!!! Makes a new folder. Will delete existing folder and all contents if run on with output_fol as existing folder
    '''

    ## tmp file paths for msbio
    print('writing tmp files')
    input_fol = f'{msbio_fol}data/tmp/input/'
    tmp_output_fol = f'{msbio_fol}data/tmp/output/'
    job_file = f'{input_fol}batch.job'
    unlock = f'sudo chown -R $USER ~/MSBio'
    sp.run(unlock, cwd = msbio_fol, shell=True)
    ### write csv files for gene lists
    for path in [input_fol, tmp_output_fol]:
        if os.path.exists(path):
            shutil.rmtree(path)
        os.makedirs(path, exist_ok=True)
    degs = pd.read_excel(deg_fn, sheet_name=None, index_col=0)
    for key in degs.keys():
        df = degs[key]
        df.head(top_n_genes).index.to_series(name='Gene').to_csv(f'{input_fol}{key}_up.csv', index=False) ## write the upregulated genes tmp files
        df.tail(top_n_genes).index.to_series(name='Gene').to_csv(f'{input_fol}{key}_down.csv', index=False) ## write the downregulated genes tmp files
    ### make job file
    json_dict = {}
    for fn in os.listdir(input_fol):
        tmp_output_fol2 = f'{tmp_output_fol}{fn.split(".")[0]}'
        os.makedirs(tmp_output_fol2, exist_ok=True)
        tmp_dict = {}
        tmp_dict['input'] = f'/{"/".join(f"{input_fol}{fn}".split("/")[-4:])}'
        tmp_dict['output'] =f'/{"/".join(f"{tmp_output_fol2}".split("/")[-4:])}'
        json_dict[fn.split('.')[0]] = tmp_dict
    with open(job_file, "w") as outfile:
        outfile.write(
            '\n'.join(json.dumps(json_dict[k]) for k in sorted(json_dict.keys()))
        )
    ## run MSBIO
    print('running MSBio')
    # commands
    up = f'bin/up.sh'
    job_run = f'bin/ms.sh /data/tmp/input/batch.job -u -S {species} -T {species}'
    down = f'bin/down.sh'
    unlock = f'sudo chown -R $USER ~/MSBio'
    # run commands
    sp.run(up, cwd = msbio_fol, shell=True)
    sp.run(job_run, cwd = msbio_fol, shell=True)
    sp.run(down, cwd = msbio_fol, shell=True)
    sp.run(unlock, cwd = msbio_fol, shell=True)
    # copy to final location
    if os.path.exists(output_fol):
        shutil.rmtree(output_fol)
    shutil.copytree(tmp_output_fol, output_fol)
    for path in [input_fol, tmp_output_fol]:
        if os.path.exists(path):
            shutil.rmtree(path)

if __name__ =='__main__':
    for fn in fns:
        ct = fn.split('/')[-1].split('.')[0]
        output_fn_fol = f'{output_fol}/{ct}'
        os.makedirs(output_fn_fol, exist_ok=True)
        deg_to_metascape(fn, output_fn_fol)



