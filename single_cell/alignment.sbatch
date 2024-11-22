#!/usr/bin/bash
#SBATCH --job-name=alignment
#SBATCH --output=alignment.%j.out
#SBATCH --error=alignment.%j.err
#SBATCH --time=2-00:00:00
#SBATCH -p normal
#SBATCH -n 64
#SBATCH --mem=128GB	

ml biology
ml cellranger #Version 7.1.0
for dir in /scratch/users/cknutsen/23_uterus_contractility/raw_data/241001_sequencing/*/
do
    dir2=${dir%*/}   
    dir2="${dir2##*/}"
    echo $dir2   
    cellranger count --id="${dir2}" \
                     --transcriptome=/home/groups/calvira/software/cellranger_base/refdata-gex-GRCh38-2020-A \
                     --fastqs=$dir \

done
