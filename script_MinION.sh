#!/bin/bash

################################################################################
##########################    CONFIGURATION    #################################
################################################################################
singularity_img="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/nanopore.simg"
#singularity shell /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/nanopore.simg
analyse="HBV"
heure=$(date +%H%M)
jour=$(date +%Y%m%d) 
rep_report="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/"
data_loc="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/"
result_loc="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/"
trim_value=500
################################################################################

################################################################################
#########################    LAUNCH SNAKEMAKE    ###############################
################################################################################
#nohup singularity exec $singularity_img snakemake --dryrun > ${rep_report}report_${jour}_${heure}.txt 
singularity exec $singularity_img snakemake \
    --config PathToData=$data_loc \
             PathToResult=$result_loc \
             TrimmParam=$trim_value  
    # --dryrun
################################################################################
###option:
# --dryrun => fait tourner le pipeline à vide pour controler
# --config => fait passer des parametres

#V2
source activate viralION
rep_report="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/"
data_loc="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/"
result_loc="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/"

snakemake -s viralION.py \
    --core all \
    --config PathToData=$data_loc \
             PathToResult=$result_loc
~/.bash_profile
NanoFilt --quality 10 --length 100 --maxlength 1500 \
    /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/MERGED/barcode05_merged.fastq > \
    /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/TRIMMED/barcode05_trimmed.fastq             