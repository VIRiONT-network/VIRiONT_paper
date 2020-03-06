#!/bin/bash

################################################################################
##########################    CONFIGURATION    #################################
################################################################################
singularity_img="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/nanopore.simg"
#singularity shell /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/nanopore.simg
analyse="HBV"
heure=$(date +%H%M)
jour=$(date +%Y%m%d) 
rep_report="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/MinIONNE_CARO_PARIS/"
data_loc="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/MinIONNE_CARO_PARIS/"
result_loc="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PARIS_PIPELINE_OUTPUT/"
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