#!/bin/bash

################################################################################
##########################    CONFIGURATION    #################################
################################################################################
singularity_img="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/nanopore.simg"
analyse="HBV"
heure=$(date +%H%M)
jour=$(date +%Y%m%d) 
rep_report="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/RAPPORT_PIPELINE/"
################################################################################


################################################################################
#########################    LAUNCH SNAKEMAKE    ###############################
################################################################################
#nohup singularity exec $singularity_img snakemake --dryrun > ${rep_report}report_${jour}_${heure}.txt 
singularity exec $singularity_img snakemake --dryrun
################################################################################
###option:
# --dryrun => fait tourner le pipeline à vide pour controler