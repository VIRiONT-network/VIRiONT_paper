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
rep_report="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/" #Define destination path for the text nohup report
data_loc="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/DATA_HBV/" #Define path where the "barcode*" are stored
result_loc="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/" #Define path where storing analysis results 
ref_loc="ref/HBV_REF.fasta" # Path to fastafile containing genotype sequences for blast
ref_table="analysis/table_analysis.csv" # Table location containing reference list for the blastn analysis. /!\ Column must correspond to the fasta headers in ref_loc.
analysis_name="HBV_REF" # Analysis Name. /!\ name must correpond to the choosen column name in ref_table.
thread_number=8 #Define number of threads to use for the analysis

snakemake -s viralION.py \
    --core $thread_number \
    --config PathToData=$data_loc \
             PathToResult=$result_loc \
             PathToReference=$ref_loc \
             AnalysisName=$analysis_name \
             AnalysisTable=$ref_table
         

 Rscript script/Blastn_analysis.R /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/BLASTN_RESULT/barcode05_fmt.txt             
 analysis/table_analysis.csv 
 HBV_REF             
 /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/BLASTN_ANALYSIS/barcode05_read-list.txt 
 /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/BLASTN_ANALYSIS/barcode05_bestref.txt 
 /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/BLASTN_ANALYSIS/barcode05_barplot.png

