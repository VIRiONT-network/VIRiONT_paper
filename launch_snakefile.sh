#!/bin/bash

################################################################################
##########################    CONFIGURATION    #################################
################################################################################
#fastq location
data_loc="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/DATA_HBV/" #Define path where the "barcode*" are stored
#output location
result_loc="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/RESULT_HBV/" #Define path where storing analysis results 
#custom reference file to use
ref_loc="ref/HBV_REF.fasta" # Path to fastafile containing genotype sequences for blast
#metadata location for the custom reference
ref_table="analysis/table_analysis.csv" # Table location containing reference list for the blastn analysis. /!\ Column must correspond to the fasta headers in ref_loc.
#core number
thread_number=8 #Define number of threads to use for the analysis
################################################################################


################################################################################
#########################    LAUNCH SNAKEMAKE    ###############################
################################################################################
snakemake -s VIRiONT.py \
    --use-conda \
    --use-singularity \
    --core $thread_number \
    --config PathToData=$data_loc \
             PathToResult=$result_loc \
             PathToReference=$ref_loc \
             AnalysisTable=$ref_table
################################################################################

#If troubles with lock:
#snakemake -s viralION.py \
#    --unlock \
#    --use-conda \
#    --core $thread_number \
#   --config PathToData=$data_loc \
#             PathToResult=$result_loc \
#             PathToReference=$ref_loc \
#             AnalysisTable=$ref_table

#Generate workflow dag
#snakemake -s viralION.py --rulegraph \
#--config PathToData=$data_loc \
#             PathToResult=$result_loc \
#             PathToReference=$ref_loc \
#             AnalysisTable=$ref_table| dot -Tpdf > documents/workflow.pdf
