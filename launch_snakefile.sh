#!/bin/bash

################################################################################
##########################    CONFIGURATION    #################################
################################################################################
data_loc="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/DATA_HDV/" #Define path where the "barcode*" are stored
result_loc="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/RESULT_HDV/" #Define path where storing analysis results 
ref_loc="ref/ICTVHDV.fasta" # Path to fastafile containing genotype sequences for blast
ref_table="analysis/table_analysis.csv" # Table location containing reference list for the blastn analysis. /!\ Column must correspond to the fasta headers in ref_loc.
analysis_name="ICTVHDV" # Analysis Name. /!\ name must correpond to the choosen column name in ref_table.
thread_number=8 #Define number of threads to use for the analysis
################################################################################


################################################################################
#########################    LAUNCH SNAKEMAKE    ###############################
################################################################################
snakemake -s viralION.py \
    --use-conda \
    --core $thread_number \
    --config PathToData=$data_loc \
             PathToResult=$result_loc \
             PathToReference=$ref_loc \
             AnalysisName=$analysis_name \
             AnalysisTable=$ref_table
################################################################################

#Generate workflow dag
#snakemake -s viralION.py --rulegraph \
#--config PathToData=$data_loc \
#             PathToResult=$result_loc \
#             PathToReference=$ref_loc \
#             AnalysisName=$analysis_name \
#             AnalysisTable=$ref_table| dot -Tpdf > documents/workflow.pdf





