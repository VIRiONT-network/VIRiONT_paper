#!/bin/bash

################################################################################
##########################    CONFIGURATION    #################################
################################################################################

################################################################################
#######################    GENERAL PARAMETERS    ###############################
################################################################################
#fastq location / Define path where the "barcode*" rep are stored
data_loc="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/DATA_HBV_PARIS" 
#output location / Define path where storing analysis results 
result_loc="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/TEST_BM/" 
#custom reference file to use /  Path to fastafile containing reference sequences for blast
ref_loc="ref/HBV_REF.fasta" 
#core number / Define number of threads to use for the analysis
thread_number=8
#memory cost in mb / Define number of threads to use for the analysis
mem_cost=32000
################################################################################

################################################################################
#################    TRIMMING/FILTERING PARAMETERS    ##########################
################################################################################
#min length for read filtering
min_length=500
#max length for read filtering
max_length=3500
#Remove N nucleotide for primer 5'
head=0
#Remove N nucleotide for primer 3'
tail=0
################################################################################

################################################################################
######################    CONSENSUS PARAMETERS    ##############################
################################################################################
#minor variant frequency
Vfreq=0.5 
#minimun coverage necessary to generate a consensus / N is called if below threshold
mincov=20
################################################################################




################################################################################
#########################    LAUNCH SNAKEMAKE    ###############################
################################################################################
snakemake -s VIRiONT_BM.py \
    --use-conda \
    --core $thread_number \
    --resources mem_mb=$mem_cost \
    --config PathToData=$data_loc \
             PathToResult=$result_loc \
             PathToReference=$ref_loc \
             AnalysisTable=$ref_table \
             Lmin=$min_length \
             Lmax=$max_length \
             headcrop=$head \
             tailcrop=$tail \
             variantfrequency=$Vfreq \
             mincov=$mincov
################################################################################

#snakemake -s VIRiONT.py --rulegraph \
#    --use-conda \
#    --core $thread_number \
#    --resources mem_mb=$mem_cost \
#    --config PathToData=$data_loc \
#             PathToResult=$result_loc \
#             PathToReference=$ref_loc \
#             AnalysisTable=$ref_table \
#             Lmin=$min_length \
#             Lmax=$max_length \
#             headcrop=$head \
#             tailcrop=$tail \
#             variantfrequency=$Vfreq | dot -Tpng > documents/workflow.png

