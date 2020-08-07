#!/bin/bash

################################################################################
##########################    CONFIGURATION    #################################
################################################################################
#fastq location
data_loc="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/DATA_HBV_PARIS/" #Define path where the "barcode*" are stored
#data_loc="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/POLIO/DATA/"
#output location
result_loc="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/RESULT_HBV_PARIS/" #Define path where storing analysis results 
#result_loc="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/POLIO/RESULTS/"
#custom reference file to use
ref_loc="ref/HBV_REF.fasta" # Path to fastafile containing genotype sequences for blast
#ref_loc="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/POLIO/POLIO.fasta"
#metadata location for the custom reference
ref_table="analysis/table_analysis.csv" # Table location containing reference list for the blastn analysis. /!\ Column must correspond to the fasta headers in ref_loc.
#ref_table="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/POLIO/table_analysis.csv" 
#min length for read filtering
min_length=1500
#max length for read filtering
max_length=3500
#Remove N nucleotide for primer 5'
head=0
#Remove N nucleotide for primer 3'
tail=0
#Multi infection threshold cutoff in percent. cutoff= Blastref/majoritaryBlastref*100
MI_cutoff=10
#minor variant frequency
Vfreq=0.5 
#core number
thread_number=10 #Define number of threads to use for the analysis
#memory cost in mb
mem_cost=32000 #Define number of threads to use for the analysis
################################################################################

################################################################################
#########################    LAUNCH SNAKEMAKE    ###############################
################################################################################
snakemake -s VIRiONT_MI1.py \
    --use-conda \
    --core $thread_number \
    --resources mem_mb=$mem_cost \
    --config PathToData=$data_loc \
             PathToResult=${result_loc}"MI_VERSION/" \
             PathToReference=$ref_loc \
             AnalysisTable=$ref_table \
             Lmin=$min_length \
             Lmax=$max_length \
             headcrop=$head \
             tailcrop=$tail \
             multiinf=$MI_cutoff

snakemake -s VIRiONT_MI2.py \
    --use-conda \
    --core $thread_number \
    --resources mem_mb=$mem_cost \
    --config PathToData=$data_loc \
             PathToResult=${result_loc}"MI_DETECTION/" \
             PathToReference=$ref_loc \
             variantfrequency=$Vfreq 

chmod 777 -R $data_loc


