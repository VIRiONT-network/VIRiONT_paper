#!/bin/bash

################################################################################
##########################    CONFIGURATION    #################################
################################################################################
#fastq location
data_loc="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/DATA_HDV/" #Define path where the "barcode*" are stored
#output location
result_loc="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/RESULT_HDV/" #Define path where storing analysis results 
#custom reference file to use
ref_loc="ref/ICTVHDV.fasta" # Path to fastafile containing genotype sequences for blast
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
#test medaka
#medaka:
medaka_variant -f /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/save_result/RESULT_HDV/REFSEQ/P4_HDV5.fasta \
     -i /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/save_result/RESULT_HDV/BAM/barcode03_sorted.bam \
     -o /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/save_result/RESULT_HDV/VCF/barcode03





#test callVarBam
clair.py callVarBam  \
    --threshold 0.5 \
    --haploid_precision \
    --chkpnt_fn "ont/model" \
    --ref_fn "/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/RESULT_HBV/REFSEQ/P3_GTJ.fasta" \
    --bam_fn "/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/RESULT_HBV/BAM/barcode02_sorted.bam" \
    --ctgName "P3_GTJ" \
    --sampleName "barcode03" \
    --call_fn "/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/RESULT_HBV/barcode02_testcutoffFA.vcf"




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
