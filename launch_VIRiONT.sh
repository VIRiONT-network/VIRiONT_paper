#!/bin/bash
#k5start -U -f /home/chu-lyon.fr/regueex/login.kt -- nohup ./launch_VIRiONT.sh &
################################################################################
##########################    CONFIGURATION    #################################
################################################################################

################################################################################
#######################    GENERAL PARAMETERS    ###############################
################################################################################
#fastq location / Define path where the "barcode*" rep are stored
data_loc="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/DATA_HDV_TEST/" 
#output location / Define path where storing analysis results 
result_loc="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CARO_PIPELINE/HDV_SGT_NEWREFWITHPRIMER_PILEUPQ20/" 
#custom reference file to use /  Path to fastafile containing reference sequences for blast
ref_loc="ref/HDV_subtype_V2.fasta" 
#core number / Define number of threads to use for the analysis
thread_number=8
#memory cost in mb / Define number of threads to use for the analysis
mem_cost=62000
################################################################################

################################################################################
#################    TRIMMING/FILTERING PARAMETERS    ##########################
################################################################################
#min length for read filtering
min_length=500
#max length for read filtering
max_length=2000
#average read quality for filtering
quality=0
#Remove N 5' nucleotides from each filtered read 
head=23
#Remove N 3' nucleotides from each filtered read 
tail=23
################################################################################

################################################################################
#################    VARIANT CALLING PARAMETERS    ##########################
################################################################################
#maximum depth for samtools mpileup
maxdepth=1000000
#base quality cutoff for samtools mpileup
basequal=20
################################################################################

################################################################################
######################    CONSENSUS PARAMETERS    ##############################
################################################################################
#minor variant frequency
Vfreq=0.50
#minimun coverage necessary to generate a consensus / N is called if below threshold
mincov=20
################################################################################

################################################################################
###################    MULTI-INFECTION PARAMETER    ############################
################################################################################
#Multi infection threshold cutoff in percent / cutoff=count(Blastref_reads)/count(majoritaryBlastref_reads)*100
MI_cutoff=30
################################################################################


################################################################################
#########################    LAUNCH SNAKEMAKE    ###############################
################################################################################

snakemake -s smk_VIRiONT_1.py -p  \
    --use-conda \
    --core $thread_number \
    --resources mem_mb=$mem_cost \
    --config PathToData=$data_loc \
             PathToResult=$result_loc \
             PathToReference=$ref_loc \
             Lmin=$min_length \
             Lmax=$max_length \
             quality=$quality \
             headcrop=$head \
             tailcrop=$tail \
             multiinf=$MI_cutoff

snakemake -s smk_VIRiONT_2.py -p \
    --use-conda \
    --core $thread_number \
    --resources mem_mb=$mem_cost \
    --config PathToData=$data_loc \
             PathToResult=$result_loc \
             PathToReference=$ref_loc \
             variantfrequency=$Vfreq \
             mincov=$mincov \
             Lmin=$min_length \
             Lmax=$max_length \
             headcrop=$head \
             tailcrop=$tail \
             depth=$maxdepth \
             basequality=$basequal \
             multiinf=$MI_cutoff

chmod -R 777 $result_loc
