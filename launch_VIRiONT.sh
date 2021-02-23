#!/bin/bash

################################################################################
#########################    LAUNCH SNAKEMAKE    ###############################
################################################################################

max_thread=`grep "thread_number"  config/params.yaml | cut -d":" --fields=2 | tr -d ' '`
#echo $max_thread
mem_cost=`grep "mem_cost"  config/params.yaml | cut -d":" --fields=2 | tr -d ' '`
#echo $mem_cost

#grep "data_loc"  config/params.yaml | cut -d":" --fields=2 | tr -d ' '


################################################################################
#########################    LAUNCH SNAKEMAKE    ###############################
################################################################################

snakemake -s smk_VIRiONT_1.py -p  \
    --use-conda \
    --cores $max_thread \
    --resources mem_mb=$mem_cost

snakemake -s smk_VIRiONT_2.py -p  \
    --use-conda \
    --cores $max_thread \
    --resources mem_mb=$mem_cost

