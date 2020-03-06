#!/bin/bash

################################################################################
##########################    CONFIGURATION    #################################
################################################################################
singularity_img="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/IMG_SINGULARITY/guppy_CPU.simg"
input_data="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/FAST5"
output_data="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/TEST_GUPPY_7T/"
################################################################################

#cp /srv/nfs/ngs-stockage/NGS_Virologie/CCharre/raw_data_run_DamienCaro/fast5_fail/* /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/FAST5/ > /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CP1.txt &
#cp /srv/nfs/ngs-stockage/NGS_Virologie/CCharre/raw_data_run_DamienCaro/fast5_pass/* /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/FAST5/ > /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/CP2.txt &

################################################################################
#########################      LAUNCH GUPPY      ###############################
################################################################################


k5start -U -f /home/chu-lyon.fr/regueex/login.kt  -- \
nohup singularity exec $singularity_img guppy_basecaller \
    --flowcell FLO-MIN106 \
    --kit SQK-PBK004 \
    --input_path $input_data \
    --save_path $output_data \
    --cpu_threads_per_caller 7 \ #/!\ A REGLER AVANT DE LANCER
    --num_callers 1 \
    --records_per_fastq 0 \
    --disable_pings \
    --recursive \
    --qscore_filtering \
    --min_qscore 7 > $output_data/rapport_nohup.txt &
