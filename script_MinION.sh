#!/bin/bash
path=$1
#Merge each fastas from each barcode rep
for barcode in `ls DATASET`
do
    cat DATASET/${barcode}/*.fastq > MERGE/${barcode}_merged.fastq
done   