#!usr/bin/en python3
import os

pathdata="/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/"

BARCODE_list=os.listdir(pathdata)
barcode=BARCODE_list

rule merge:
    input:
        fastq_files = pathdata+'DATASET/{barcode}/*.fastq'
    wildcard_constraints:
        barcode=BARCODE_list       
    output:
        merged_fq = pathdata+'MERGED_FQ/{barcode}_merged.fastq'    
    shell:
        "cat {input} > {output} "           