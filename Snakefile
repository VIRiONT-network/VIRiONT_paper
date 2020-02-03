#!usr/bin/en python3
#singularity shell /srv/nfs/ngs-stockage/NGS_Virologie/NEXTSTRAIN/nextstrainV3.simg
#cp ~/git/MinION_HBV/Snakefile  /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR

import os

#Data repository location
workdir : "/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/"

#get all barcodes in a list after demultiplexing
(BARCODE) = os.listdir('DATASET/')

#make a dictionary storing key=barcode, value=reads
list_fastq={}
for barcode in os.listdir('DATASET/'):
    list_fastq[barcode]=os.listdir('DATASET/'+barcode)

#final output
rule all:
    input:
        merged_file = expand('MERGED/{barcode}_merged.fastq',barcode=BARCODE)

#concatenate all fastq files 
rule merge:
    input: 
        lambda wildcards: expand("DATASET/{barcode}/{read}", read=list_fastq[wildcards.barcode],barcode=BARCODE)
    output: 
        merged_fastq = "MERGED/{barcode}_merged.fastq"
    shell: 
        """
        cat DATASET/{wildcards.barcode}/* > {output}
        """