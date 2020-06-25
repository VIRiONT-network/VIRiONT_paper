#!usr/bin/en python3
import os

configfile : "config.yaml"

datapath=config['PathToData']
resultpath=config['PathToResult']

#get all barcodes in a list after demultiplexing
(BARCODE) = os.listdir(datapath)

#make a dictionary storing key=barcode, value=reads
list_fastq={}
for barcode in os.listdir(datapath):
    list_fastq[barcode]=os.listdir(datapath+barcode)

#final output
rule all:
    input:
        merged_file = expand(resultpath+'MERGED/{barcode}_merged.fastq',barcode=BARCODE),
       
#concatenate all fastq files 
rule merge:
    input: 
        lambda wildcards: expand(datapath+"{barcode}", barcode=BARCODE)
    output: 
        merged_fastq = resultpath+"MERGED/{barcode}_merged.fastq"  
    params:
        path = datapath
    shell: 
        """
        cat {params.path}{wildcards.barcode}/* > {output}
        """
