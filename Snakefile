#!usr/bin/en python3
#singularity shell /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/nanopore.simg
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
        merged_file = expand('MERGED/{barcode}_merged.fastq',barcode=BARCODE),
        trimmed_file = expand('TRIMMED/{barcode}_trimmed.fastq',barcode=BARCODE),
        converted_fastq = expand("FASTA/{barcode}.fasta", barcode=BARCODE)
    params:
        trimmomatic =  "~/git/MinION_HBV/tool/Trimmomatic-0.39/trimmomatic-0.39.jar" ,
        trim_param = 500 

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

#trimming
rule trimming:
    input:
        merged_fastq = rules.merge.output.merged_fastq
    output:
        trimmed_fastq = "TRIMMED/{barcode}_trimmed.fastq"
   
    shell:
        """
        java -jar {rules.all.params.trimmomatic} SE -phred33 {input} {output} MINLEN:{rules.all.params.trim_param}
        """

rule converting:
    input:
        trimmed_fastq = rules.trimming.output.trimmed_fastq
    output:
        converted_fastq = "FASTA/{barcode}.fasta" 
    shell:
        """
        seqtk seq -A {input} > {output}
        """               
        