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
        converted_fastq = expand("FASTA/{barcode}.fasta", barcode=BARCODE),
        R_data = expand("BLASTN/{barcode}_fmt.txt" ,barcode=BARCODE)
    params:
        trimmomatic =  "TOOL/Trimmomatic-0.39/trimmomatic-0.39.jar" ,
        trim_param = 500 ,
        ref_HVB = "DATA/HBV_REF.fasta",
        DB_HBV = "HBV_REF"

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
#convert fastq into fasta
rule converting:
    input:
        trimmed_fastq = rules.trimming.output.trimmed_fastq
    output:
        converted_fastq = "FASTA/{barcode}.fasta" 
    shell:
        """
        seqtk seq -A {input} > {output}
        """               

#build the blast database
rule make_db_HBV:
    input:
        ref_HVB = "DATA/HBV_REF.fasta"
    output:
        database = expand("DB/HBV_REF.{ext}", ext=["nhr", "nin", "nsq"])
    shell:
        """
        makeblastdb -in {input} -out DB/HBV_REF -input_type fasta -dbtype nucl
        """

#execute blastn
rule blastn:
    input: 
        fasta_file = rules.converting.output.converted_fastq,
        database = rules.make_db_HBV.output.database
    output:
        R_data = "BLASTN/{barcode}_fmt.txt"
    shell:
        """
        blastn -db DB/{rules.all.params.DB_HBV} -query {input.fasta_file} -outfmt 6 -out {output}
        """                        