#!usr/bin/en python3
import os

configfile : "config.yaml"

datapath=config['PathToData']
resultpath=config['PathToResult']

#get all barcodes in a list after demultiplexing
(BARCODE) = os.listdir(datapath)
#(GENOTYPE)={}
#GENOTYPE['genotype'] =("GTA","GTB","GTC","GTD","GTE","GTF","GTG","GTH","GTI","GTJ")

#make a dictionary storing key=barcode, value=reads
list_fastq={}
for barcode in os.listdir(datapath):
    list_fastq[barcode]=os.listdir(datapath+barcode)

#final output
rule all:
    input:
        merged_file = expand(resultpath+'MERGED/{barcode}_merged.fastq',barcode=BARCODE),
        trimmed_file = expand(resultpath+'TRIMMED/{barcode}_trimmed.fastq',barcode=BARCODE),
        converted_fastq = expand(resultpath+"FASTA/{barcode}.fasta", barcode=BARCODE),
        R_data = expand(resultpath+"BLASTN/{barcode}_fmt.txt" ,barcode=BARCODE),
        #Rresults = expand(datapath+"R_RESULT/{barcode}_list_{genotype}.txt",barcode=BARCODE,genotype=GENOTYPE['genotype'])
        Rplot = expand(resultpath+"RDATA/{barcode}_barplot.png" ,barcode=BARCODE),
    params:
        trimmomatic =  "tool/Trimmomatic-0.39/trimmomatic-0.39.jar" ,
        trim_param = 500 ,
        ref_HVB = "DATA/HBV_REF.fasta",
        DB_HBV = "HBV_REF"



#concatenate all fastq files 
rule merge:
    input: 
        lambda wildcards: expand(datapath+"{barcode}/{read}", read=list_fastq[wildcards.barcode],barcode=BARCODE)
    output: 
        merged_fastq = resultpath+"MERGED/{barcode}_merged.fastq"  
    params:
        path = datapath
    shell: 
        """
        cat {params.path}{wildcards.barcode}/* > {output}
        """

#trimming
rule trimming:
    input:
        merged_fastq = rules.merge.output.merged_fastq
    output:
        trimmed_fastq = resultpath+"TRIMMED/{barcode}_trimmed.fastq"
   
    shell:
        """
        java -jar {rules.all.params.trimmomatic} SE -phred33 {input} {output} MINLEN:{rules.all.params.trim_param}
        """

#convert fastq into fasta
rule converting:
    input:
        trimmed_fastq = rules.trimming.output.trimmed_fastq
    output:
        converted_fastq = resultpath+"FASTA/{barcode}.fasta" 
    shell:
        """
        seqtk seq -A {input} > {output}
        """               

#build the blast database
rule make_db_HBV:
    input:
        ref_HVB = "ref/HBV_REF.fasta"
    output:
        database = expand(resultpath+"DB/HBV_REF.{ext}", ext=["nhr", "nin", "nsq"])
    params:
        path = resultpath         
    shell:
        """
        makeblastdb -in {input} -out {params.path}/DB/HBV_REF -input_type fasta -dbtype nucl
        """

#execute blastn
rule blastn:
    input: 
        fasta_file = rules.converting.output.converted_fastq,
        database = rules.make_db_HBV.output.database
    output:
        R_data = resultpath+"BLASTN/{barcode}_fmt.txt"
    params:
        path = resultpath        
    shell:
        """
        blastn -db {params.path}/DB/{rules.all.params.DB_HBV} -query {input.fasta_file} -outfmt 6 -out {output}
        """                        

rule R_HBV_analysis:
    input:
        R_data = rules.blastn.output.R_data
    output:
        Rplot = resultpath+"RDATA/{barcode}_barplot.png"
    params:
        path = resultpath
    shell:
        """

        if [ ! -d {params.path}RDATA ];then
            mkdir {params.path}RDATA
	    fi
        if [ ! -d {params.path}R_RESULT ];then
            mkdir {params.path}R_RESULT
	    fi
        Rscript script/HBV_analysis.R {input} {params.path}
        """        
