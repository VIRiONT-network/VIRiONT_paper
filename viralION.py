#!usr/bin/en python3
import os
import glob

configfile : "config.yaml"

datapath=config['PathToData']
resultpath=config['PathToResult']
refpath=config['PathToReference']

#get database name
filename=os.path.basename(refpath)
list_split=filename.split(".")
database_name=list_split[0]


#get all barcodes in a list after demultiplexing
barcode_list = glob.glob(datapath+"barcode*")
BARCODE=[]
for BC in barcode_list:
    barcode=BC[-9:]
    BARCODE.append(barcode)

#final output
rule all:
    input:
        merged_file = expand(resultpath+'MERGED/{barcode}_merged.fastq',barcode=BARCODE),
        trimmed_file = expand(resultpath+'TRIMMED/{barcode}_trimmed.fastq',barcode=BARCODE),
        converted_fastq = expand(resultpath+"FASTA/{barcode}.fasta", barcode=BARCODE),
        ref_rep=resultpath+"REFSEQ/",
        database = expand(resultpath+"DB/"+database_name+".{ext}", ext=["nhr", "nin", "nsq"]),
        R_data = expand(resultpath+"BLASTN_RESULT/{barcode}_fmt.txt" ,barcode=BARCODE),


#concatenate all fastq files 
rule merge:
    message:
        "Merging fastq: {barcode}/*.fastq ==> MERGED/{barcode}_merged.fastq "
    input: 
        lambda wildcards: expand(resultpath+"{barcode}", barcode=BARCODE)
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
    message:
        "Filtering fastq using NanoFilt."
    input:
        merged_fastq = rules.merge.output.merged_fastq
    output:
        trimmed_fastq = resultpath+"TRIMMED/{barcode}_trimmed.fastq"
   
    shell:
        "NanoFilt --quality 10 --length 100 --maxlength 1500 {input} > {output} "

#convert fastq into fasta
rule converting:
    message:
        "Converting {barcode}_trimmed.fastq ==> {barcode}.fasta using seqkt"
    input:
        trimmed_fastq = rules.trimming.output.trimmed_fastq
    output:
        converted_fastq = resultpath+"FASTA/{barcode}.fasta" 
    shell:
        """
        seqtk seq -A {input} > {output}
        """    

#Separate each fasta from ref into each unique fasta sequence.
rule split_reference:
    message:
        "Spliting reference for isolate each genotype sequence."
    input:
        ref_file= refpath
    output:
        ref_rep=resultpath+"REFSEQ/"
    shell:
        "script/split_reference.py {input} {output} "

#build the blast database
rule make_db_HBV:
    message:
        "build blast database from reference using blastn."
    input:
        ref_fasta_file = refpath 
    output:
        database = expand(resultpath+"DB/"+database_name+".{ext}", ext=["nhr", "nin", "nsq"])
    params:
        database_path = resultpath+"DB/"+database_name        
    shell:
        """
        makeblastdb -in {input.ref_fasta_file} -out {params.database_path} -input_type fasta -dbtype nucl
        """

#execute blastn
rule blastn_geno:
    message:
        "Blasting {barcode}.fasta on the custom database."
    input: 
        fasta_file = rules.converting.output.converted_fastq,
        database = rules.make_db_HBV.output.database ,
    output:
        R_data = resultpath+"BLASTN_RESULT/{barcode}_fmt.txt" ,
    params:
        database_path = resultpath+"DB/"+database_name         
    shell:
        """
        blastn -db {params.database_path} -query {input.fasta_file} -outfmt 6 -out {output.R_data}
        """   