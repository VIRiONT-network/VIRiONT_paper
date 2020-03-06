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
        read_list = expand(resultpath+"R_RESULT/{barcode}_list.txt",barcode=BARCODE),
        merged_filtered = expand(resultpath+"FILTER/{barcode}_bestgeno.fastq",barcode=BARCODE), 
        spliced_data = expand(resultpath+"BAM/{barcode}_spliced.bam" ,barcode=BARCODE),
        sorted_bam =  expand(resultpath+"BAM/{barcode}_spliced_sorted.bam" ,barcode=BARCODE), 
        depth_file = expand(resultpath+"DEPTH/{barcode}_spliced_sorted_indexed.depth"  ,barcode=BARCODE), 
        vcf = expand(resultpath+"VCF/{barcode}_varcall"  ,barcode=BARCODE),
        fasta_cons = expand(resultpath+"FINAL_OUTPUT/{barcode}_cons.fasta",barcode=BARCODE),

    params:
        trimmomatic =  "tool/Trimmomatic-0.39/trimmomatic-0.39.jar" ,
        trim_param = config['TrimmParam'] ,
        seqkit = "tool/seqkit/seqkit",
        minimap = "tool/minimap2/minimap2",
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
        read_list = resultpath+"R_RESULT/{barcode}_list.txt",
        best_geno = resultpath+"R_RESULT/{barcode}_bestgeno.txt"
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

rule extract_read_from_merge:
    input:
        read_list = rules.R_HBV_analysis.output.read_list,
        trim_fastq = rules.trimming.output.trimmed_fastq
    output:
        merged_filtered = resultpath+"FILTER/{barcode}_bestgeno.fastq" 
    shell:
        """
        {rules.all.params.seqkit} grep --pattern-file {input.read_list} {input.trim_fastq} > {output}
        """              
rule alignemnt_splice:
    input:   
        best_geno = rules.R_HBV_analysis.output.best_geno ,
        merged_filtered = rules.extract_read_from_merge.output.merged_filtered
    output:
        spliced_bam = resultpath+"BAM/{barcode}_spliced.bam"  
    shell:
        """
        bestgeno=`cat {input.best_geno}`
        {rules.all.params.minimap} -ax splice ref/${{bestgeno}}.fa {input.merged_filtered} > {output.spliced_bam}
        """              

rule bam_sorting:
    input:
        spliced_bam = rules.alignemnt_splice.output.spliced_bam
    output:
        sorted_bam =  resultpath+"BAM/{barcode}_spliced_sorted.bam"     
    shell:
        """
        samtools sort {input} > {output}
        """  

rule bam_indexing:
    input:
        sorted_bam = rules.bam_sorting.output.sorted_bam
    output:
        depth_file = resultpath+"DEPTH/{barcode}_spliced_sorted_indexed.depth"   
    shell:
        """
        samtools index {input}
        samtools depth -m 200000 {input} > {output}
        """                  

rule bam_mpileup:
    input:
        sorted_bam = rules.bam_sorting.output.sorted_bam   ,
        best_geno = rules.R_HBV_analysis.output.best_geno ,
        depth_file = rules.bam_indexing.output.depth_file
    output:
        vcf = resultpath+"VCF/{barcode}_varcall"
    shell:
        """
        bestgeno=`cat {input.best_geno}`
        pathref=`echo "ref/$bestgeno.fa"`
        samtools mpileup -d 200000 -f $pathref {input.sorted_bam} > {output}
        """             

rule script_varcaller:
    input:
        vcf = rules.bam_mpileup.output.vcf
    output:
        fasta_cons = resultpath+"FINAL_OUTPUT/{barcode}_cons.fasta"
    shell:
        """
        perl script/pathogen_varcaller_MINION.PL {input} 0.5 {output}
        """                