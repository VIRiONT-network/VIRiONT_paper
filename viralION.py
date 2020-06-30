#!usr/bin/en python3
import os
import glob

configfile : "config.yaml"

datapath=config['PathToData']
resultpath=config['PathToResult']
refpath=config['PathToReference']
analysis=config['AnalysisName']
analysis_table=config['AnalysisTable']

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
        best_ref = expand(resultpath+"BLASTN_ANALYSIS/{barcode}_bestref.txt",barcode=BARCODE),
        merged_filtered = expand(resultpath+"FILTERED/{barcode}_bestref.fastq" ,barcode=BARCODE),
        spliced_bam = expand(resultpath+"BAM/{barcode}_spliced.bam"  ,barcode=BARCODE),
        sorted_bam = expand(resultpath+"BAM/{barcode}_sorted.bam" ,barcode=BARCODE),
        coverage = expand(resultpath+"COVERAGE/{barcode}.cov" ,barcode=BARCODE), 
        vcf_bcftools = expand(resultpath+"VCF/{barcode}.vcf" ,barcode=BARCODE), 
        vcf_norm = expand(resultpath+"VCF_NORM/{barcode}_norm.vcf"  ,barcode=BARCODE), 
        vcf_filter = expand(resultpath+"VCF_FILTER/{barcode}_filter.vcf.gz"  ,barcode=BARCODE),      
        cons = expand(resultpath+"CONS/{barcode}.fasta" ,barcode=BARCODE),    


#concatenate all fastq files 
rule merge:
    message:
        "Merging fastq: {barcode}/*.fastq ==> MERGED/{barcode}_merged.fastq "
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

#trimming
rule trimming:
    message:
        "Filtering fastq using NanoFilt."
    input:
        merged_fastq = rules.merge.output.merged_fastq
    output:
        trimmed_fastq = resultpath+"TRIMMED/{barcode}_trimmed.fastq"
    conda:
        "env/nanofilt.yaml"
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
    conda:
        "env/seqtk.yaml"        
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
rule make_db:
    message:
        "build blast database from reference using blastn."
    input:
        ref_fasta_file = refpath 
    output:
        database = expand(resultpath+"DB/"+database_name+".{ext}", ext=["nhr", "nin", "nsq"])
    params:
        database_path = resultpath+"DB/"+database_name
    conda:
        "env/blast.yaml"                  
    shell:
        """
        makeblastdb -in {input.ref_fasta_file} -out {params.database_path} -input_type fasta -dbtype nucl
        """

#execute blastn
rule blastn_ref:
    message:
        "Blasting {barcode}.fasta on the custom database."
    input: 
        fasta_file = rules.converting.output.converted_fastq,
        database = rules.make_db.output.database ,
    output:
        R_data = resultpath+"BLASTN_RESULT/{barcode}_fmt.txt" 
    threads: 4
    params:
        database_path = resultpath+"DB/"+database_name    
    conda:
        "env/blast.yaml"               
    shell:
        """
        blastn -db {params.database_path} -query {input.fasta_file} -outfmt 6 -out {output.R_data} -num_threads {threads}
        """   

rule blastn_analysis:
    message:
        "computing the majoritary reference using R."
    input:
        R_data = rules.blastn_ref.output.R_data ,
        AnalTable = analysis_table ,
    output:
        read_list = resultpath+"BLASTN_ANALYSIS/{barcode}_read-list.txt",
        best_ref = resultpath+"BLASTN_ANALYSIS/{barcode}_bestref.txt",
        ref_count_plot = resultpath+"BLASTN_ANALYSIS/{barcode}_barplot.png"
    params:
        anal=analysis
    conda:
        "env/seqtk.yaml"          
    shell:
        """
        Rscript script/Blastn_analysis.R {input.R_data} \
            {input.AnalTable} {params.anal} \
            {output.read_list} {output.best_ref} {output.ref_count_plot}
        """  

rule extract_read_from_merge:
    message:
        "filtrate read from majoritary reference using seqkit."
    input:
        read_list = rules.blastn_analysis.output.read_list,
        trim_fastq = rules.trimming.output.trimmed_fastq
    output:
        merged_filtered = resultpath+"FILTERED/{barcode}_bestref.fastq" 
    conda:
        "env/seqkit.yaml"          
    shell:
        """
        seqkit grep --pattern-file {input.read_list} {input.trim_fastq} > {output}
        """           

rule alignemnt:
    message:
        "alignment on the majoritary reference using minimap2."
    input:   
        best_ref = rules.blastn_analysis.output.best_ref ,
        merged_filtered = rules.extract_read_from_merge.output.merged_filtered ,
        split_ref_path = rules.split_reference.output.ref_rep ,
    output:
        spliced_bam = resultpath+"BAM/{barcode}_spliced.bam"  
    conda:
        "env/minimap2.yaml"              
    shell:
        """
        bestref=`cat {input.best_ref}`
        minimap2 -ax splice {input.split_ref_path}${{bestref}}.fasta {input.merged_filtered} > {output.spliced_bam}
        """  
rule sort_index:
    message:
        "sorting and indexing bam files using samtools."
    input:
        bam = rules.alignemnt.output.spliced_bam ,
    output:
        sorted_bam = resultpath+"BAM/{barcode}_sorted.bam" ,
        index_file = resultpath+"BAM/{barcode}_sorted.bam.bai"
    conda:
        "env/samtools.yaml"          
    shell:
        """
        samtools sort {input.bam} > {output.sorted_bam}
        samtools index {output.sorted_bam}
        """        
rule coverage:
    message:
        "compute coverage from bam using bedtools."
    input:
        sorted_bam = rules.sort_index.output.sorted_bam
    output:
        coverage = resultpath+"COVERAGE/{barcode}.cov" 
    conda:
        "env/bedtools.yaml"          
    shell:
        "bedtools genomecov -ibam {input} -d > {output}"

rule variant_calling:
    message:
        "compute genetic likehood and variant calling using bcftools."
    input:
        split_ref_path = rules.split_reference.output.ref_rep ,
        best_ref = rules.blastn_analysis.output.best_ref ,
        sorted_bam = rules.sort_index.output.sorted_bam
    output:
        vcf_bcftools = resultpath+"VCF/{barcode}.vcf" 
    conda:
        "env/bcftools.yaml"          
    shell:
        """
        bestref=`cat {input.best_ref}`
        bcftools mpileup -f {input.split_ref_path}${{bestref}}.fasta {input.sorted_bam} | bcftools call -mv -Ov -o {output.vcf_bcftools} 
        """   

rule VC_norm:
    message:
        "normalize indel from vcf using bcftools."
    input:
        split_ref_path = rules.split_reference.output.ref_rep ,
        best_ref = rules.blastn_analysis.output.best_ref ,
        vcf = rules.variant_calling.output.vcf_bcftools
    output:
        vcf_norm = resultpath+"VCF_NORM/{barcode}_norm.vcf" 
    conda:
        "env/bcftools.yaml"             
    shell:
        """
        bestref=`cat {input.best_ref}`
        bcftools norm -f {input.split_ref_path}${{bestref}}.fasta {input.vcf} -Ov -o {output.vcf_norm}
        """

rule VC_filter:
    message:
        "filter adjacent indel from vcf using bcftools."
    input:
        vcf_norm = rules.VC_norm.output.vcf_norm
    output:
        vcf_filter = resultpath+"VCF_FILTER/{barcode}_filter.vcf.gz"   
    conda:
        "env/bcftools.yaml"           
    shell:
        "bcftools filter --IndelGap 5 {input.vcf_norm} -Oz -o {output.vcf_filter}"   

rule create_cons:
    message:
        "generating consensus from filter.vcf using bcftools."
    input:
        split_ref_path = rules.split_reference.output.ref_rep ,
        best_ref = rules.blastn_analysis.output.best_ref ,
        vcf_filter = rules.VC_filter.output.vcf_filter
    output:
        cons = resultpath+"CONS/{barcode}.fasta" 
    conda:
        "env/bcftools.yaml"         
    shell:
        """
        bestref=`cat {input.best_ref}`
        bcftools index {input.vcf_filter}
        cat {input.split_ref_path}${{bestref}}.fasta | bcftools consensus {input.vcf_filter} > {output.cons}
        """       


