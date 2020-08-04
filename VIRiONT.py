#!usr/bin/en python3
import os
import glob

configfile : "config.yaml"

datapath=config['PathToData']
resultpath=config['PathToResult']
refpath=config['PathToReference']
analysis_table=config['AnalysisTable']
trim_min=config['Lmin']
trim_max=config['Lmax']

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
rule pipeline_ending:
    input:
        #merged_file = expand(resultpath+'MERGED/{barcode}_merged.fastq',barcode=BARCODE),
        #trimmed_file = expand(resultpath+'TRIMMED/{barcode}_trimmed.fastq',barcode=BARCODE),
        #human_bam = expand(resultpath+"DEHOSTING/{barcode}_human.bam",barcode=BARCODE),
        #viral_bam = expand(resultpath+"VIRAL/{barcode}_viral.bam",barcode=BARCODE),
        #viral_fastq = expand(resultpath + 'VIRAL/{barcode}_viral.fastq',barcode=BARCODE),
        #converted_fastq = expand(resultpath+"FASTA/{barcode}.fasta", barcode=BARCODE),
        #ref_rep=resultpath+"REFSEQ/",
        #database = expand(resultpath+"DB/"+database_name+".{ext}", ext=["nhr", "nin", "nsq"]),
        #R_data = expand(resultpath+"BLASTN_RESULT/{barcode}_fmt.txt" ,barcode=BARCODE),
        #best_ref = expand(resultpath+"BLASTN_ANALYSIS/{barcode}_bestref.txt",barcode=BARCODE),
        #merged_filtered = expand(resultpath+"FILTERED/{barcode}_bestref.fastq" ,barcode=BARCODE),
        #conv_filtered_fastq = expand(resultpath+"METRIC/{barcode}_readseq.fasta" ,barcode=BARCODE),
        #spliced_bam = expand(resultpath+"BAM/{barcode}_spliced.bam"  ,barcode=BARCODE),
        #sorted_bam = expand(resultpath+"BAM/{barcode}_sorted.bam" ,barcode=BARCODE),
        #coverage = expand(resultpath+"COVERAGE/{barcode}.cov" ,barcode=BARCODE), 
        cov_plot = resultpath+"COVERAGE/cov_plot.pdf" ,
        fasta_cons = expand(resultpath+"CONSENSUS/{barcode}_cons.fasta",barcode=BARCODE), 
        #QC
        metric_sum = resultpath+"QC_ANALYSIS/POST-ASSIGN_METRICS/metric_summary.tsv",
        metric_sum_dehost = resultpath+"QC_ANALYSIS/DEHOSTING/metric_summary.tsv" ,
        fastqc_results = expand(resultpath+"QC_ANALYSIS/FASTQ_RAW/{barcode}/",barcode=BARCODE), 
        flagstat = expand(resultpath+"QC_ANALYSIS/DEHOSTING/{barcode}_human.txt",barcode=BARCODE), 


rule merge:
    message:
        "Merging fastq: path/to/data/*.fastq ==> path/to/results/MERGED/{barcode}_merged.fastq "
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

rule fastqc:
    input:
        viral_fastq = rules.merge.output.merged_fastq
    output:
        fastqc_results = directory(resultpath+"QC_ANALYSIS/FASTQ_RAW/{barcode}/")
    conda:
        "env/fastqc.yaml"
    shell:
        """
        fastqc --outdir {output.fastqc_results} -f fastq {input.viral_fastq}
        """

rule get_hg19:
    message:
        "download if necessary the hg19 reference genome."
    output:
        hg19_ref =  "ref/hg19.fa"
    shell:
        """
        wget -P ref/ http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
        gunzip ref/hg19.fa.gz       
        """
rule index_hg19:
    input:
        hg19 = rules.get_hg19.output.hg19_ref
    output:
        hg19_index = "ref/hg19.mmi"
    conda:
        "env/minimap2.yaml"
    shell:
        "minimap2 -d {output} {input}"

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
        "NanoFilt --quality 10 --length {trim_min} --maxlength {trim_max} {input} > {output} "

rule mapp_dehost:
    input:
        trimmed_fastq = rules.trimming.output.trimmed_fastq ,
        ref_file= rules.index_hg19.output.hg19_index
    output:
        human_bam = temp(resultpath+"DEHOSTING/{barcode}_human.bam")
    conda:
        "env/minimap2.yaml" 
    threads: 4
    shell:
        """
        minimap2 -t {threads} -ax splice {input.ref_file} {input.trimmed_fastq}  | samtools view -b > {output.human_bam}
        """    

rule samtools_flagstat:
    input:
        human_bam = rules.mapp_dehost.output.human_bam
    output:
        flagstat = resultpath+"QC_ANALYSIS/DEHOSTING/{barcode}_human.txt"
    conda:
        "env/samtools.yaml"
    shell:
        "samtools flagstat {input.human_bam}  > {output.flagstat}"

rule nonhuman_extract:
    input:
        human_bam = rules.mapp_dehost.output.human_bam
    output:
        nonhuman_bam = resultpath+"NONHUMAN/{barcode}_nonhuman.bam"
    conda:
        "env/samtools.yaml" 
    shell: 
        "samtools view -b -f 4 {input.human_bam} > {output.nonhuman_bam} "   

rule bamtofastq:
    input:
        nonhuman_bam = rules.nonhuman_extract.output.nonhuman_bam 
    output:
        nonhuman_fastq = resultpath + 'NONHUMAN/{barcode}_nonhuman.fastq',
    conda:
        "env/bedtools.yaml"
    shell:
        """
        bedtools bamtofastq  -i {input.nonhuman_bam} -fq {output.nonhuman_fastq} 
        """

rule converting:
    message:
        "Converting fastq==>fasta using seqkt for blastn research"
    input:
        nonhuman_fastq = rules.bamtofastq.output.nonhuman_fastq
    output:
        converted_fastq = temp(resultpath+"FASTA/{barcode}.fasta" ),
        converted_fastq_temp = temp(resultpath+"FASTA/{barcode}_temp.fasta" )
 
    conda:
        "env/seqtk.yaml"        
    shell:
        """
        seqtk seq -A {input.nonhuman_fastq} > {output.converted_fastq}
        sed '/^>/d' {output.converted_fastq} > {output.converted_fastq_temp}
        """    

rule compute_metrics_dehosting_R:
    input:
        fasta_file = rules.converting.output.converted_fastq_temp,
    output:
        metric = temp(resultpath+"QC_ANALYSIS/DEHOSTING/{barcode}.tsv")
    conda:
        "env/Renv.yaml"   
    shell:
        """
        Rscript script/read_metrics.R {input.fasta_file} {wildcards.barcode} "DEHOSTED" {output.metric}
        """
rule concat_metrics:
    input:
        metric = expand(rules.compute_metrics_dehosting_R.output.metric,barcode=BARCODE)
    output:
        metric_sum = resultpath+"QC_ANALYSIS/DEHOSTING/metric_summary.tsv"
    run:    
        listfile=glob.glob(resultpath + "QC_ANALYSIS/DEHOSTING/*.tsv")
        sumfile=open(resultpath+"QC_ANALYSIS/DEHOSTING/metric_summary.tsv",'w')
        sumfile.write("SAMPLE\tSTATUS\tNB_READ\tMEAN\tMEDIAN\n")
        for file in listfile:
            readfile=open(file,'r')
            for line in readfile:
                sumfile.write(line+"\n")
            readfile.close()
        sumfile.close()

rule split_reference:
    message:
        "Spliting reference for isolate each genotype sequence."
    input:
        ref_file= refpath
    output:
        ref_rep=directory(resultpath+"REFSEQ/")
    shell:
        "script/split_reference.py {input} {output} "

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

rule blastn_ref:
    message:
        "Blasting {barcode}.fasta on the custom database."
    input: 
        fasta_file = rules.converting.output.converted_fastq,
        database = rules.make_db.output.database ,
    output:
        R_data = resultpath+"BLASTN_ANALYSIS/{barcode}_fmt.txt" 
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
        analyse=database_name
    conda:
        "env/Renv.yaml"          
    shell:
        """
        Rscript script/Blastn_analysis.R {input.R_data} \
            {input.AnalTable} {params.analyse} \
            {output.read_list} {output.best_ref} {output.ref_count_plot}
        """  

rule extract_read_from_merge:
    message:
        "filtrate read from majoritary reference using seqkit."
    input:
        read_list = rules.blastn_analysis.output.read_list,
        nonhuman_fastq = rules.bamtofastq.output.nonhuman_fastq
    output:
        merged_filtered = resultpath+"FILTERED/{barcode}_bestref.fastq" 
    conda:
        "env/seqkit.yaml"          
    shell:
        """
        seqkit grep --pattern-file {input.read_list} {input.nonhuman_fastq} > {output}
        """           

rule converting_filtered:
    message:
        "Converting fastq==>fasta using seqkt for blastn research"
    input:
        merged_filtered = rules.extract_read_from_merge.output.merged_filtered ,
    output:
        converted_fastq = temp(resultpath+"METRIC/{barcode}_readseq.fasta") ,
        converted_fastq_temp = temp(resultpath+"METRIC/{barcode}_filtered.fasta") ,
    conda:
        "env/seqtk.yaml"        
    shell:
        """
        seqtk seq -A {input.merged_filtered} > {output.converted_fastq_temp}
        sed '/^>/d' {output.converted_fastq_temp} > {output.converted_fastq}
        """  

rule compute_metrics_R:
    input:
        converted_fastq = rules.converting_filtered.output.converted_fastq ,
        best_ref = rules.blastn_analysis.output.best_ref ,
        split_ref_path = rules.split_reference.output.ref_rep ,
    output:
        metric = temp(resultpath+"QC_ANALYSIS/POST-ASSIGN_METRICS/{barcode}.tsv")
    conda:
        "env/Renv.yaml"   
    shell:
        """
        bestref=`cat {input.best_ref}`
        Rscript script/read_metrics.R {input.converted_fastq} {wildcards.barcode} $bestref {output.metric}
        """    

rule concat_metrics_dehost:
    input:
        metric = expand(rules.compute_metrics_R.output.metric,barcode=BARCODE)
    output:
        metric_sum = resultpath+"QC_ANALYSIS/POST-ASSIGN_METRICS/metric_summary.tsv"
    run:    
        listfile=glob.glob(resultpath + "QC_ANALYSIS/POST-ASSIGN_METRICS/*.tsv")
        sumfile=open(resultpath+"QC_ANALYSIS/POST-ASSIGN_METRICS/metric_summary.tsv",'w')
        sumfile.write("SAMPLE\tREFERENCE\tNB_READ\tMEAN\tMEDIAN\n")
        for file in listfile:
            readfile=open(file,'r')
            for line in readfile:
                sumfile.write(line+"\n")
            readfile.close()
        sumfile.close()
        
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

rule alignqc:
    input:
        split_ref_path = rules.split_reference.output.ref_rep ,
        sorted_bam = rules.sort_index.output.sorted_bam ,
        best_ref = rules.blastn_analysis.output.best_ref ,
    output:
        alignqc_results = directory(resultpath+"QC_ANALYSIS/ALIGNQC/{barcode}/")
    conda:
        "env/alignqc.yaml"
    shell:
        """
        rm -r {output.alignqc_results} #prevent alignqc to fail du to an already existing directory
        bestref=`cat {input.best_ref}`
        alignqc analyze {input.sorted_bam} -g {input.split_ref_path}${{bestref}}.fasta \
            --no_transcriptome  --output_folder {output.alignqc_results}
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
        """
        bedtools genomecov -ibam {input} -d | sed 's/$/ {wildcards.barcode}/' > {output.coverage} 
        """
rule plot_cov:
    input:
        cov = expand(rules.coverage.output.coverage,barcode=BARCODE)
    output:
        cov_sum = temp(resultpath+"COVERAGE/cov_sum.cov") ,
        cov_plot = resultpath+"COVERAGE/cov_plot.pdf"
    conda:
        "env/Renv.yaml" 
    shell:
        """
        cat {input.cov} > {output.cov_sum}
        Rscript script/plot_cov.R {output.cov_sum} {output.cov_plot}
        """


rule bam_mpileup:
    input:
        split_ref_path = rules.split_reference.output.ref_rep ,
        sorted_bam = rules.sort_index.output.sorted_bam ,
        best_ref = rules.blastn_analysis.output.best_ref ,
    output:
        vcf = resultpath+"VCF/{barcode}_sammpileup.vcf"
    conda:
        "env/samtools.yaml"  
    shell:
        """
        bestref=`cat {input.best_ref}`
        samtools mpileup -d 20000 -f {input.split_ref_path}${{bestref}}.fasta {input.sorted_bam} -Q 7 > {output}
        """             

rule script_varcaller:
    input:
        vcf = rules.bam_mpileup.output.vcf
    output:
        fasta_cons = resultpath+"CONSENSUS/{barcode}_cons.fasta"
    shell:
        """
        perl script/pathogen_varcaller_MINION.PL {input} 0.5 {output} 
        """     