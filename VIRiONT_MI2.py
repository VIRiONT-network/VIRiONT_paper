#!usr/bin/en python3
import os
import glob
import pandas as pd

configfile : "config.yaml"

datapath=config['PathToData']
if (datapath[-1] != "/"):
	datapath=datapath+"/"
resultpath=config['PathToResult']
if (resultpath[-1] != "/"):
	resultpath=resultpath+"/"
refpath=config['PathToReference']
analysis_table=config['AnalysisTable']
trim_min=config['Lmin']
trim_max=config['Lmax']
trim_head=config['headcrop']
trim_tail=config['tailcrop']
variant_frequency=config['variantfrequency']


#get database name
filename=os.path.basename(refpath)
list_split=filename.split(".")
database_name=list_split[0]


#get all barcodes in a list after demultiplexing
barcode_list = glob.glob(datapath+"barcode*")
BARCODE=[]
for BC in barcode_list:
	barcode=str(os.path.basename(BC))
	BARCODE.append(barcode)

#BAC A SABLE
data_multiinf = resultpath+"04_BLASTN_ANALYSIS/SUMMARY_Multi_Infection.tsv"
multiinf_table = pd.read_csv(resultpath+"04_BLASTN_ANALYSIS/SUMMARY_Multi_Infection.tsv",sep="\t",names=["barcode", "reference", "ratio"])
sample_list=list(multiinf_table['barcode'])
reference_list=list(multiinf_table['reference'])
assoc_sample_ref_number=len(sample_list)

#produce read list
readlist=[]
for i in range(0,assoc_sample_ref_number):
	readlist.append(resultpath +"READLIST/"+sample_list[i]+"/"+reference_list[i]+"_readlist.txt")

#produce filtered fastq
fastq_filtered=[]
for i in range(0,assoc_sample_ref_number):
	fastq_filtered.append(resultpath +"05_REFFILTERED_FASTQ/"+sample_list[i]+"/"+reference_list[i]+"_filtered.fastq")

#produce bam
bamfile=[]
for i in range(0,assoc_sample_ref_number):
	bamfile.append(resultpath +"06_BAM/"+sample_list[i]+"/"+reference_list[i]+"_sorted.bam")

#produce cov files
covfile=[]
for i in range(0,assoc_sample_ref_number):
	covfile.append(resultpath +"07_COVERAGE/"+sample_list[i]+"/"+reference_list[i]+".cov")

#produce vcf files
vcffile=[]
for i in range(0,assoc_sample_ref_number):
	vcffile.append(resultpath +"08_VCF/"+sample_list[i]+"/"+reference_list[i]+"_sammpileup.vcf")

#produce consensus files
consfile=[]
for i in range(0,assoc_sample_ref_number):
	consfile.append(resultpath +"09_CONSENSUS/"+sample_list[i]+"/"+reference_list[i]+"_cons.fasta")

#produce filterseq files
filtseqfile=[]
for i in range(0,assoc_sample_ref_number):
	filtseqfile.append(resultpath +"10_QC_ANALYSIS/"+sample_list[i]+"/"+reference_list[i]+"_filterseq.txt")


rule all:
    input:
        #readlist,
        #fastq_filtered,
        #bamfile,
        #covfile,
        #vcffile,
        consfile,
        #filtseqfile,
        #raw_read = expand(resultpath+"10_QC_ANALYSIS/{barcode}_rawseq.txt",barcode=BARCODE),
        #trimmed_read = expand(resultpath+"10_QC_ANALYSIS/{barcode}_trimmseq.txt",barcode=BARCODE),
        #dehosted_read = expand(resultpath+"10_QC_ANALYSIS/{barcode}_dehostseq.txt",barcode=BARCODE),
        #raw_all = resultpath+"10_QC_ANALYSIS/CAT_rawseq.tab",
        #trimmed_all = resultpath+"10_QC_ANALYSIS/CAT_trimmseq.tab",
        #dehosted_all = resultpath+"10_QC_ANALYSIS/CAT_dehostseq.tab",
        #bestmatched_all = resultpath+"10_QC_ANALYSIS/CAT_BMseq.tab",
        cov_plot = resultpath+"07_COVERAGE/cov_plot.pdf",
        #summ_table = resultpath+"10_QC_ANALYSIS/METRIC_summary_tableTEMP.csv",
        full_summ_table = resultpath+"10_QC_ANALYSIS/METRIC_summary_table.csv",
        #allseq = resultpath+"10_QC_ANALYSIS/TREE/allseq.fasta",
        #align_seq = resultpath+"10_QC_ANALYSIS/TREE/align_seq.fasta",
        #NWK_tree = resultpath+"10_QC_ANALYSIS/TREE/IQtree_analysis.treefile",
        tree_pdf = resultpath+"10_QC_ANALYSIS/RADIAL_tree.pdf"



rule getfastqlist:
    input:
        blastn_result = resultpath+"04_BLASTN_ANALYSIS/{barcode}_blastnR.tsv"
    output:
        fastqlist = temp(resultpath +"READLIST/{barcode}/{reference}_readlist.txt" )
    shell:
        """
        Rscript script/get_read_list.R {input.blastn_result} {wildcards.reference} {output.fastqlist}
        """

rule extract_matching_read:
	message:
		"filtrate read from majoritary reference."
	input:
		read_list = rules.getfastqlist.output.fastqlist,
		nonhuman_fastq = resultpath + '03_DEHOSTING/{barcode}_nonhuman.fastq'      
	output:
		merged_filtered = resultpath +"05_REFFILTERED_FASTQ/{barcode}/{reference}_filtered.fastq" 
	conda:
		"env/seqkit.yaml"          
	shell:
		"""
		seqkit grep --pattern-file {input.read_list} {input.nonhuman_fastq} > {output}
		"""   

rule filtered_fastq_alignemnt:
    message:
        "alignment on the majoritary reference using minimap2."
    input:   
        merged_filtered = rules.extract_matching_read.output.merged_filtered ,
        split_ref_path = resultpath+"00_SUPDATA/REFSEQ/" ,
    output:
        spliced_bam = resultpath+"06_BAM/{barcode}/{reference}_sorted.bam" 
    conda:
        "env/minimap2.yaml"              
    shell:
        """
        minimap2 -ax splice {input.split_ref_path}{wildcards.reference}.fasta {input.merged_filtered} | samtools sort > {output.spliced_bam}
        samtools index {output.spliced_bam}  
        """

rule compute_coverage:
    message:
        "compute coverage from bam using bedtools."
    input:
        sorted_bam = rules.filtered_fastq_alignemnt.output.spliced_bam
    output:
        coverage = resultpath+"07_COVERAGE/{barcode}/{reference}.cov" 
    conda:
        "env/bedtools.yaml"          
    shell:
        """
        bedtools genomecov -ibam {input} -d -split | sed 's/$/\t{wildcards.barcode}\tMINION/' > {output.coverage} 
        """

rule plot_coverage:
    message:
        "Plot coverage summary."
    input:
        cov = covfile
    output:
        cov_sum = temp(resultpath+"07_COVERAGE/cov_sum.cov") ,
        cov_plot = resultpath+"07_COVERAGE/cov_plot.pdf"
    conda:
        "env/Renv.yaml" 
    shell:
        """
        cat {input.cov} > {output.cov_sum}
        Rscript script/plot_cov_MI.R {output.cov_sum} {output.cov_plot}
        """

rule variant_calling:
    message:
        "Variant calling on reference aligned bam."
    input:
        split_ref_path = resultpath+"00_SUPDATA/REFSEQ/"  ,
        sorted_bam = rules.filtered_fastq_alignemnt.output.spliced_bam ,
    output:
        vcf = resultpath+"08_VCF/{barcode}/{reference}_sammpileup.vcf"
    conda:
        "env/samtools.yaml"  
    shell:
        """
        samtools mpileup -d 20000 -f {input.split_ref_path}{wildcards.reference}.fasta {input.sorted_bam} -Q 7 > {output}
        """   

rule generate_consensus:
    message:
        "Generate consensus sequence from variant calling file."
    input:
        vcf = rules.variant_calling.output.vcf
    output:
        fasta_cons_temp = temp(resultpath+"09_CONSENSUS/{barcode}/{reference}_cons_temp.fasta") ,
        fasta_cons = resultpath+"09_CONSENSUS/{barcode}/{reference}_cons.fasta"
    shell:
        """
        perl script/pathogen_varcaller_MINION.PL {input} {variant_frequency} {output.fasta_cons_temp} 200
        sed  's/>.*/>{wildcards.barcode}_{wildcards.reference}_VIRiONT/' {output.fasta_cons_temp} > {output.fasta_cons}
        """
        
rule read_metric:
    message:
        "Extract read sequence from fastq for metric computation."
    input:
        raw_fastq = resultpath+"01_MERGED/{barcode}_merged.fastq"  ,
        trimmed_fastq = resultpath+"02_TRIMMED/{barcode}_trimmed.fastq" ,
        dehosted_fastq = resultpath + '03_DEHOSTING/{barcode}_nonhuman.fastq' 
    output:
        raw_read = temp(resultpath+"10_QC_ANALYSIS/{barcode}_rawseq.txt"),
        trimmed_read = temp(resultpath+"10_QC_ANALYSIS/{barcode}_trimmseq.txt"),
        dehosted_read = temp(resultpath+"10_QC_ANALYSIS/{barcode}_dehostseq.txt"),
    shell:
        """
        awk '(NR%4==2)' {input.raw_fastq} | sed 's/$/ {wildcards.barcode}/' > {output.raw_read}
        awk '(NR%4==2)' {input.trimmed_fastq} | sed 's/$/ {wildcards.barcode}/' > {output.trimmed_read}
        awk '(NR%4==2)' {input.dehosted_fastq} | sed 's/$/ {wildcards.barcode}/' > {output.dehosted_read}
        """

rule read_metric_MI:
    message:
        "Extract read sequence from fastq for metric computation."
    input:
        bestmatched_fastq = rules.extract_matching_read.output.merged_filtered ,
    output:
        bestmatched_read = temp(resultpath+"10_QC_ANALYSIS/{barcode}/{reference}_filterseq.txt"),
    shell:
        """
        awk '(NR%4==2)' {input.bestmatched_fastq} | sed 's/$/ {wildcards.barcode} {wildcards.reference}/' > {output.bestmatched_read}
        """

rule summarize_metric:
    message:
        "Compute metric and write the summary table."
    input:
        raw_read = expand(rules.read_metric.output.raw_read ,barcode=BARCODE),
        trimmed_read = expand(rules.read_metric.output.trimmed_read,barcode=BARCODE),
        dehosted_read = expand(rules.read_metric.output.dehosted_read,barcode=BARCODE),
        bestmatched_read = expand(filtseqfile),
    output:
        raw_all = temp(resultpath+"10_QC_ANALYSIS/CAT_rawseq.tab"),
        trimmed_all = temp(resultpath+"10_QC_ANALYSIS/CAT_trimmseq.tab"),
        dehosted_all = temp(resultpath+"10_QC_ANALYSIS/CAT_dehostseq.tab"),
        bestmatched_all = temp(resultpath+"10_QC_ANALYSIS/CAT_BMseq.tab"),
        summ_table = temp(resultpath+"10_QC_ANALYSIS/METRIC_summary_tableTEMP.csv"),
        full_summ_table = resultpath+"10_QC_ANALYSIS/METRIC_summary_table.csv",
    shell:
        """
        cat {input.raw_read} > {output.raw_all}
        cat {input.trimmed_read} > {output.trimmed_all}
        cat {input.dehosted_read} > {output.dehosted_all}
        cat {input.bestmatched_read} > {output.bestmatched_all}
        Rscript script/summarize_metric_MI1.R {output.raw_all} {output.trimmed_all} \
             {output.dehosted_all}  {output.summ_table}
        Rscript script/summarize_metric_MI2.R {output.bestmatched_all}\
             {output.summ_table} {output.full_summ_table}
        """

rule prepareSEQ:
    input:
        #matching_refs = multiinf_table,
        allcons = expand(consfile),
        reference_file = refpath
    output:
        filtered_refseq_list = temp(resultpath+"10_QC_ANALYSIS/matching_ref.txt"),
        filtered_refseq = temp(resultpath+"10_QC_ANALYSIS/matching_ref.fasta"),
        cons_seq = resultpath+"09_CONSENSUS/all_cons.fasta",
        allseq = temp(resultpath+"10_QC_ANALYSIS/allseq.fasta"),
    conda:
        "env/seqkit.yaml"   
    shell:
        """
        awk '{{print $2}}' {data_multiinf} | sort | uniq > {output.filtered_refseq_list}
		seqkit grep --pattern-file {output.filtered_refseq_list} {input.reference_file} > {output.filtered_refseq}
        cat {input.allcons} > {output.cons_seq}
        cat {output.filtered_refseq} {output.cons_seq} > {output.allseq}
        """

rule sequenceAlign:
    input:
        allseq = rules.prepareSEQ.output.allseq
    output:
        align_seq = temp(resultpath+"10_QC_ANALYSIS/align_seq.fasta")
    conda:
        "env/muscle.yaml"
    shell:
        "muscle -in {input} -out {output} -maxiters 2"

rule buildTree:
    input:
        align_seq = rules.sequenceAlign.output.align_seq,
    output:
        iqtree = temp(expand(resultpath+"10_QC_ANALYSIS/IQtree_analysis"+".{ext}", ext=["bionj","ckp.gz","iqtree","log","mldist","model.gz","treefile"]))
    conda:
        "env/iqtree.yaml"
    shell:
        " iqtree -s {input.align_seq} --prefix {resultpath}10_QC_ANALYSIS/IQtree_analysis "

rule plotTree:
    input:
        iqtree = rules.buildTree.output.iqtree ,
        NWK_data = resultpath+"10_QC_ANALYSIS/IQtree_analysis.treefile"
    output:
        tree_pdf = resultpath+"10_QC_ANALYSIS/RADIAL_tree.pdf"
    conda:
        "env/ETE3.yaml"
    shell:
        "python3 script/makeTREE.py {input.NWK_data} {output.tree_pdf} "