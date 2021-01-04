#!usr/bin/en python3
import os
import glob
import pandas as pd

configfile : "config/config.yaml"

datapath=config['PathToData']
if (datapath[-1] != "/"):
	datapath=datapath+"/"
resultpath=config['PathToResult']
if (resultpath[-1] != "/"):
	resultpath=resultpath+"/"
refpath=config['PathToReference']
trim_min=config['Lmin']
trim_max=config['Lmax']
trim_head=config['headcrop']
trim_tail=config['tailcrop']
quality_read=config['quality']
MI_cutoff=config['multiinf']


#NanoFilt Command construction with input parameters
if (trim_min>0):
	lengmin="--length "+str(trim_min)
else:
	lengmin=""
if (trim_max>0):
	lengmax="--maxlength "+str(trim_max)
else:
	lengmax=""
if (trim_head>0):
	head="--headcrop "+str(trim_head)
else:
	head=""
if (trim_tail>0):
	tail="--tailcrop "+str(trim_tail)
else:
	tail=""
if (quality_read>0):
	qfiltering="--quality "+str(quality_read)
else:
	qfiltering=""


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

#final output
rule pipeline_ending:
	input:
		######## INTERMEDIATE FILES ########
		#database = expand(resultpath+"00_SUPDATA/DB/"+database_name+".{ext}", ext=["nhr", "nin", "nsq"]),
		#merged_file = expand(resultpath+'01_MERGED/{barcode}_merged.fastq',barcode=BARCODE),
		#human_bam = expand(resultpath+"02_DEHOSTING/{barcode}_human.bam",barcode=BARCODE),
		#meta_fastq = expand(resultpath + '02_DEHOSTING/{barcode}_meta.fastq',barcode=BARCODE),
		#trimmed_file = expand(resultpath+'03_FILTERED_TRIMMED/{barcode}_filtered_trimmed.fastq',barcode=BARCODE),
		#converted_fastq = expand(resultpath+"FASTA/{barcode}.fasta", barcode=BARCODE),
		#ref_rep=resultpath+"00_SUPDATA/REFSEQ/",
		#R_data = expand(resultpath+"04_BLASTN_ANALYSIS/{barcode}_fmt.txt" ,barcode=BARCODE),
		#merged_data = resultpath+"04_BLASTN_ANALYSIS/ALL_refcount.tsv"
		######## FINAL OUTPUTS ########		
		#blastn_result=expand(resultpath+"04_BLASTN_ANALYSIS/{barcode}_blastnR.tsv",barcode=BARCODE),
		summ_multiinf = resultpath+"04_BLASTN_ANALYSIS/SUMMARY_Multi_Infection.tsv",

rule merging_fastq:
	message:
		"Merging fastq/fastq.gz files from /datapath/{wildcards.barcode}/ to the /outputpath/01_MERGED/ folder."
	input: 
		barcode_rep = datapath+"{barcode}/",
	output: 
		merged_fastq = resultpath+"01_MERGED/{barcode}_merged.fastq"  
	shell: 
		"""
		gz_file_num=`find {input.barcode_rep} -name "*.gz" | wc -l`
		if [ $gz_file_num -gt 0 ]
		then
        	zcat {input.barcode_rep}* > {output}
		else
        	cat {input.barcode_rep}* > {output}
		fi
		#cat {input.barcode_rep}* > {output}
		"""

rule get_hg19:
	message:
		"download if necessary the hg19 reference genome. Executed once per VIRiONT installation."
	output:
		hg19_ref =  "ref/hg19.fa"
	shell:
		"""
		wget -P ref/ http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
		gunzip ref/hg19.fa.gz       
		"""

rule index_hg19:
	message:
		"Indexing the hg19 reference genome for a quicker dehosting. Executed once per VIRiONT installation."
	input:
		hg19 = rules.get_hg19.output.hg19_ref
	output:
		hg19_index = "ref/hg19.mmi"
	conda:
		"env/minimap2.yaml"
	shell:
		"minimap2 -d {output} {input}"

rule hg19_dehosting:
	message:
		"Aligning /outputpath/01_MERGED/{wildcards.barcode}_merged.fastq on human genome for identifying host reads using minimap2."
	input:
		merged_fastq = rules.merging_fastq.output.merged_fastq ,
		ref_file= rules.index_hg19.output.hg19_index
	output:
		human_bam = temp(resultpath+"02_DEHOSTING/{barcode}_human.bam")
	conda:
		"env/minimap2.yaml" 
	threads: 4
	shell:
		"""
		minimap2 -t {threads} -ax splice {input.ref_file} {input.merged_fastq}  | samtools view -b > {output.human_bam}
		"""    

rule nonhuman_read_extract:
	message:
		"Extract unaligned reads from {wildcards.barcode}_human.bam using samtools."
	input:
		human_bam = rules.hg19_dehosting.output.human_bam
	output:
		nonhuman_bam = temp(resultpath+"02_DEHOSTING/{barcode}_meta.bam")
	conda:
		"env/samtools.yaml" 
	shell: 
		"samtools view -b -f 4 {input.human_bam} > {output.nonhuman_bam} "   

rule converting_bam_fastq:
	message:
		"Convert {wildcards.barcode}_nonhuman.bam to {wildcards.barcode}_nonhuman.fastq using bedtools."
	input:
		nonhuman_bam = rules.nonhuman_read_extract.output.nonhuman_bam 
	output:
		nonhuman_fastq = resultpath + '02_DEHOSTING/{barcode}_meta.fastq',
	conda:
		"env/bedtools.yaml"
	shell:
		"""
		bedtools bamtofastq  -i {input.nonhuman_bam} -fq {output.nonhuman_fastq} 
		"""

rule trimming_fastq:
	message:
		"Filtering and trimming {wildcards.barcode}_meta.fastq using NanoFilt with input parameters."
	input:
		meta_fastq = rules.converting_bam_fastq.output.nonhuman_fastq
	output:
		trimmed_fastq = resultpath+"03_FILTERED_TRIMMED/{barcode}_filtered_trimmed.fastq"
	conda:
		"env/nanofilt.yaml"
	shell:
		"NanoFilt {lengmin} {lengmax} {head} {tail} {qfiltering} {input.meta_fastq} > {output.trimmed_fastq} "

rule converting_fastq_fasta:
	message:
		"Converting {wildcards.barcode}_meta.fastq into {wildcards.barcode}.fasta for blastn research using seqtk."
	input:
		nonhuman_fastq = rules.trimming_fastq.output.trimmed_fastq
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

rule split_reference:
	message:
		"Spliting the input reference for isolate each genotype sequence if multiple reference are given."
	input:
		ref_file= refpath
	output:
		ref_rep=directory(resultpath+"00_SUPDATA/REFSEQ/")
	shell:
		"script/split_reference.py {input} {output} "

rule make_db:
	message:
		"Build blast database from input reference using makeblastdb."
	input:
		ref_fasta_file = refpath 
	output:
		database = expand(resultpath+"00_SUPDATA/DB/"+database_name+".{ext}", ext=["nhr", "nin", "nsq"])
	params:
		database_path = resultpath+"00_SUPDATA/DB/"+database_name
	conda:
		"env/blast.yaml"                  
	shell:
		"""
		makeblastdb -in {input.ref_fasta_file} -out {params.database_path} -input_type fasta -dbtype nucl
		"""

rule blastn_ref:
	message:
		"Blasting {wildcards.barcode}.fasta on the custom database."
	input: 
		fasta_file = rules.converting_fastq_fasta.output.converted_fastq,
		database = rules.make_db.output.database ,
	output:
		R_data = resultpath+"04_BLASTN_ANALYSIS/{barcode}_fmt.txt" 
	threads: 4
	params:
		database_path = resultpath+"00_SUPDATA/DB/"+database_name    
	conda:
		"env/blast.yaml"               
	shell:
		"""
		blastn -db {params.database_path} -query {input.fasta_file} \
			-outfmt "6 qseqid sseqid bitscore slen qlen length pident" \
			-out {output.R_data} -num_threads {threads}
		"""   

rule count_refmatching:
	input:
		R_data = rules.blastn_ref.output.R_data ,
		ref_table = rules.split_reference.output.ref_rep,
	output:
		ref_count = temp(resultpath+"04_BLASTN_ANALYSIS/{barcode}_refcount.tsv"),
		blastn_result = resultpath+"04_BLASTN_ANALYSIS/{barcode}_blastnR.tsv"
	conda:
		"env/Renv.yaml"     
	shell:
		"""
		Rscript script/count_ref.R {input.R_data} \
			{input.ref_table}R_table_analysis.csv \
			{wildcards.barcode} \
			{output.ref_count} \
			{output.blastn_result}
		"""  

rule MI_analysis:
	input:
		count_ref_data = expand(rules.count_refmatching.output.ref_count,barcode=BARCODE)
	output:
		merged_data = temp(resultpath+"04_BLASTN_ANALYSIS/ALL_refcount.tsv"),
		plot_pdf = resultpath+"04_BLASTN_ANALYSIS/read_repartition.pdf",
		summ_multiinf = resultpath+"04_BLASTN_ANALYSIS/SUMMARY_Multi_Infection.tsv"
	conda:
		"env/Renv.yaml"   
	shell:
		"""
		cat {input.count_ref_data} > {output.merged_data}
		Rscript script/MI_analysis.R {output.merged_data} \
			{MI_cutoff} \
			{output.plot_pdf} \
			{output.summ_multiinf}
		"""