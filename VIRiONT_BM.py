#!usr/bin/en python3
import os
import glob

configfile : "config.yaml"

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
variant_frequency=config['variantfrequency']
mincov_cons=config['mincov']

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
        #merged_file = expand(resultpath+''01_MERGED/{barcode}_merged.fastq',barcode=BARCODE),
        #trimmed_file = expand(resultpath+'02_TRIMMED/{barcode}_trimmed.fastq',barcode=BARCODE),
        #human_bam = expand(resultpath+"03_DEHOSTING/{barcode}_human.bam",barcode=BARCODE),
		#viral_fastq = expand(resultpath + '03_DEHOSTING/{barcode}_nonhuman.fastq',barcode=BARCODE),
        #converted_fastq = expand(resultpath+"FASTA/{barcode}.fasta", barcode=BARCODE),
        #ref_rep=resultpath+"00_SUPDATA/REFSEQ/",
        #database = expand(resultpath+"00_SUPDATA/DB/"+database_name+".{ext}", ext=["nhr", "nin", "nsq"]),
        #R_data = expand(resultpath+"04_BLASTN_ANALYSIS/{barcode}_fmt.txt" ,barcode=BARCODE),
        #best_ref = expand(resultpath+"04_BLASTN_ANALYSIS/{barcode}_bestref.txt",barcode=BARCODE),
        #merged_filtered = expand(resultpath+"FILTERED/{barcode}_bestref.fastq" ,barcode=BARCODE),
        #conv_filtered_fastq = expand(resultpath+"METRIC/{barcode}_readseq.fasta" ,barcode=BARCODE),
        #spliced_bam = expand(resultpath+"BAM/{barcode}_spliced.bam"  ,barcode=BARCODE),
        #sorted_bam = expand(resultpath+"BAM/{barcode}_sorted.bam" ,barcode=BARCODE),
        ######## COVERAGE ANALYSIS ########  
        #coverage = expand(resultpath+"COVERAGE/{barcode}.cov" ,barcode=BARCODE), 
        cov_plot = resultpath+"07_COVERAGE/cov_plot.pdf" ,
        ######## CONSENSUS FILES ########   
        #fasta_cons = expand(resultpath+"09_CONSENSUS/{barcode}_cons.fasta",barcode=BARCODE), 
        cons_seq = resultpath+"09_CONSENSUS/all_cons.fasta",
        ######## QC METRIC FILES ######## 
        summ_table = resultpath+"10_QC_ANALYSIS/METRIC_summary_table.csv",
        ######## QC PHYLOGENETIC TREE FILES ########  
        #allseq = resultpath+"11_PHYLOGENETIC_TREE/allseq.fasta",
        #align_seq = resultpath+"11_PHYLOGENETIC_TREE/align_seq.fasta",
        #NWK_tree = resultpath+"11_PHYLOGENETIC_TREE/IQtree_analysis.treefile",
        tree_pdf = resultpath+"11_PHYLOGENETIC_TREE/RADIAL_tree.pdf"

rule merging_fastq:
    message:
        "Merging fastq into the {wildcards.barcode}/ folder if needed. "
    input: 
        lambda wildcards: expand(datapath+"{barcode}", barcode=BARCODE)
    output: 
        merged_fastq = resultpath+"01_MERGED/{barcode}_merged.fastq"  
    params:
        path = datapath
    shell: 
        """
        cat {params.path}{wildcards.barcode}/* > {output}
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

rule trimming_fastq:
    message:
        "Filtering and trimming {wildcards.barcode}_merged.fastq using NanoFilt with input parameters."
    input:
        merged_fastq = rules.merging_fastq.output.merged_fastq
    output:
        trimmed_fastq = resultpath+"02_TRIMMED/{barcode}_trimmed.fastq"
    conda:
        "env/nanofilt.yaml"
    shell:
        "NanoFilt --quality 10 --length {trim_min} --maxlength {trim_max} {input} > {output} "

rule hg19_dehosting:
    message:
        "Aligning {wildcards.barcode}_trimmed.fastq on human genome for identifying host reads using minimap2."
    input:
        trimmed_fastq = rules.trimming_fastq.output.trimmed_fastq ,
        ref_file= rules.index_hg19.output.hg19_index
    output:
        human_bam = temp(resultpath+"03_DEHOSTING/{barcode}_human.bam")
    conda:
        "env/minimap2.yaml" 
    threads: 4
    shell:
        """
        minimap2 -t {threads} -ax splice {input.ref_file} {input.trimmed_fastq}  | samtools view -b > {output.human_bam}
        """    

rule nonhuman_read_extract:
    message:
        "Extract unaligned reads from {wildcards.barcode}_human.bam using samtools."
    input:
        human_bam = rules.hg19_dehosting.output.human_bam
    output:
        nonhuman_bam = resultpath+"03_DEHOSTING/{barcode}_nonhuman.bam"
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
        nonhuman_fastq = resultpath + '03_DEHOSTING/{barcode}_nonhuman.fastq',
    conda:
        "env/bedtools.yaml"
    shell:
        """
        bedtools bamtofastq  -i {input.nonhuman_bam} -fq {output.nonhuman_fastq} 
        """

rule converting_fastq_fasta:
    message:
        "Converting {wildcards.barcode}_nonhuman.fastq into {wildcards.barcode}.fasta for blastn research using seqtk."
    input:
        nonhuman_fastq = rules.converting_bam_fastq.output.nonhuman_fastq
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
        blastn -db {params.database_path} -query {input.fasta_file} -outfmt 6 -out {output.R_data} -num_threads {threads}
        """   

rule blastn_analysis:
    message:
        "Computing the majoritary reference for {wildcards.barcode}_nonhuman.fastq."
    input:
        R_data = rules.blastn_ref.output.R_data ,
        ref_table = rules.split_reference.output.ref_rep,
    output:
        read_list = resultpath+"04_BLASTN_ANALYSIS/{barcode}_read-list.txt",
        best_ref = resultpath+"04_BLASTN_ANALYSIS/{barcode}_bestref.txt",
        ref_count_plot = resultpath+"04_BLASTN_ANALYSIS/{barcode}_count_barplot.png",
        ref_ratio_plot = resultpath+"04_BLASTN_ANALYSIS/{barcode}_ratio_plot.png",
    params:
        analyse=database_name
    conda:
        "env/Renv.yaml"          
    shell:
        """
        Rscript script/Blastn_analysis.R {input.R_data} \
            {input.ref_table}R_table_analysis.csv  \
            {output.read_list} {output.best_ref} {output.ref_ratio_plot} \
            {output.ref_count_plot}
        """  

rule extract_matching_read:
    message:
        "Filtrate read from {wildcards.barcode}_nonhuman.fastq into {wildcards.barcode}_filtered.fastq using seqkit."
    input:
        read_list = rules.blastn_analysis.output.read_list,
        nonhuman_fastq = rules.converting_bam_fastq.output.nonhuman_fastq
       
    output:
        merged_filtered = resultpath+"05_REFFILTERED/{barcode}_bestref.fastq" 
    conda:
        "env/seqkit.yaml"          
    shell:
        """
        seqkit grep --pattern-file {input.read_list} {input.nonhuman_fastq} > {output}
        """           

rule filtered_fastq_alignemnt:
    message:
        "Align {wildcards.barcode}_filtered.fastq on the bestreference.fasta using minimap2."
    input:   
        best_ref = rules.blastn_analysis.output.best_ref ,
        merged_filtered = rules.extract_matching_read.output.merged_filtered ,
        split_ref_path = rules.split_reference.output.ref_rep ,
    output:
        spliced_bam = resultpath+"06_BAM/{barcode}_sorted.bam" 
    conda:
        "env/minimap2.yaml"              
    shell:
        """
        bestref=`cat {input.best_ref}`
        minimap2 -ax splice {input.split_ref_path}${{bestref}}.fasta {input.merged_filtered} | samtools sort > {output.spliced_bam}
        samtools index {output.spliced_bam}  
        """  

rule compute_coverage:
    message:
        "Compute coverage from {wildcards.barcode}_sorted.bam using bedtools."
    input:
        sorted_bam = rules.filtered_fastq_alignemnt.output.spliced_bam
    output:
        coverage = resultpath+"07_COVERAGE/{barcode}.cov" 
    conda:
        "env/bedtools.yaml"          
    shell:
        """
        bedtools genomecov -ibam {input} -d -split | sed 's/$/ {wildcards.barcode}/' > {output.coverage} 
        """

rule plot_coverage:
    message:
        "Plot coverage summary using R."
    input:
        cov = expand(rules.compute_coverage.output.coverage,barcode=BARCODE)
    output:
        cov_sum = temp(resultpath+"07_COVERAGE/cov_sum.cov") ,
        cov_plot = resultpath+"07_COVERAGE/cov_plot.pdf"
    conda:
        "env/Renv.yaml" 
    shell:
        """
        cat {input.cov} > {output.cov_sum}
        Rscript script/plot_cov.R {output.cov_sum} {output.cov_plot}
        """

rule variant_calling:
    message:
        "Variant calling on {wildcards.barcode}_sorted.bam aligned bam using samtools."
    input:
        split_ref_path = rules.split_reference.output.ref_rep ,
        sorted_bam = rules.filtered_fastq_alignemnt.output.spliced_bam ,
        best_ref = rules.blastn_analysis.output.best_ref ,
    output:
        vcf = resultpath+"08_VCF/{barcode}_sammpileup.vcf"
    conda:
        "env/samtools.yaml"  
    shell:
        """
        bestref=`cat {input.best_ref}`
        samtools mpileup -d 20000 -f {input.split_ref_path}${{bestref}}.fasta {input.sorted_bam} -Q 7 > {output}
        """             

rule read_metric:
    message:
        "Extract read sequence from  MERGED/TRIMMED/DEHOSTED {wildcards.barcode}.fastq for metric computation."
    input:
        raw_fastq = rules.merging_fastq.output.merged_fastq ,
        trimmed_fastq = rules.trimming_fastq.output.trimmed_fastq ,
        dehosted_fastq = rules.converting_bam_fastq.output.nonhuman_fastq ,
        bestmatched_fastq = rules.extract_matching_read.output.merged_filtered ,
    output:
        raw_read = temp(resultpath+"10_QC_ANALYSIS/{barcode}_rawseq.txt"),
        trimmed_read = temp(resultpath+"10_QC_ANALYSIS/{barcode}_trimmseq.txt"),
        dehosted_read = temp(resultpath+"10_QC_ANALYSIS/{barcode}_dehostseq.txt"),
        bestmatched_read = temp(resultpath+"10_QC_ANALYSIS/{barcode}_BMseq.txt"),
    shell:
        """
        awk '(NR%4==2)' {input.raw_fastq} | sed 's/$/ {wildcards.barcode}/' > {output.raw_read}
        awk '(NR%4==2)' {input.trimmed_fastq} | sed 's/$/ {wildcards.barcode}/'> {output.trimmed_read}
        awk '(NR%4==2)' {input.dehosted_fastq} | sed 's/$/ {wildcards.barcode}/'> {output.dehosted_read}
        awk '(NR%4==2)' {input.bestmatched_fastq} | sed 's/$/ {wildcards.barcode}/' > {output.bestmatched_read}
        """

rule summarize_metric:
    message:
        "Compute metric and write the summary table using R."
    input:
        raw_read = expand(rules.read_metric.output.raw_read ,barcode=BARCODE),
        trimmed_read = expand(rules.read_metric.output.trimmed_read ,barcode=BARCODE),
        dehosted_read = expand(rules.read_metric.output.dehosted_read ,barcode=BARCODE),
        bestmatched_read = expand(rules.read_metric.output.bestmatched_read ,barcode=BARCODE),
    output:
        raw_all = temp(resultpath+"10_QC_ANALYSIS/CAT_rawseq.tab"),
        trimmed_all = temp(resultpath+"10_QC_ANALYSIS/CAT_trimmseq.tab"),
        dehosted_all = temp(resultpath+"10_QC_ANALYSIS/CAT_dehostseq.tab"),
        bestmatched_all = temp(resultpath+"10_QC_ANALYSIS/CAT_BMseq.tab"),
        summ_table = resultpath+"10_QC_ANALYSIS/METRIC_summary_table.csv"
    shell:
        """
        cat {input.raw_read} > {output.raw_all}
        cat {input.trimmed_read} > {output.trimmed_all}
        cat {input.dehosted_read} > {output.dehosted_all}
        cat {input.bestmatched_read} > {output.bestmatched_all}
        Rscript script/summarize_metric.R {output.raw_all} {output.trimmed_all} \
             {output.dehosted_all} {output.bestmatched_all} {output.summ_table}
        """

rule generate_consensus:
    message:
        "Generate consensus sequence from {wildcards.barcode}_sammpileup.vcf using perl script."
    input:
        vcf = rules.variant_calling.output.vcf,    
    output:
        fasta_cons_temp = temp(resultpath+"09_CONSENSUS/{barcode}_cons_temp.fasta") ,
        fasta_cons = resultpath+"09_CONSENSUS/{barcode}_cons.fasta"
    shell:
        """
        perl script/pathogen_varcaller_MINION.PL {input.vcf} {variant_frequency} {output.fasta_cons_temp} {mincov_cons}
        sed  's/>.*/>{wildcards.barcode}_VIRiONT/' {output.fasta_cons_temp} > {output.fasta_cons}
        """    

rule prepareSEQ:
    message:
        "Concatenate all consensus sequences with matched references."
    input:
        allcons = expand(rules.generate_consensus.output.fasta_cons,barcode=BARCODE),
        ref_matched = expand(rules.blastn_analysis.output.best_ref,barcode=BARCODE)
    output:
        filtered_refseq_list = temp(resultpath+"11_PHYLOGENETIC_TREE/matching_ref.txt"),
        filtered_refseq = temp(resultpath+"11_PHYLOGENETIC_TREE/matching_ref.fasta"),
        cons_seq = resultpath+"09_CONSENSUS/all_cons.fasta",
        allseq = temp(resultpath+"11_PHYLOGENETIC_TREE/allseq.fasta"),
    conda:
        "env/seqkit.yaml"   
    shell:
        """
        cat {input.ref_matched} | sort | uniq > {output.filtered_refseq_list}
		seqkit grep --pattern-file {output.filtered_refseq_list} {refpath} > {output.filtered_refseq}
        cat {input.allcons} > {output.cons_seq}
        cat {output.filtered_refseq} {output.cons_seq} > {output.allseq}
        """
    
rule sequenceAlign:
    message:
        "Multiple Alignment on all sequences using Muscle."
    input:
        allseq = rules.prepareSEQ.output.allseq
    output:
        align_seq = temp(resultpath+"11_PHYLOGENETIC_TREE/align_seq.fasta")
    conda:
        "env/muscle.yaml"
    shell:
        "muscle -in {input} -out {output} -maxiters 2"

rule buildTree:
    message:
        "Build Newick tree from alignment using iqtree."
    input:
        align_seq = rules.sequenceAlign.output.align_seq,
    output:
        iqtree = expand(resultpath+"11_PHYLOGENETIC_TREE/IQtree_analysis"+".{ext}", ext=["contree","bionj","ckp.gz","iqtree","log","mldist","splits.nex","treefile"])
    conda:
        "env/iqtree.yaml"
    shell:
        "iqtree -s {input.align_seq} -m K80 -B 1000 -T AUTO --prefix {resultpath}11_PHYLOGENETIC_TREE/IQtree_analysis "

rule plotTree:
    message:
        "Plot radial tree using ETE 3."
    input:
        iqtree = rules.buildTree.output.iqtree ,
        NWK_data = resultpath+"11_PHYLOGENETIC_TREE/IQtree_analysis.treefile"
    output:
        tree_pdf = resultpath+"11_PHYLOGENETIC_TREE/RADIAL_tree.pdf"
    conda:
        "env/ETE3.yaml"
    shell:
        "python3 script/makeTREE.py {input.NWK_data} {output.tree_pdf} "
