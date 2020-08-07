# VIRiONT Pipeline (VIRal Oxford Nanopore Technologies sequencing Pipeline)

# Quick presentation

VIRiONT pipeline is designed to analyse data from amplicon approach based Nanopore sequencing. It was primary developed to analyse data from hepatitis B virus (HBV) complete genome long read sequencing after an amplification step. It was secondary extended to hepatitis delta virus (HDV). However, this pipeline can be adapted for other viruses and pathogens as well. You just need a custom reference dataset specific to the pathogen and to the amplicon design used prior your ONT sequencing.

This pipeline was primary designed for studying HBV full genome long read sequencing (about 3kb in lenght). If you plan to use this workflow on shorter or longer amplicons, be sure to chose and adapt parameters used for the filtering step, see below.

# Workflow

The pipeline takes as input uncompressed fastq from demultiplexed nanopore data by Guppy, usually stored like this:  
```
barcoding/barcode01/*.fastq
barcoding/barcode02/*.fastq
...
barcoding/barcode24/*.fastq
``` 
For each barcode, you can find herein the global workflow:  
**Step1** => merging (if several fastq in the barcode rep) and trimming fastq using given parameters.  
**Step2** => removing human read from fastq files.  
**Step3** => blastn research on each reference for get the best maching reference.     
**Step4** => build consensus sequence on the best matching reference after generated alignement files with minimap2 (option splice)    

# workflow dag

![image info](./documents/workflow.png)

# Requirements & Tools

Installation and use of the tools required for VIRiONT is fully managed by conda and snakemake.  
Environment files and software versions are available in the *env/* folder.  

*filtering and trimming* : NanoFilt v2.7.1. See => https://github.com/wdecoster/nanofilt <=  
*dehosting and mapping* : minimap2 v2.17. See => https://github.com/lh3/minimap2 <=  
*bam management* : samtools v1.3.1. See => http://samtools.github.io/ <=  
*blastn analysis* : blast v2.5.0. See => https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download <=  
*fastq and fasta management* : seqkit v0.12.1. See => https://github.com/shenwei356/seqkit <=  
*fastq and fasta management* : seqtk v1.3. See => https://github.com/lh3/seqtk <=   
*coverage computation* : bedtools v2.29.2. See => https://github.com/arq5x/bedtools <=   
*statistics and plotting* : R v4.0.2. See => https://cran.r-project.org/ <=  

# Multiple Infection case

Two workflow are available:  
-One workflow *VIRiONT_BM.py* generating one consensus sequence per barcode, called using the script *VIRiONT_BM.sh*. The major reference in read matching count is used to produce the consensus sequence.  
-The second workflow *VIRiONT_MI.py* called using the script *VIRiONT_MI.sh*. Using an input given threshold, selected references are used to produde as many consensus sequences as selected references for each barcode. A reference is selected when the current_reference_read_count/major_reference*100 is equal or above to the input threshold

# Quick using steps

Step 1 : Get and install Anaconda here if needed => https://www.anaconda.com/products/individual <=  
You may need to start new terminal for starting using conda.   
```
wget https://repo.anaconda.com/archive/Anaconda3-2020.07-Linux-x86_64.sh
chmod +x Anaconda3-2020.07-Linux-x86_64.sh
./Anaconda3-2020.07-Linux-x86_64.sh #Simply follow screen instructions for installing conda.
```
Step 2 : make sure snakemake is installed on your computer. 
```
conda -h
conda -V
``` 
Snakemake 3.9.0 version or above is required for conda interaction.  
pandas python library is also required for dataframe reading.
You can quikly create a new conda environment with the pandas library and the latest available snakemake version by using this command(currently the 5.20.1 version):  
```
conda create -c bioconda -c conda-forge -n VIRiONT_env snakemake-minimal pandas
```
Step 3 : download latest version of the pipeline using git command:  
```
git clone https://github.com/Stygiophobic/VIRiONT.git
```
Step 4 : launch the pipeline by executing:  
```
conda activate VIRiONT_env
cd VIRiONT
./VIRiONT_BM.sh # Generate consensus from best matching reference only 
./VIRiONT_MI.sh # Generate consensus for all references over the given threshold
```

# Pipeline ouputs

All pipeline output are stored in the choosen path (see below).For each analysis, the following output are produced:  

**00_SUPDATA/BD** folder, containing the database built using blast and individual reference sequences.  
**01_MERGED** folder, containing merged fastq.  
**02_TRIMMED** folder, containing trimmed fastq by NanoFilt with the given parameters.  
**03_DEHOSTING** folder, containing dehosted (non human) bam.  
**04_BLASTN_ANALYSIS** folder, containing the blastn analysis using R, including the list of reads matching with the best reference, with barplot of reference repartition for each barcode.  
**05_BESTREF_FILTERED** folder, containing fastq with read corresponding with best reference only using seqkit.  
BAM folder, containing sorted and indexed bam from each filtered fastq, and mapped on their respecting best matching reference.  
**06_BAM** folder, containing bam aligned on the best reference.  
**07_COVERAGE** folder, containing coverage table from each bam using bedtools,and coverage plots compiled in one pdf.  
**08_VCF** foler, containing mpileup and variant calling results.  
**09_CONSENSUS** folder, containing consencus sequences generated by a custom perl script.  
**10_QC_ANALYSIS** folder, containing usefull read metrics at several pipeline steps.  


# Input and configuration

Here is a view on the parameters to check before launching analysis. To change parameters, you have to open VIRiONT/VIRiONT.sh with a text editor. All parameters are located in the ###### CONFIGURATION ####### section.

**data_loc** : path where fastq data are stored. Be sure this path leads on all barcode folders you want to analyse.Currently, only fastq repositories marked as "barcode*" are interpreted as repository data. If needed, rename your rep as "barcode*".  
**result_loc** : path leading to the output folders produced by the analysis. NB: the pipeline will recursively create the path, so a previous mkdir is unnecessary.  
**ref_loc** : path to the file containing all references sequences used for the blastn analysis. We you need to create a new one, check examples in *ref/* folder.  
**ref_table** : path to the Table (currently a csv is required) containing reference list for the blastn analysis. /!\ Column must correspond to the fasta name in ref_loc for correctly computing the best reference, and each row must correspond to each fasta header in the reference fasta file. Please check analysis/table_analysis if needed for an example.  
**min_length** : minimal read size required for passing the filtering step.  
**max_length** : maximal read size required for passing the filtering step.  
**head** : nucleotid length to trimm in 5'.  
**tail** : nucleotid length to trimm in 3'.  
**Vfreq** : minor variant frequency threshold.  
**thread_number** : Define number of threads to use for the analysis.  
**mem_cost** : define the memory amount in mb to use for the analysis.
