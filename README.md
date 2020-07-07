# MinION_HBV

# Quick presentation

This pipeline is designed to analyse viral data from nanopore sequencing. It was primary developed to analyse HBV data, but was extended also to HDV. However, this pipeline can be adapted for other virus aswell with custom preparations specific to the virus.

# Workflow

The pipeline take as input fastq from demultiplexed nanopore data.
For each barcode, here is the global workflow:
Step1 => merging fastq
Step2 => filtering merged fastq
Step3 => blast research on each reference given
Step4 => build consensus sequence on the majoritary reference

# Requirements & Tools

This pipeline use several tools.
Instalation and use of these tool is managed by conda and snakemake.

# Quick using steps

Step 1 : Get and install Anaconda here if needed => https://www.anaconda.com/products/individual <=  
Step 2 : make sure snakemake is installed on your computer. You can quikly create a new conda environment with snakemake by using this command:  
```
conda create -c bioconda -c conda-forge -n snakemake snakemake-minimal
```
Step 3 : download latest version of the pipeline using git command:  
```
git clone https://github.com/Stygiophobic/VIRiONT.git
```
Step 4 : launch the pipeline by executing:  
```
conda activate snakemake
cd VIRiONT
nohup ./launch_snakefile.sh > Report_Analysis.txt & 
```

# Pipeline ouputs

All pipeline output are stored in the choosen path (see below).
For each analysis, the following output are produced:
MERGED folder, containing merged fastq.
TRIMMED folder, containing trimmed fastq by NanoFilt with the given parameters.
REFSEQ folder , containing all the splited references from the main reference fasta given in input.
BD folder, containing the database built using blast.
FASTA folder, containing fasta converted fastq using seqtk. Useless?
BLASTN_RESULT folder, containing the blastn results for each fastq.
BLASTN_ANALYSIS folder, containing the blastn analysis using R, including the list of reads matching with the best reference, with barplot of reference repartition for each barcode.
FILTERED folder, containing fastq with read corresponding with best reference only using seqkit.
BAM folder, containing sorted and indexed bam from each filtered fastq, and mapped on their respecting best matching reference.
COVERAGE folder, containing coverage table from each bam using bedtools.
VCF / VCF_NORM / VCF_FILTER folders, containing variant calling normalized and filtered results using bcftools. To improve!
CONS folder, containing consensus sequence generated for each fastq using bcftools.


# Input and configuration

Here is a view on the parameters to check before launching analysis. To change parameters, you have to open MyREF/XXXXX.sh with a text editor. All parameters are located in the ###### CONFIGURATION ####### section.

data_loc : here is the path where fastq data are stored. Be sure this path leads on all barcode folders you want to analyse. Currently, only fastq repositories marked as "barcode*" are interpreted as repository data. If needed, rename your rep as "barcode*"
result_loc : path leading to the output folders produced by the analysis. NB: snakemake will recursively create the path, so a previous mkdir is unnecessary.
ref_loc : path to the file containing all references sequences used for the blastn analysis. We you need to create a new one, check examples in ref/ folder.
ref_table : path to the Table (currently a csv is required) containing reference list for the blastn analysis. /!\ Column must correspond to the fasta name in ref_loc for correctly computing the best reference, and each row must correspond to each fasta header in the reference fasta file. Please check analysis/table_analysis if needed.
thread_number : Define number of threads to use for the analysis.