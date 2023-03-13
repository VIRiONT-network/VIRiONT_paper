# VIRiONT (VIRal in-house Oxford Nanopore Technologies sequencing) Pipeline

# UPDATE (2023-03-08)
The following repository was dedicated to the work below:
>Charre C, Regue H, Dény P, Josset L, Chemin I, Zoulim F, Scholtes C. Improved Hepatitis Delta Virus genome characterization by single molecule full-length genome sequencing combined with VIRiONT pipeline. J Med Virol. 2023 Mar 6.  
>[https://doi.org/10.1002/jmv.28634](https://doi.org/10.1002/jmv.28634)

# Quick presentation

VIRiONT pipeline is designed to analyze data from amplicon-approach-based Nanopore sequencing. It was primarily developed to analyze data from hepatitis B virus (HBV) complete genome long read sequencing after an amplification step of the full genome in one fragment (about 3.2kb in lenght). It was secondarily extended to hepatitis delta virus (HDV, about 1.7kb in lenght). However, this pipeline can be adapted for other viruses and pathogens as well. You just need to upload a custom reference dataset specific to the pathogen and adapted to the amplicon design used prior your ONT sequencing.  
Of note, if you plan to use this workflow to analyze shorter or longer amplicons, be sure to choose and adapt parameters of the read lenght filtering steps, see below.

# Workflow

The pipeline takes as input uncompressed or compressed fastq from demultiplexed nanopore data by Guppy, usually stored like this:  
```
barcoding/barcode01/*.fastq
barcoding/barcode02/*.fastq
...
barcoding/barcode24/*.fastq
``` 

For each barcode, you can find herein the global workflow:  
**Step1** => removing human reads from fastq files namely the dehosting step.  
**Step2** => trimming fastq using given parameters (primer removal and read lenght filtering).  
**Step3** => blastn analysis leads to the selection of the best matching reference(s) among the uploaded custom reference dataset, based on a best bitscore mapping read count.  
**Step4** => generation of a final consensus sequence using a custom perl script with a tunable minimal variant frequency :  
(i) a first Minimap2 (option splice) alignment guided by the best matching reference (selected at the blastn step) leads to an intermediate pre-consensus sequence.  
(ii) This latter sequence is used as the sample’s own mapping reference for a second realignment of the reads that enables to generate to the final consensus sequence. Of note, consensus sequence can be called at a tunable minimum depth set up per default at 20X (usually applied among Nanopore community).  
**Step5** => generation of a phylognenetic tree including consensus and reference sequences of the custom dataset using a Maximum-likelihood statistical method (1000 bootstrap replicates) after a MUSCLE-based-alignment.  

![image info](./documents/Figure_Workflow_VIRiONT.png)

# Requirements & Tools

Installation and use of the tools required by VIRiONT pipeline is fully managed by conda and snakemake.  
Environment files and software versions are available in the *env/* folder.  

*filtering and trimming* : NanoFilt v2.7.1. See => [NanoFilt](https://github.com/wdecoster/nanofilt) <=  
*dehosting and mapping* : minimap2 v2.17. See => [minimap2](https://github.com/lh3/minimap2) <=  
*bam management and pileup* : samtools v1.12. See => [samtools](http://samtools.github.io/) <=  
*blastn analysis* : blast v2.5.0. See => [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) <=  
*fastq and fasta management* : seqkit v0.12.1. See => [seqkit](https://github.com/shenwei356/seqkit) <=  
*fastq and fasta management* : seqtk v1.3. See => [seqtk](https://github.com/lh3/seqtk) <=   
*coverage computation* : bedtools v2.29.2. See => [bedtools](https://github.com/arq5x/bedtools) <=   
*statistics and plotting* : R v4.0.3. See => [R](https://cran.r-project.org/) <=  
*consensus and reference sequence alignment before phylognetic tree generation* : muscle v3.8.1551. See => [muscle](http://www.drive5.com/muscle) <=  
*newick tree computation* : iqtree v2.0.3. See => [iqtree](http://www.iqtree.org/) <=  
*tree drawing* : ete3 v3.1.2. See => [ete3](http://etetoolkit.org/) <=  
*read correction* : CANU v2.1.1. See => [CANU](https://canu.readthedocs.io/en/latest/) <=  


# Multiple Infection case

VIRiONT can detect co-infection by several genotypes or sub-genotypes for instance or contamination case when processing to the blastn analysis step.  
Based on the best bitscore value, we proceed to a readcount matching for each reference sequence.   
The best-matching reference is used to percent normalize each readcount.  
A reference is selected if the percent normalized count is equal or above the input tunable threshold.  
The selected references are used to independently run the pipeline and produce as many consensus sequences as selected references for each barcode.  

*note: if you're not confident with this option, it is posssible to set the **MI_cutoff** to 100 for getting only the major best-matching reference.*

# Few words about CANU and read correction  

The error rate associated with ONT technologies remain one of the main issues. However several long read correction methods are now available to significantly decrease the error rate. Here, we choose to use for VIRiONT the [CANU assembler](https://canu.readthedocs.io/en/latest/), which can also correct reads instead of de novo assembly.  
If a read correction is envisaged for your data, please keep in mind these informations:  
-Dont hesitate to visit the official documentation to learn more about CANU [https://canu.readthedocs.io/en/latest/](https://canu.readthedocs.io/en/latest/).  
-The "maxInputCoverage" parameter is tunable with VIRiONT. This parameter subsample data if too low, and significantly increase computing time if to high. Depending on your read load and computing capacities, choose this parameter carefully.  
-The "correctedErrorRate" parameter is also tunable with VIRiONT. Change this parameter can impact read correction and increase/decrease computing time.  
-Read correction with VIRiONT is done after read assignation on the the major best-matching references selected.  
-The correction step redirect all resources (threads and memory) on CANU.  
-The corrected reads from CANU are in a fasta format output. Due to a conversion fasta=>fastq , the basequality associated information is lost during the conversion process. VIRiONT set default quality to "?" (Phred score 30 in the ASCI table). Some variant calling parameters should become useless.  

# Concerning the mutation research module  

Currently, this step is strictly optionnal in the analysis, because data table are specificly adapted to our work, using specific amplicon design and custom references.  
Thus, we do not recommand using this sub-analysis even if dealing with HBV data, unless you are able to adapt table and mutation screening script to your experimental design.   
If you are interested in this feature, results should be published soon.  

# Quick using steps

**Step 1 :** Get and install Anaconda herein if needed.  
More informations about conda here => https://www.anaconda.com/products/individual <=  

```
wget https://repo.anaconda.com/archive/Anaconda3-2020.07-Linux-x86_64.sh # download the linux installing script.
chmod +x Anaconda3-2020.07-Linux-x86_64.sh #give rights to execute the script.
./Anaconda3-2020.07-Linux-x86_64.sh #Simply follow screen instructions for installing conda.
#Check if conda is installed. You'll probably need to refresh/open a new terminal for starting using conda.
conda -h #print conda commands.
conda -V #print installed conda version.
```

**Step 2 :** Make sure snakemake is installed on your computer.  

Snakemake 3.9.0 version or above is required for conda interaction.  
pandas python library is also required for dataframe reading.  
You can quickly create a new conda environment with the pandas library and the latest available snakemake version by using this command (currently the 5.20.1 version).  
All test and pipeline implementation have been done under the snakemake 5.20.1 version. (Changes or errors can occur when using versions above):  

```
conda create -c bioconda -c conda-forge -n VIRiONT_env snakemake-minimal=5.20.1 pandas
```

**Step 3 :** Download latest version of the pipeline using git command:  

```
git clone https://github.com/VIRiONT-network/VIRiONT.git
```

**Step 4 :** After setting your parameters in the launch_VIRiONT.sh script launch the pipeline by executing:  

```
conda activate VIRiONT_env #only if you previously created this environment for VIRiONT use.  
cd VIRiONT  
./launch_VIRiONT.sh  
```

# Input and configuration

Herein is a general overview of the tunable parameters to set before launching analysis. Currently, to change parameters, you have to open  *config/params.yaml* with a text editor.  
All parameters are located in the ###### CONFIGURATION ####### section.  

**GENERAL PARAMETERS:**  
**data_loc** : path where fastq data are stored. Be sure this path leads on all barcode folders you want to analyze. Currently, only fastq repositories marked as "barcode*" are interpreted as repository data. If needed, rename your rep as "barcode*". "barcode12" as an example.  
**result_loc** : path leading to the output folders produced by the analysis. NB: the pipeline will recursively create the path, so a previous mkdir is unnecessary.  
**ref_loc** : path to the custom reference sequence dataset fasta file used by the pipeline especially for the blastn analysis. If you need to create a new one, check examples in *ref/* folder.  
**thread_number** : define number of threads to be used for the analysis.  
**mem_cost** : define the memory amount in mb to be used for the analysis.  

**TRIMMING/FILTERING PARAMETERS:**  
**min_length** : minimal read lenght required for passing the filtering step.  
**max_length** : maximal read lenght required for passing the filtering step.  
**min_qual_ONT** : minimal read quality for passing the filtering step. Based on the ONT quality score.  
**head_trim** : number of nucleotides to be trimmed at the 5'end (forward primer removal).  
**tail_trim** : number of nucleotides to be trimmed at the 3'end (reverse primer removal).  

**READ CORRECTION**  
If you want to correct your read, please check the CANU documentation for more information: => [documentation](https://canu.readthedocs.io/)
**correction** : enable read correction. TRUE/FALSE to enable/disable correction.   
**cov_correction** : average coverage for read correction. Reads above this coverage will be removed from the analysis. Corresponding to the CANU's "maxInputCoverage" parameter.  
**error_rate** : CANU parameter for ONT read correction. Corresponding to the CANU's "correctedErrorRate" parameter.     

**VARIANT CALLING PARAMETERS:**  
**max_depth** : maximum depth for samtools mpileup (sub-sampling of the total amount of reads is possible for a faster analysis).  
**basequality** : base quality cutoff for samtools mpileup (default parameter of samtools mpileup is set up at 13).  

**CONSENSUS PARAMETERS:**   
**Vfreq** : minor variant frequency threshold to call the consensus sequence. You can generate consensus sequences through a simple 50% majority rule or with a lower minor variant frequency to get degenerated bases.  
**min_cov** : minimum coverage expected to generate consensus sequence. N is called if the depth at the position is below this threshold. Of note, a minium depth of 20x is recommended by the Nanopore community.

**MULTI-INFECTION PARAMETER:**  
**MI_cutoff** : Percentage cutoff for detecting a multiple infection case. See more above in the **Multiple Infection case** section.  

**HBV MUTATION RESEARCH:**  
Optionnal part.  
**HBV_mut** : set this paramater to TRUE will initiate optionnal mutation research for specitif HBV dataset. FALSE to disable.     
**path_table** : location of the mutation table.  
**freq_min** : variant frequency under the threshold will be removed from filtered table. Should be between 0 and 100.  
**window_pos** : If using reference without primers, correct position in the mutation table. Could be Positive or negative.  


# Pipeline ouputs

All VIRiONT pipeline outputs are stored in the choosen path (see above) indicated in the VIRiONT_MI.sh script, in the form of the following folders and files:  

**param_file.txt** : a text file, summarizing parameters used for your current analysis.   
**00_SUPDATA/DB** : a folder, containing the database built by makeblastdb from your uploaded custom dataset of reference sequences fasta file.  
**00_SUPDATA/REFSEQ** : a folder, containing the individual reference sequences.  
**01_MERGED** : a folder, containing merged demultiplexed PASS fastq files.  
**02_DEHOSTING** : a folder, containing dehosted (non human) fastq files.  
**03_FILTERED_TRIMMED** : a folder, containing filtered and trimmed fastq by NanoFilt with the given parameters.  
**04_BLASTN_ANALYSIS** : a folder, containing the blastn analysis results, including the list of reads matching with the best reference and a barplot of reference repartition for each barcode, available as a pdf file.  
**05_REFILTERED_FASTQ** : a folder, containing fastq files with only reads mapping best-matching references.   
**06_PRECONSENSUS** : a folder, containing intermediate data related to the generated pre-consensus sequence.  
**07_BAM** : a folder, containing bam files generated from the alignment guided by the pre-consensus sequence.  
**08_VCF** : a folder, containing mpileup and variant calling files.  
**09_CONSENSUS** : a folder, containing consensus sequences in fasta files generated by a custom perl script.  
**10_QC_ANALYSIS** : a folder, containing useful read metrics after the different steps of filtering.  
**11_PHYLOGENETIC_TREE** : a folder, containing phylogenetic tree of consensus sequences with associated matched reference sequences.  
**12_COVERAGE** : a folder, containing coverage table from each bam file generated using bedtools as well as coverage plots for each sample compiled into one pdf file.  

