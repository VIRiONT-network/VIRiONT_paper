#!/bin/bash
path=$1
#Merge each fastas from each barcode rep
for barcode in `ls DATASET`
do
    cat DATASET/${barcode}/*.fastq > MERGE/${barcode}_merged.fastq
done 


##################################################################
#pour faire un blast en local									#
#pour générer une banque de données							#
#compiler toutes les séquences de référence par GT dans un		# 
#même fichier fasta											#
##################################################################

# Basecalling #

guppy_basecaller --flowcell FLO-MIN106 --kit SQK-PBK004 --input_path /chemin/fast5 --save_path /chemin/Guppy_basecall_output --cpu_threads_per_caller nombre de coeur --
num_callers 1 --records_per_fastq 0 --disable_pings --recursive --qscore_filtering --min_qscore 7

# Demultiplexing #

guppy_barcoder  -i /chemin/Guppy_basecall_output/Merge_fastq -s /chemin/Guppy_demultiplexing_output-merge --barcode_kits  SQK-PBK004 -t 8 --recursive --trim_barcodes  #Blastnlocal

#se placer dans le fichier du code barre BCn

# cat *.fastq> merged.fastq

#longueur de reads à éliminer e.g. <x 500 ou 1500
Trimmomatic SE -threads 8 -phred33 merged.fastq output.fastq MINLEN:x 


java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 MERGED/BC01_merged.fastq HadrienR/TRIMMED/BC01_trimmed.fastq MINLEN:x 

#transformation des fichier fastq à blaster en fichier fasta

seqtk seq -A TxBCn.fastq>TxBCn.fasta

#blast local 

blastn -db Ref(tous les gT).fasta -query TxBC11.fasta -outfmt 6 -out TxBC11fmt6.txt  #=> à analyser via R


#Alignement	
												
#Alignement Nanopore
minimap2 -ax map-ont P3_GTE.fasta TmergedBCn.fastq>SpliceBC06.bam

#Alignement pour slices
minimap2 -ax splice P3_GTE.fasta TmergedBCn.fastq>SpliceBC06.bam

samtools sort SpliceBC06.bam>SSpliceBC06.bam

samtools index SSpliceBC06.bam

samtools depth -m 200000 SSpliceBC06.bam>SSpliceBC06.depth


#essais pour Fleur
cutadapt -a A{100} -o output.fastq mergedBC06.fastq

grep -c "AAAAAAAAAAA" input.fastq
#nb de lignes (reads) avec le motif AAAAAAAAAAA

grep -B1 -A2 "AAAAAAAAAAAAA" input.fastq | sed '/--/d' > output.fastq

more input.fastq #permet d’afficher les premières lignes d’un fichier

Ancre de Fleur pour 3’RACE : 5’CGATACGCTACGTAACGGCATGAC 3’

#consensusJosset

samtools mpileup -d 200000 -f P3_GTE.fasta ST500BC02.bam> varcalerBC02T500

perl pathogen_varcaler_MINION.PL varcalerBC02T500 0.5 output.fasta

essai
minimap2 -ax splice Ref-D.fasta mergedBC06.fastq|samtools sort 


 


