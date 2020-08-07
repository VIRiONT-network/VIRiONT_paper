argv <- commandArgs(TRUE)

blastn_result<-as.character(argv[1])
reference<-as.character(argv[2])
readlist<-as.character(argv[3])

blastn_result<-read.table(blastn_result,sep="\t",header=T)

subtable<-subset(blastn_result,GENOTYPE==reference)
list_fastq<-as.character(subtable$READ)
write(list_fastq,readlist)
