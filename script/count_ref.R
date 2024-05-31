argv <- commandArgs(TRUE)

#parse arguments
blastn_data<-read.table(argv[1],header=F,sep="\t")
colnames(blastn_data)<-c("qseqid","sseqid","bitscore","slen","qlen","length","pindent")
blastn_data<-blastn_data[,c(1:3)]
colnames(blastn_data)<-c("READ","GENOTYPE","SCORE")

analysis_table<-read.csv2(argv[2])
samplename<-as.character(argv[3])
output_ref_count<-as.character(argv[4])
output_blastR<-as.character(argv[5])
bitscore_min<-as.numeric(argv[6])

blastn_data<-subset(blastn_data,SCORE>=bitscore_min)


if(nrow(blastn_data)==0){
    write.table(blastn_data,output_ref_count,col.names = F,row.names = F,quote = F,sep="\t")
    f <- file(output_blastR, open="w")
    truncate(f)
    quit(save = "no")
}

#get the bitscore max for each read
MAX_SCORE<-aggregate(blastn_data$SCORE,by=list(blastn_data$READ), function(x) max(x))
colnames(MAX_SCORE)<-c("READ","MAX")

#get reference matching for the bitscore max, for each read
TABLE_BEST_GENO<-merge(blastn_data,MAX_SCORE,by = "READ")
TABLE_BEST_GENO<-subset(TABLE_BEST_GENO,SCORE==MAX)

write.table(TABLE_BEST_GENO,output_blastR,sep="\t",row.names=F,quote=F)

TABLE_BEST_GENO$count<-1

#reference matching count
COUNT_GENO<-aggregate(TABLE_BEST_GENO$count,by=list(TABLE_BEST_GENO$GENOTYPE),function(x) sum(x))
colnames(COUNT_GENO)<-c("reference","count")
COUNT_GENO$sample<-samplename

write.table(COUNT_GENO,output_ref_count,col.names = F,row.names = F,quote = F,sep="\t")
