#setwd("C:/Users/regueex/Desktop/VIRiONT_R")

argv <- commandArgs(TRUE)

BM_file<-as.character(argv[1])
half_summ<-as.character(argv[2])
output_summ_table<-as.character(argv[3])

half_summ_table<-read.csv2(half_summ)
##############################################

BM_read_table<-read.table(BM_file,header=F)
colnames(BM_read_table)<-c("sequence","sample","reference")
BM_read_table$length<-nchar(as.character(BM_read_table$sequence))
BM_read_table$count<-1
BM_read_table<-BM_read_table[,-1]

bc_ref_assign<-unique(BM_read_table[,c("sample","reference")])


nbread_BM<-setNames(aggregate(BM_read_table$count,list(BM_read_table$sample),
                                  function(x) sum(x)),c("sample","readnumber"))
min_BM<-setNames(aggregate(BM_read_table$length,list(BM_read_table$sample),
                               function(x) min(x)),c("sample","minlength"))
max_BM<-setNames(aggregate(BM_read_table$length,list(BM_read_table$sample),
                               function(x) max(x)),c("sample","maxlength"))
mean_BM<-setNames(aggregate(BM_read_table$length,list(BM_read_table$sample),
                                function(x) mean(x)),c("sample","meanlength"))
median_BM<-setNames(aggregate(BM_read_table$length,list(BM_read_table$sample),
                                  function(x) median(x)),c("sample","medianlength"))

BM_metric<-merge(nbread_BM,min_BM,by="sample")
BM_metric<-merge(BM_metric,max_BM,by="sample")
BM_metric<-merge(BM_metric,mean_BM,by="sample")
BM_metric<-merge(BM_metric,median_BM,by="sample")
BM_metric<-merge(BM_metric,bc_ref_assign,by="sample")
colnames(BM_metric)[7]<-"status"

full_summ_metric<-rbind(half_summ_table,BM_metric)


write.csv2(full_summ_metric,output_summ_table,row.names = F)
