#setwd("C:/Users/regueex/Desktop/VIRiONT_R")

argv <- commandArgs(TRUE)

raw_file<-as.character(argv[1])
trimm_file<-as.character(argv[2])
dehost_file<-as.character(argv[3])
BM_file<-as.character(argv[4])
output_summ_table<-as.character(argv[5])



raw_read_table<-read.table(raw_file,header=F)
colnames(raw_read_table)<-c("sequence","sample")
raw_read_table$length<-nchar(as.character(raw_read_table$sequence))
raw_read_table$count<-1
raw_read_table<-raw_read_table[,-1]

nbread_raw<-setNames(aggregate(raw_read_table$count,list(raw_read_table$sample),
                               function(x) sum(x)),c("sample","readnumber"))
min_raw<-setNames(aggregate(raw_read_table$length,list(raw_read_table$sample),
                function(x) min(x)),c("sample","minlength"))
max_raw<-setNames(aggregate(raw_read_table$length,list(raw_read_table$sample),
                            function(x) max(x)),c("sample","maxlength"))
mean_raw<-setNames(aggregate(raw_read_table$length,list(raw_read_table$sample),
                            function(x) mean(x)),c("sample","meanlength"))
median_raw<-setNames(aggregate(raw_read_table$length,list(raw_read_table$sample),
                            function(x) median(x)),c("sample","medianlength"))

raw_metric<-merge(nbread_raw,min_raw,by="sample")
raw_metric<-merge(raw_metric,max_raw,by="sample")
raw_metric<-merge(raw_metric,mean_raw,by="sample")
raw_metric<-merge(raw_metric,median_raw,by="sample")
raw_metric$status<-"RAW"

##############################################

trimm_read_table<-read.table(trimm_file,header=F)
colnames(trimm_read_table)<-c("sequence","sample")
trimm_read_table$length<-nchar(as.character(trimm_read_table$sequence))
trimm_read_table$count<-1
trimm_read_table<-trimm_read_table[,-1]

nbread_trimm<-setNames(aggregate(trimm_read_table$count,list(trimm_read_table$sample),
                               function(x) sum(x)),c("sample","readnumber"))
min_trimm<-setNames(aggregate(trimm_read_table$length,list(trimm_read_table$sample),
                            function(x) min(x)),c("sample","minlength"))
max_trimm<-setNames(aggregate(trimm_read_table$length,list(trimm_read_table$sample),
                            function(x) max(x)),c("sample","maxlength"))
mean_trimm<-setNames(aggregate(trimm_read_table$length,list(trimm_read_table$sample),
                             function(x) mean(x)),c("sample","meanlength"))
median_trimm<-setNames(aggregate(trimm_read_table$length,list(trimm_read_table$sample),
                               function(x) median(x)),c("sample","medianlength"))

trimm_metric<-merge(nbread_trimm,min_trimm,by="sample")
trimm_metric<-merge(trimm_metric,max_trimm,by="sample")
trimm_metric<-merge(trimm_metric,mean_trimm,by="sample")
trimm_metric<-merge(trimm_metric,median_trimm,by="sample")
trimm_metric$status<-"TRIMM"

##############################################

dehost_read_table<-read.table(dehost_file,header=F)
colnames(dehost_read_table)<-c("sequence","sample")
dehost_read_table$length<-nchar(as.character(dehost_read_table$sequence))
dehost_read_table$count<-1
dehost_read_table<-dehost_read_table[,-1]

nbread_dehost<-setNames(aggregate(dehost_read_table$count,list(dehost_read_table$sample),
                                  function(x) sum(x)),c("sample","readnumber"))
min_dehost<-setNames(aggregate(dehost_read_table$length,list(dehost_read_table$sample),
                               function(x) min(x)),c("sample","minlength"))
max_dehost<-setNames(aggregate(dehost_read_table$length,list(dehost_read_table$sample),
                               function(x) max(x)),c("sample","maxlength"))
mean_dehost<-setNames(aggregate(dehost_read_table$length,list(dehost_read_table$sample),
                                function(x) mean(x)),c("sample","meanlength"))
median_dehost<-setNames(aggregate(dehost_read_table$length,list(dehost_read_table$sample),
                                  function(x) median(x)),c("sample","medianlength"))

dehost_metric<-merge(nbread_dehost,min_dehost,by="sample")
dehost_metric<-merge(dehost_metric,max_dehost,by="sample")
dehost_metric<-merge(dehost_metric,mean_dehost,by="sample")
dehost_metric<-merge(dehost_metric,median_dehost,by="sample")
dehost_metric$status<-"DEHOST"

##############################################

BM_read_table<-read.table(BM_file,header=F)
colnames(BM_read_table)<-c("sequence","sample")
BM_read_table$length<-nchar(as.character(BM_read_table$sequence))
BM_read_table$count<-1
BM_read_table<-BM_read_table[,-1]

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
BM_metric$status<-"BESTMATCH_FILTERED"

summ_metric<-rbind(raw_metric,trimm_metric,dehost_metric,BM_metric)


write.csv2(summ_metric,output_summ_table,row.names = F)
