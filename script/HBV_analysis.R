library(ggplot2)

argv <- commandArgs(TRUE)

table<-read.table(argv[1],header=F,sep="\t")
path=as.character(argv[2])

string<-as.character(argv[1])
split1<-strsplit(string,"/")
filename<-split1[[1]][length(split1[[1]])]
split2<-strsplit(filename,"_")
barcode<-split2[[1]][1]

list_read<-as.character(table$V1)
list_read_uniq<-unique(list_read)

size_table<-as.numeric(length(list_read_uniq))
nb_value<-size_table*11# 10 genotypes+read column

table_result<-matrix(1:nb_value,nrow = size_table ,ncol = 11) #
table_result<-as.data.frame(table_result)

colnames(table_result)<-c("READ","P3_GTA","P3_GTB","P3_GTC","P3_GTD","P3_GTE","P3_GTF","P3_GTG","P3_GTH","P3_GTI","P3_GTJ")
table_result[,1]<-as.character(list_read_uniq)
for (i in 2:11) table_result[,i]<-as.numeric(table_result[,i])
table_result[,2:11]<-0

getBestGeno<-function(readName){
  sub_table<-table[table$V1==readName,]
  score_max<-max(sub_table$V12)
  best_P3_GT<-sub_table[sub_table$V12==score_max,]
  best_geno_list<-as.character(best_P3_GT$V2)
  return(best_geno_list)
}

for (i in 1:nrow(table_result)){ #nrow(table_result)
  read<-as.character(table_result[i,1])
  best_geno_list<-getBestGeno(read)
  if (("P3_GTA" %in% best_geno_list)==TRUE) table_result[i,2]<-1
  if (("P3_GTB" %in% best_geno_list)==TRUE) table_result[i,3]<-1
  if (("P3_GTC" %in% best_geno_list)==TRUE) table_result[i,4]<-1
  if (("P3_GTD" %in% best_geno_list)==TRUE) table_result[i,5]<-1
  if (("P3_GTE" %in% best_geno_list)==TRUE) table_result[i,6]<-1
  if (("P3_GTF" %in% best_geno_list)==TRUE) table_result[i,7]<-1
  if (("P3_GTG" %in% best_geno_list)==TRUE) table_result[i,8]<-1
  if (("P3_GTH" %in% best_geno_list)==TRUE) table_result[i,9]<-1
  if (("P3_GTI" %in% best_geno_list)==TRUE) table_result[i,10]<-1
  if (("P3_GTJ" %in% best_geno_list)==TRUE) table_result[i,11]<-1
}

label<-c("P3_GTA","P3_GTB","P3_GTC","P3_GTD","P3_GTE","P3_GTF","P3_GTG","P3_GTH","P3_GTI","P3_GTJ")
count_Geno<-c()
count_Geno[1]<-sum(as.numeric(table_result[,2]))
count_Geno[2]<-sum(as.numeric(table_result[,3]))
count_Geno[3]<-sum(as.numeric(table_result[,4]))
count_Geno[4]<-sum(as.numeric(table_result[,5]))
count_Geno[5]<-sum(as.numeric(table_result[,6]))
count_Geno[6]<-sum(as.numeric(table_result[,7]))
count_Geno[7]<-sum(as.numeric(table_result[,8]))
count_Geno[8]<-sum(as.numeric(table_result[,9]))
count_Geno[9]<-sum(as.numeric(table_result[,10]))
count_Geno[10]<-sum(as.numeric(table_result[,11]))

hist_data<-as.data.frame(count_Geno)
hist_data$Genotype<-label

png(filename = paste0(path,"RDATA/",barcode,"_barplot.png"))
ggplot(data=hist_data, aes(x=Genotype, y=count_Geno)) +
  geom_bar(stat="identity",color="black",fill="steelblue")+
  geom_text(aes(label=count_Geno), vjust=1.6, color="white", size=3.5)+
  ggtitle("Genotype repartition per read")+
  labs(y= "read count", x = "genotype")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

hist_data<-hist_data[order(count_Geno,decreasing = T),]

major_geno<-as.character(hist_data[1,2])


if (major_geno=="P3_GTA") {
  read_P3_GTA_majo<-table_result[table_result$P3_GTA==1,]
  list_fastq_GTA<-as.character(read_P3_GTA_majo$READ)
  write(list_fastq_GTA,paste0(path,"R_RESULT/",barcode,"_list.txt"))
  write(major_geno,paste0(path,"R_RESULT/",barcode,"_bestgeno.txt"))
}

if (major_geno=="P3_GTB") {
  read_P3_GTB_majo<-table_result[table_result$P3_GTB==1,]
  list_fastq_GTB<-as.character(read_P3_GTB_majo$READ)
  write(list_fastq_GTB,paste0(path,"R_RESULT/",barcode,"_list.txt"))
  write(major_geno,paste0(path,"R_RESULT/",barcode,"_bestgeno.txt"))
}

if (major_geno=="P3_GTC") {
  read_P3_GTC_majo<-table_result[table_result$P3_GTC==1,]
  list_fastq_GTC<-as.character(read_P3_GTC_majo$READ)
  write(list_fastq_GTC,paste0(path,"R_RESULT/",barcode,"_list.txt"))
  write(major_geno,paste0(path,"R_RESULT/",barcode,"_bestgeno.txt"))
}

if (major_geno=="P3_GTD") {
  read_P3_GTD_majo<-table_result[table_result$P3_GTD==1,]
  list_fastq_GTD<-as.character(read_P3_GTD_majo$READ)
  write(list_fastq_GTD,paste0(path,"R_RESULT/",barcode,"_list.txt"))
  write(major_geno,paste0(path,"R_RESULT/",barcode,"_bestgeno.txt"))
}

if (major_geno=="P3_GTE") {
  read_P3_GTE_majo<-table_result[table_result$P3_GTE==1,]
  list_fastq_GTE<-as.character(read_P3_GTE_majo$READ)
  write(list_fastq_GTE,paste0(path,"R_RESULT/",barcode,"_list.txt"))
  write(major_geno,paste0(path,"R_RESULT/",barcode,"_bestgeno.txt"))
}

if (major_geno=="P3_GTF") {
  read_P3_GTF_majo<-table_result[table_result$P3_GTF==1,]
  list_fastq_GTF<-as.character(read_P3_GTF_majo$READ)
  write(list_fastq_GTF,paste0(path,"R_RESULT/",barcode,"_list.txt"))
  write(major_geno,paste0(path,"R_RESULT/",barcode,"_bestgeno.txt"))
}

if (major_geno=="P3_GTG") {
  read_P3_GTG_majo<-table_result[table_result$P3_GTG==1,]
  list_fastq_GTG<-as.character(read_P3_GTG_majo$READ)
  write(list_fastq_GTG,paste0(path,"R_RESULT/",barcode,"_list.txt"))
  write(major_geno,paste0(path,"R_RESULT/",barcode,"_bestgeno.txt"))
}

if (major_geno=="P3_GTH") {
  read_P3_GTH_majo<-table_result[table_result$P3_GTH==1,]
  list_fastq_GTH<-as.character(read_P3_GTH_majo$READ)
  write(list_fastq_GTH,paste0(path,"R_RESULT/",barcode,"_list.txt"))
  write(major_geno,paste0(path,"R_RESULT/",barcode,"_bestgeno.txt"))
}

if (major_geno=="P3_GTI") {
  read_P3_GTI_majo<-table_result[table_result$P3_GTI==1,]
  list_fastq_GTI<-as.character(read_P3_GTI_majo$READ)
  write(list_fastq_GTI,paste0(path,"R_RESULT/",barcode,"_list.txt"))
  write(major_geno,paste0(path,"R_RESULT/",barcode,"_bestgeno.txt"))
}

if (major_geno=="P3_GTJ") {
  read_P3_GTJ_majo<-table_result[table_result$P3_GTJ==1,]
  list_fastq_GTJ<-as.character(read_P3_GTJ_majo$READ)
  write(list_fastq_GTJ,paste0(path,"R_RESULT/",barcode,"_list.txt"))
  write(major_geno,paste0(path,"R_RESULT/",barcode,"_bestgeno.txt"))
}

save.image(paste0(path,"RDATA/",barcode,"_rsave.RData"))