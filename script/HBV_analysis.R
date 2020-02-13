library(ggplot2)

argv <- commandArgs(TRUE)

table<-read.table(argv[1],header=F,sep="\t")
path=as.character(argv[2])

string<-as.character(argv[1])
split1<-strsplit(string,"/")
filename<-split1[[1]][length(split1[[1]])]
split2<-strsplit(filename,"_")
barcode<-split2[[1]][1]

table_f<-table[,-c(3:11)]
colnames(table_f)<-c("READ","GENOTYPE","SCORE")

list_read<-as.character(table_f$READ)
list_read_uniq<-unique(list_read)


MAX_SCORE<-aggregate(table_f$SCORE, list(table_f$READ), function(x) max(x))
colnames(MAX_SCORE)<-c("READ","MAX")

SUMMARY_TABLE<-merge(table_f,MAX_SCORE,by = "READ")

TABLE_BEST_GENO<-subset(SUMMARY_TABLE,SCORE==MAX)

TABLE_GTA<-subset(TABLE_BEST_GENO,GENOTYPE=="P3_GTA")
TABLE_GTB<-subset(TABLE_BEST_GENO,GENOTYPE=="P3_GTB")
TABLE_GTC<-subset(TABLE_BEST_GENO,GENOTYPE=="P3_GTC")
TABLE_GTD<-subset(TABLE_BEST_GENO,GENOTYPE=="P3_GTD")
TABLE_GTE<-subset(TABLE_BEST_GENO,GENOTYPE=="P3_GTE")
TABLE_GTF<-subset(TABLE_BEST_GENO,GENOTYPE=="P3_GTF")
TABLE_GTG<-subset(TABLE_BEST_GENO,GENOTYPE=="P3_GTG")
TABLE_GTH<-subset(TABLE_BEST_GENO,GENOTYPE=="P3_GTH")
TABLE_GTI<-subset(TABLE_BEST_GENO,GENOTYPE=="P3_GTI")
TABLE_GTJ<-subset(TABLE_BEST_GENO,GENOTYPE=="P3_GTJ")

label<-c("P3_GTA","P3_GTB","P3_GTC","P3_GTD","P3_GTE","P3_GTF","P3_GTG","P3_GTH","P3_GTI","P3_GTJ")
count_Geno<-c()
count_Geno[1]<-length(TABLE_GTA$READ)
count_Geno[2]<-length(TABLE_GTB$READ)
count_Geno[3]<-length(TABLE_GTC$READ)
count_Geno[4]<-length(TABLE_GTD$READ)
count_Geno[5]<-length(TABLE_GTE$READ)
count_Geno[6]<-length(TABLE_GTF$READ)
count_Geno[7]<-length(TABLE_GTG$READ)
count_Geno[8]<-length(TABLE_GTH$READ)
count_Geno[9]<-length(TABLE_GTI$READ)
count_Geno[10]<-length(TABLE_GTJ$READ)


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
  list_fastq_GTA<-as.character(TABLE_GTA$READ)
  write(list_fastq_GTA,paste0(path,"R_RESULT/",barcode,"_list.txt"))
  write(major_geno,paste0(path,"R_RESULT/",barcode,"_bestgeno.txt"))
}
if (major_geno=="P3_GTB") {
  list_fastq_GTB<-as.character(TABLE_GTB$READ)
  write(list_fastq_GTB,paste0(path,"R_RESULT/",barcode,"_list.txt"))
  write(major_geno,paste0(path,"R_RESULT/",barcode,"_bestgeno.txt"))
}
if (major_geno=="P3_GTC") {
  list_fastq_GTC<-as.character(TABLE_GTC$READ)
  write(list_fastq_GTC,paste0(path,"R_RESULT/",barcode,"_list.txt"))
  write(major_geno,paste0(path,"R_RESULT/",barcode,"_bestgeno.txt"))
}
if (major_geno=="P3_GTD") {
  list_fastq_GTD<-as.character(TABLE_GTD$READ)
  write(list_fastq_GTD,paste0(path,"R_RESULT/",barcode,"_list.txt"))
  write(major_geno,paste0(path,"R_RESULT/",barcode,"_bestgeno.txt"))
}
if (major_geno=="P3_GTE") {
  list_fastq_GTE<-as.character(TABLE_GTE$READ)
  write(list_fastq_GTE,paste0(path,"R_RESULT/",barcode,"_list.txt"))
  write(major_geno,paste0(path,"R_RESULT/",barcode,"_bestgeno.txt"))
}
if (major_geno=="P3_GTF") {
  list_fastq_GTF<-as.character(TABLE_GTF$READ)
  write(list_fastq_GTF,paste0(path,"R_RESULT/",barcode,"_list.txt"))
  write(major_geno,paste0(path,"R_RESULT/",barcode,"_bestgeno.txt"))
}
if (major_geno=="P3_GTG") {
  list_fastq_GTG<-as.character(TABLE_GTG$READ)
  write(list_fastq_GTG,paste0(path,"R_RESULT/",barcode,"_list.txt"))
  write(major_geno,paste0(path,"R_RESULT/",barcode,"_bestgeno.txt"))
}
if (major_geno=="P3_GTH") {
  list_fastq_GTH<-as.character(TABLE_GTH$READ)
  write(list_fastq_GTH,paste0(path,"R_RESULT/",barcode,"_list.txt"))
  write(major_geno,paste0(path,"R_RESULT/",barcode,"_bestgeno.txt"))
}
if (major_geno=="P3_GTI") {
  list_fastq_GTI<-as.character(TABLE_GTI$READ)
  write(list_fastq_GTI,paste0(path,"R_RESULT/",barcode,"_list.txt"))
  write(major_geno,paste0(path,"R_RESULT/",barcode,"_bestgeno.txt"))
}
if (major_geno=="P3_GTJ") {
  list_fastq_GTJ<-as.character(TABLE_GTI$READ)
  write(list_fastq_GTJ,paste0(path,"R_RESULT/",barcode,"_list.txt"))
  write(major_geno,paste0(path,"R_RESULT/",barcode,"_bestgeno.txt"))
}

#write(list_fastq_GTA,paste0(path,"RDATA/",barcode,"_GTA_list.txt"))
#write(list_fastq_GTB,paste0(path,"RDATA/",barcode,"_GTB_list.txt"))
#write(list_fastq_GTC,paste0(path,"RDATA/",barcode,"_GTC_list.txt"))
#write(list_fastq_GTD,paste0(path,"RDATA/",barcode,"_GTD_list.txt"))
#write(list_fastq_GTE,paste0(path,"RDATA/",barcode,"_GTE_list.txt"))
#write(list_fastq_GTF,paste0(path,"RDATA/",barcode,"_GTF_list.txt"))
#write(list_fastq_GTG,paste0(path,"RDATA/",barcode,"_GTG_list.txt"))
#write(list_fastq_GTH,paste0(path,"RDATA/",barcode,"_GTH_list.txt"))
#write(list_fastq_GTI,paste0(path,"RDATA/",barcode,"_GTI_list.txt"))
#write(list_fastq_GTJ,paste0(path,"RDATA/",barcode,"_GTJ_list.txt"))

save.image(paste0(path,"RDATA/",barcode,"_rsave.RData"))
