library(ggplot2)

argv <- commandArgs(TRUE)

blastn_data<-read.table(argv[1],header=F,sep="\t")
#blastn_data<-read.table("C:/Users/regueex/Desktop/test_ViralION/barcode05_fmt.txt")

analysis_table<-read.csv2(argv[2])
#analysis_table<-read.csv2("C:/Users/regueex/Desktop/test_ViralION/table_analysis.csv")

analysis<-as.character(argv[3])
#analysis<-"HBV_REF"

output_readlist<-as.character(argv[4])
output_bestref<-as.character(argv[5])
output_plot<-as.character(argv[6])

reference_list<-as.character(analysis_table[,which(colnames(analysis_table)==analysis)])
reference_list<-subset(reference_list,reference_list!="" & !(is.na(reference_list)))


table_f<-blastn_data[,-c(3:11)]
colnames(table_f)<-c("READ","GENOTYPE","SCORE")

list_read<-as.character(table_f$READ)
list_read_uniq<-unique(list_read)

MAX_SCORE<-aggregate(table_f$SCORE, list(table_f$READ), function(x) max(x))
colnames(MAX_SCORE)<-c("READ","MAX")

SUMMARY_TABLE<-merge(table_f,MAX_SCORE,by = "READ")

TABLE_BEST_GENO<-subset(SUMMARY_TABLE,SCORE==MAX)

count_Geno<-c()
label<-c()
for (ref in reference_list) {
  subtable<-subset(TABLE_BEST_GENO,GENOTYPE==ref)
  label<-c(label,ref)
  count_Geno<-c(count_Geno,nrow(subtable))
}

hist_data<-as.data.frame(count_Geno)
hist_data$Genotype<-label

png(filename = output_plot)
ggplot(data=hist_data, aes(x=Genotype, y=count_Geno)) +
  geom_bar(stat="identity",color="black",fill="steelblue")+
  geom_text(aes(label=count_Geno), vjust=1.6, color="white", size=3.5)+
  ggtitle("Reference repartition per read")+
  labs(y= "read count", x = "genotype")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

hist_data<-hist_data[order(count_Geno,decreasing = T),]

major_geno<-as.character(hist_data[1,2])
subtable<-subset(TABLE_BEST_GENO,GENOTYPE==major_geno)
list_fastq<-as.character(subtable$READ)
write(list_fastq,output_readlist)
write(major_geno,output_bestref)
