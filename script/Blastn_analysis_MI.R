library(ggplot2)

argv <- commandArgs(TRUE)

blastn_data<-read.table(argv[1],header=F,sep="\t")
#blastn_data<-read.table("C:/Users/regueex/Desktop/test_ViralION/barcode05_fmt.txt")

analysis_table<-read.csv2(argv[2])
#analysis_table<-read.csv2("C:/Users/regueex/Desktop/test_ViralION/table_analysis.csv")

analysis<-as.character(argv[3])
#analysis<-"HBV_REF"

output_blastnresult<-as.character(argv[4])

output_plot<-as.character(argv[5])

cutoff<-as.numeric(argv[6])
#cutoff<-10

samplename<-as.character(argv[7])
#samplename<-"BC01"

multiinf_table_name<-as.character(argv[8])

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

write.table(TABLE_BEST_GENO,output_blastnresult,sep="\t",row.names=F,quote=F)

count_Geno<-c()
label<-c()
for (ref in reference_list) {
  subtable<-subset(TABLE_BEST_GENO,GENOTYPE==ref)
  label<-c(label,ref)
  count_Geno<-c(count_Geno,nrow(subtable))
}

hist_data<-as.data.frame(count_Geno)
hist_data$Genotype<-label

#########    DETECT  MULTI INFECTIONS    ###################

hist_data$Ratio_Bestref<-hist_data$count_Geno/max(hist_data$count_Geno)*100

#filtrate result depending on cutoff
multiinf_table<-subset(hist_data,Ratio_Bestref>=cutoff )
multiinf_table$sample<-samplename
multiinf_table<-multiinf_table[,c("sample","Genotype","Ratio_Bestref")]

write.table(multiinf_table,multiinf_table_name,sep="\t",row.names = F,col.names = F,quote=F)

############################################################
hist_data$Ratio_Bestref<-trunc(hist_data$Ratio_Bestref)
hist_data<-subset(hist_data,count_Geno>0) #Subset for plot with plenty ref

png(filename = output_plot)
ggplot(data=hist_data, aes(x=Genotype, y=Ratio_Bestref)) +
  geom_bar(stat="identity",color="black",fill="steelblue")+
  geom_text(aes(label=Ratio_Bestref), vjust=-1, color="black", size=4)+
  ggtitle("Reference repartition per read")+
  labs(y= "percentage reference/best_reference", x = "reference")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
