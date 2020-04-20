library(ggplot2)

argv <- commandArgs(TRUE)

table<-read.table(argv[1],header=F,sep="\t")
output=as.character(argv[2])

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

TABLE_A1<-subset(TABLE_BEST_GENO,GENOTYPE=="A1")
TABLE_A2<-subset(TABLE_BEST_GENO,GENOTYPE=="A2")
TABLE_A3<-subset(TABLE_BEST_GENO,GENOTYPE=="A3")
TABLE_A4<-subset(TABLE_BEST_GENO,GENOTYPE=="A4")
TABLE_A5<-subset(TABLE_BEST_GENO,GENOTYPE=="A5")
TABLE_A6<-subset(TABLE_BEST_GENO,GENOTYPE=="A6")
TABLE_B1<-subset(TABLE_BEST_GENO,GENOTYPE=="B1")
TABLE_B2<-subset(TABLE_BEST_GENO,GENOTYPE=="B2")
TABLE_B3<-subset(TABLE_BEST_GENO,GENOTYPE=="B3")
TABLE_B4<-subset(TABLE_BEST_GENO,GENOTYPE=="B4")
TABLE_B6<-subset(TABLE_BEST_GENO,GENOTYPE=="B6")
TABLE_C1<-subset(TABLE_BEST_GENO,GENOTYPE=="C1")
TABLE_C2<-subset(TABLE_BEST_GENO,GENOTYPE=="C2")
TABLE_C4<-subset(TABLE_BEST_GENO,GENOTYPE=="C4")
TABLE_C5<-subset(TABLE_BEST_GENO,GENOTYPE=="C5")
TABLE_C6<-subset(TABLE_BEST_GENO,GENOTYPE=="C6")
TABLE_C8<-subset(TABLE_BEST_GENO,GENOTYPE=="C8")
TABLE_C10<-subset(TABLE_BEST_GENO,GENOTYPE=="C10")
TABLE_C11<-subset(TABLE_BEST_GENO,GENOTYPE=="C11")
TABLE_D1<-subset(TABLE_BEST_GENO,GENOTYPE=="D1")
TABLE_D2<-subset(TABLE_BEST_GENO,GENOTYPE=="D2")
TABLE_D3<-subset(TABLE_BEST_GENO,GENOTYPE=="D3")
TABLE_D4<-subset(TABLE_BEST_GENO,GENOTYPE=="D4")
TABLE_D5<-subset(TABLE_BEST_GENO,GENOTYPE=="D5")
TABLE_D7<-subset(TABLE_BEST_GENO,GENOTYPE=="D7")
TABLE_E<-subset(TABLE_BEST_GENO,GENOTYPE=="E")
TABLE_F1<-subset(TABLE_BEST_GENO,GENOTYPE=="F1")
TABLE_F2<-subset(TABLE_BEST_GENO,GENOTYPE=="F2")
TABLE_F3<-subset(TABLE_BEST_GENO,GENOTYPE=="F3")
TABLE_F4<-subset(TABLE_BEST_GENO,GENOTYPE=="F4")
TABLE_G<-subset(TABLE_BEST_GENO,GENOTYPE=="G")
TABLE_H<-subset(TABLE_BEST_GENO,GENOTYPE=="H")
TABLE_I1<-subset(TABLE_BEST_GENO,GENOTYPE=="I1")
TABLE_I2<-subset(TABLE_BEST_GENO,GENOTYPE=="I2")
TABLE_J<-subset(TABLE_BEST_GENO,GENOTYPE=="J")

label<-c("A1","A2","A3","A4","A5","A6","B1","B2","B3","B4","B6","C1","C2","C4","C5","C6","C8","C10","C11","D1","D2","D3","D4","D5","D7","E","F1","F2","F3","F4","G","H","I1","I2","J")
count_Geno<-c()
count_Geno[1]<-length(TABLE_A1$READ)
count_Geno[2]<-length(TABLE_A2$READ)
count_Geno[3]<-length(TABLE_A3$READ)
count_Geno[4]<-length(TABLE_A4$READ)
count_Geno[5]<-length(TABLE_A5$READ)
count_Geno[6]<-length(TABLE_A6$READ)
count_Geno[7]<-length(TABLE_B1$READ)
count_Geno[8]<-length(TABLE_B2$READ)
count_Geno[9]<-length(TABLE_B3$READ)
count_Geno[10]<-length(TABLE_B4$READ)
count_Geno[11]<-length(TABLE_B6$READ)
count_Geno[12]<-length(TABLE_C1$READ)
count_Geno[13]<-length(TABLE_C2$READ)
count_Geno[14]<-length(TABLE_C4$READ)
count_Geno[15]<-length(TABLE_C5$READ)
count_Geno[16]<-length(TABLE_C6$READ)
count_Geno[17]<-length(TABLE_C8$READ)
count_Geno[18]<-length(TABLE_C10$READ)
count_Geno[19]<-length(TABLE_C11$READ)
count_Geno[20]<-length(TABLE_D1$READ)
count_Geno[21]<-length(TABLE_D2$READ)
count_Geno[22]<-length(TABLE_D3$READ)
count_Geno[23]<-length(TABLE_D4$READ)
count_Geno[24]<-length(TABLE_D5$READ)
count_Geno[25]<-length(TABLE_D7$READ)
count_Geno[26]<-length(TABLE_E$READ)
count_Geno[27]<-length(TABLE_F1$READ)
count_Geno[28]<-length(TABLE_F2$READ)
count_Geno[29]<-length(TABLE_F3$READ)
count_Geno[30]<-length(TABLE_F4$READ)
count_Geno[31]<-length(TABLE_G$READ)
count_Geno[32]<-length(TABLE_H$READ)
count_Geno[33]<-length(TABLE_I1$READ)
count_Geno[34]<-length(TABLE_I2$READ)
count_Geno[35]<-length(TABLE_J$READ)

hist_data<-as.data.frame(count_Geno)
hist_data$Genotype<-label

#PLOT
png(filename = output)
ggplot(data=hist_data, aes(x=Genotype, y=count_Geno)) +
  geom_bar(stat="identity",color="black",fill="steelblue")+
  geom_text(aes(label=count_Geno), vjust=1.6, color="black", size=3.5)+
  ggtitle("Genotype repartition per read longer than 1500 nt")+
  labs(y= "Reads count", x = "Genotypes")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()






