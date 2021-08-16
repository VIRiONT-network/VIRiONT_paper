library(ggplot2)
library(gridExtra)

argv <- commandArgs(TRUE)

ref_count<-read.table(argv[1],header=F,sep="\t",col.names =c("reference","count","sample"))
MI_cutoff<-as.numeric(argv[2])
pdf_output_file<-as.character(argv[3])
output_table_MI<-as.character(argv[4])

list_sample<-sort(unique(as.character(ref_count$sample)))

#plot results
pdf(pdf_output_file,paper = "a4r")
for (samplename in list_sample){
  subset_refcount<-subset(ref_count,sample==samplename)
  subset_refcount$ratio_MI<-subset_refcount$count/max(subset_refcount$count)*100
  subset_refcount$ratio_MI<-format(subset_refcount$ratio_MI,scientific = FALSE)
  subset_refcount$ratio_MI<-round(as.numeric(subset_refcount$ratio_MI),1)

  raw_plot<-ggplot(data=subset_refcount, aes(x=reference, y=count)) +
    geom_bar(stat="identity",color="black",fill="steelblue")+
    geom_text(aes(label=count), vjust=-0.5, color="black", size=2)+
    ggtitle(paste0("Reference repartition per read for sample:",samplename))+
    labs(y= "raw read count", x = "reference")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))

  ratio_plot<-ggplot(data=subset_refcount, aes(x=reference, y=ratio_MI)) +
    geom_bar(stat="identity",color="black",fill="steelblue")+
    geom_hline(yintercept = MI_cutoff,linetype="dashed",color="red",size=1)+
    geom_text(aes(label=ratio_MI), vjust=-0.5, color="black", size=2)+
    ggtitle(paste0("matching ratio for sample:",samplename))+
    labs(y= "percentage reference/best_reference", x = "reference")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))
  
  grid.arrange(raw_plot,ratio_plot, nrow=2,ncol=1)
}
dev.off()

#get MI results for VIRiONT
max_refcount<-aggregate(ref_count$count,by=list(ref_count$sample),function(x) max(x))
colnames(max_refcount)<-c("sample","maxcount")

all_data<-merge(ref_count,max_refcount,by="sample")
all_data$MI_ratio<-all_data$count/all_data$maxcount*100
all_data$MI_ratio<-format(all_data$MI_ratio,scientific = FALSE)
all_data$MI_ratio<-round(as.numeric(all_data$MI_ratio),1)

#subset data
matchingMI_table<-subset(all_data,MI_ratio>=MI_cutoff)
matchingMI_table<-matchingMI_table[,c("sample","reference","count","MI_ratio")]

#minimum read number to be processed in the 2nd workflow
nbread_cutoff<-50

matchingMI_table<-subset(matchingMI_table,count>=nbread_cutoff)

write.table(matchingMI_table,output_table_MI,row.names = F,quote=F,col.names = F,sep="\t")
