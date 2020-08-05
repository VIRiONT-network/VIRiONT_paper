library(ggplot2)

argv <- commandArgs(TRUE)

covdata<-read.table(argv[1],header=F)
colnames(covdata)<-c("REF","POS","COV","SAMPLE")

covplot<-ggplot(covdata, aes(x=POS, y=COV,fill=SAMPLE))+
  geom_bar(stat="identity",width=1)+
  ylab("Coverage") + xlab("Genomic position")+
  theme_bw()+ theme(legend.position="none")+
  facet_wrap(.~SAMPLE, scales="free",ncol=2)

ggsave(filename=argv[2],plot=covplot,path=NULL,width=10,height=2*(length(unique(covdata$SAMPLE))),limitsize = FALSE)
