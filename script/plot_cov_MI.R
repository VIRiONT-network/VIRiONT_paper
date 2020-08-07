library(ggplot2)

argv <- commandArgs(TRUE)

covdata<-read.table(argv[1],header=F)
colnames(covdata)<-c("REF","POS","COV","SAMPLE","REFERENCE")

covplot<-ggplot(covdata, aes(x=POS, y=COV),fill=REFERENCE)+
  geom_line(aes(color=REFERENCE),show.legend = FALSE,alpha=.5) +
  ylab("Coverage") + xlab("Genomic position")+
  theme(legend.position="none")+
  facet_wrap(vars(SAMPLE,REFERENCE),scales = "free",ncol = 2)


ggsave(filename=argv[2],plot=covplot,path=NULL,width=10,height=2*(length(unique(covdata$SAMPLE))),limitsize = FALSE)
