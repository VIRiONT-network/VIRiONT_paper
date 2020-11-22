argv <- commandArgs(TRUE)

rawtablepath<-as.character(argv[1])
dehosttablepath<-as.character(argv[2])
trimtablepath<-as.character(argv[3])
reftablepath<-as.character(argv[4])
outputtablepath<-as.character(argv[5])

raw_count<-read.csv2(rawtablepath,header=F)
dehost_count<-read.csv2(dehosttablepath,header=F)
trimm_count<-read.csv2(trimtablepath,header=F)
ref_count<-read.csv2(reftablepath,header=F)

ALLDATA<-rbind(raw_count,dehost_count,trimm_count,ref_count)
colnames(ALLDATA)<-c("readlength","sample","step","assignedref")
ALLDATA$count<-1

countread_data<-aggregate(ALLDATA$count,by=list(ALLDATA$sample,ALLDATA$step,ALLDATA$assignedref),sum)
colnames(countread_data)<-c("sample","step","assignedref","read_count")

minlengthread_data<-aggregate(ALLDATA$readlength,by=list(ALLDATA$sample,ALLDATA$step,ALLDATA$assignedref),min)
colnames(minlengthread_data)<-c("sample","step","assignedref","minlengthread")

maxlengthread_data<-aggregate(ALLDATA$readlength,by=list(ALLDATA$sample,ALLDATA$step,ALLDATA$assignedref),max)
colnames(maxlengthread_data)<-c("sample","step","assignedref","maxlengthread")

meanread_data<-aggregate(ALLDATA$readlength,by=list(ALLDATA$sample,ALLDATA$step,ALLDATA$assignedref),mean)
colnames(meanread_data)<-c("sample","step","assignedref","meanread_length")

medianread_data<-aggregate(ALLDATA$readlength,by=list(ALLDATA$sample,ALLDATA$step,ALLDATA$assignedref),median)
colnames(medianread_data)<-c("sample","step","assignedref","medianread_length")

METRIC_data<-merge(countread_data,minlengthread_data,by=c("sample","step","assignedref"))
METRIC_data<-merge(METRIC_data,maxlengthread_data,by=c("sample","step","assignedref"))
METRIC_data<-merge(METRIC_data,meanread_data,by=c("sample","step","assignedref"))
METRIC_data<-merge(METRIC_data,medianread_data,by=c("sample","step","assignedref"))

METRIC_data<-METRIC_data[order(METRIC_data$sample, METRIC_data$step),]

write.csv2(METRIC_data,outputtablepath,row.names = F,quote=F)
