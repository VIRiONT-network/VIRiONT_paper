argv <- commandArgs(TRUE)

read_list<-read.table(argv[1],header=F)
sample<-argv[2]
ref<-argv[3]
sum_file<-argv[4]

read_list$V1<-as.character(read_list$V1)
read_list$LEN<-nchar(read_list$V1)

nb_read<-nrow(read_list)
mean_len<-mean(read_list$LEN)
median_len<-median(read_list$LEN)

vec_result<-c(sample,ref,nb_read,mean_len,median_len)

writeLines(vec_result,sum_file, sep = "\t", useBytes = FALSE)
