
setwd("C:/Ubuntu/longmatchquery")
library("Rcpp")
sourceCpp(file="hashmap.cpp")
setwd("C:/Ubuntu/longmatchquery/bigeg2")
#map=read.table('p.map')[,3:4]
#gd=map[,2]/100000000
refindex=rep(c(0,1,2,3),each=2000)

starttime=proc.time()
paintings=paintingalldense(refindex,nind=500,targetfrac=1.0,fixrho=FALSE,method="Viterbi",
                          minmatchfrac = 0.0005,L_initial=320,L_minmatch=20)
proc.time()-starttime

cc=as.data.frame(paintings[[1]])
image(as.matrix(cc[,-1]))

cc2=as.data.frame(paintings[[1]])
image(as.matrix(cc2[,-1]))

diff=abs(cc-cc2)[,-1]
plot(y=rowMeans(diff),x=cc[,1])

painting=paintingalldense(refindex,500,fixrho=TRUE,indfrac=0.001,minmatchfrac = 0.0001)
painting=paintingalldense(refindex,500,0.1,ite_time=10,method="EM",
                          indfrac=0.1,EMsnpfrac=0.3,minmatchfrac = 0.001,L_initial=320,L_minmatch=20)

coancestry=chunklengthall(refindex,minmatchfrac = 0.001)


aa=as.data.frame(paintings[[22]])


allscore=vector()
for(i in 1:1000){
  print(i)
  allscore<-c(allscore,as.numeric(t(as.data.frame(paintings[[i]]))[1,]))
}



## comparison of painting results
setwd("C:/Ubuntu/longmatchquery")
library("Rcpp")
sourceCpp(file="hashmap.cpp")
setwd("C:/Ubuntu/longmatchquery/bigeg2")
map=read.table('p.map')[,3:4]
gd=map[,2]/100000000
refindex=rep(c(0,1,2,3),each=2000)

compare <- function(minmatchfrac,L_minmatch){
  setwd("C:/Ubuntu/longmatchquery/bigeg2")
  starttime=proc.time()
  paintings=paintingalldense(refindex,nind=50,targetfrac=1.0,fixrho=FALSE,method="Viterbi",
                             minmatchfrac = minmatchfrac,L_initial=160,L_minmatch=L_minmatch)
  time=(proc.time()-starttime)[3]
  diff_all<-vector()
  reliability <- vector()
  matchfrac <- vector()
  entropy_all <- vector()
  entropy_match <- vector()
  entropy_unmatch <- vector()
  setwd("C:/Ubuntu/longmatchquery/bigeg2/chromopainter")
  for(i in 1:100){
    print(i)
    p_chromo=read.table(paste0("test_",i,".txt"))[,-1]
    p_pbwt=as.data.frame(paintings[[i]])[,-1]
    #reliability=c(reliability,p_pbwt[,1])
    diff=unlist(abs(p_chromo-p_pbwt))
    diff_all=c(diff_all,diff)
    maxinfo1=apply(p_chromo,1,which.max)
    maxinfo2=apply(p_pbwt,1,which.max)
    matchfrac <- c(matchfrac,sum(diag(table(data.frame(maxinfo1,maxinfo2))))/2091)
    entropy=apply(p_pbwt+1e-20,1,function(x)-sum(x*log(x)))
    entropy_all=c(entropy_all,entropy)
    entropy_unmatch <- c(entropy_unmatch,mean(entropy[which((maxinfo1==maxinfo2)==FALSE)]))
    entropy_match <- c(entropy_match,mean(entropy[which((maxinfo1==maxinfo2)==TRUE)]))
  }
  return(list(diff_all,matchfrac,time,entropy_match,entropy_unmatch,entropy_all))
}




#absdiff_0.001matches_Lmin10=compare(0.001,10)
#absdiff_0.001matches_Lmin20=compare(0.001,20)
#absdiff_0.001matches_Lmin40=compare(0.001,40)
#absdiff_0.002matches_Lmin10=compare(0.002,10)
#absdiff_0.002matches_Lmin20=compare(0.002,20)
#absdiff_0.002matches_Lmin40=compare(0.002,40)
absdiff_0.005matches_Lmin5=compare(0.005,5)
absdiff_0.005matches_Lmin10=compare(0.005,10)
absdiff_0.005matches_Lmin20=compare(0.005,20)
absdiff_0.01matches_Lmin5=compare(0.01,5)
absdiff_0.01matches_Lmin10=compare(0.01,10)
absdiff_0.01matches_Lmin20=compare(0.01,20)
absdiff_0.025matches_Lmin5=compare(0.025,5)
absdiff_0.025matches_Lmin10=compare(0.025,10)
absdiff_0.025matches_Lmin20=compare(0.025,20)
absdiff_0.05matches_Lmin5=compare(0.05,5)
absdiff_0.05matches_Lmin10=compare(0.05,10)
absdiff_0.05matches_Lmin20=compare(0.05,20)
absdiff_0.1matches_Lmin5=compare(0.1,5)
absdiff_0.1matches_Lmin10=compare(0.1,10)
absdiff_0.1matches_Lmin20=compare(0.1,20)

mean_use2<-vector()
data2=absdiff_0.025matches_Lmin10[[1]]
lengthall=length(absdiff_0.025matches_Lmin10[[6]])
for(i in 1:lengthall){
  mean_use2[i]=(data2[i]+data2[i+lengthall]+data2[i+lengthall*2]+data2[i+lengthall*3])/4
}
g1=ggplot(data=NULL,aes(x=absdiff_0.05matches_Lmin10[[6]],y=mean_use))+
  geom_point()+xlab("entropy")+ylab("absolute difference")+ggtitle("Sparsity=5%")
g2=ggplot(data=NULL,aes(x=absdiff_0.025matches_Lmin10[[6]],y=mean_use2))+
  geom_point()+xlab("entropy")+ylab("absolute difference")+ggtitle("Sparsity=2.5%")
library(gridExtra)
grid.arrange(g1,g2,nrow=1)


data_all=data.frame(mean=c(mean(absdiff_0.005matches_Lmin5[[1]]),
                           mean(absdiff_0.005matches_Lmin10[[1]]),
                           mean(absdiff_0.005matches_Lmin20[[1]]),
                           mean(absdiff_0.01matches_Lmin5[[1]]),
                           mean(absdiff_0.01matches_Lmin10[[1]]),
                           mean(absdiff_0.01matches_Lmin20[[1]]),
                           mean(absdiff_0.025matches_Lmin5[[1]]),
                           mean(absdiff_0.025matches_Lmin10[[1]]),
                           mean(absdiff_0.025matches_Lmin20[[1]]),
                           mean(absdiff_0.05matches_Lmin5[[1]]),
                           mean(absdiff_0.05matches_Lmin10[[1]]),
                           mean(absdiff_0.05matches_Lmin20[[1]]),
                           mean(absdiff_0.1matches_Lmin5[[1]]),
                           mean(absdiff_0.1matches_Lmin10[[1]]),
                           mean(absdiff_0.1matches_Lmin20[[1]])),
                    median=c(median(absdiff_0.005matches_Lmin5[[1]]),
                             median(absdiff_0.005matches_Lmin10[[1]]),
                             median(absdiff_0.005matches_Lmin20[[1]]),
                             median(absdiff_0.01matches_Lmin5[[1]]),
                             median(absdiff_0.01matches_Lmin10[[1]]),
                             median(absdiff_0.01matches_Lmin20[[1]]),
                             median(absdiff_0.025matches_Lmin5[[1]]),
                             median(absdiff_0.025matches_Lmin10[[1]]),
                             median(absdiff_0.025matches_Lmin20[[1]]),
                             median(absdiff_0.05matches_Lmin5[[1]]),
                             median(absdiff_0.05matches_Lmin10[[1]]),
                             median(absdiff_0.05matches_Lmin20[[1]]),
                             median(absdiff_0.1matches_Lmin5[[1]]),
                             median(absdiff_0.1matches_Lmin10[[1]]),
                             median(absdiff_0.1matches_Lmin20[[1]])),
                    quantile99=c(quantile(absdiff_0.005matches_Lmin5[[1]],0.99),
                                 quantile(absdiff_0.005matches_Lmin10[[1]],0.99),
                                 quantile(absdiff_0.005matches_Lmin20[[1]],0.99),
                                 quantile(absdiff_0.01matches_Lmin5[[1]],0.99),
                                 quantile(absdiff_0.01matches_Lmin10[[1]],0.99),
                                 quantile(absdiff_0.01matches_Lmin20[[1]],0.99),
                                 quantile(absdiff_0.025matches_Lmin5[[1]],0.99),
                                 quantile(absdiff_0.025matches_Lmin10[[1]],0.99),
                                 quantile(absdiff_0.025matches_Lmin20[[1]],0.99),
                                 quantile(absdiff_0.05matches_Lmin5[[1]],0.99),
                                 quantile(absdiff_0.05matches_Lmin10[[1]],0.99),
                                 quantile(absdiff_0.05matches_Lmin20[[1]],0.99),
                                 quantile(absdiff_0.1matches_Lmin5[[1]],0.99),
                                 quantile(absdiff_0.1matches_Lmin10[[1]],0.99),
                                 quantile(absdiff_0.1matches_Lmin20[[1]],0.99)),
                    maxmatchfrac=c(mean(absdiff_0.005matches_Lmin5[[2]]),
                                   mean(absdiff_0.005matches_Lmin10[[2]]),
                                   mean(absdiff_0.005matches_Lmin20[[2]]),
                                   mean(absdiff_0.01matches_Lmin5[[2]]),
                                   mean(absdiff_0.01matches_Lmin10[[2]]),
                                   mean(absdiff_0.01matches_Lmin20[[2]]),
                                   mean(absdiff_0.025matches_Lmin5[[2]]),
                                   mean(absdiff_0.025matches_Lmin10[[2]]),
                                   mean(absdiff_0.025matches_Lmin20[[2]]),
                                   mean(absdiff_0.05matches_Lmin5[[2]]),
                                   mean(absdiff_0.05matches_Lmin10[[2]]),
                                   mean(absdiff_0.05matches_Lmin20[[2]]),
                                   mean(absdiff_0.1matches_Lmin5[[2]]),
                                   mean(absdiff_0.1matches_Lmin10[[2]]),
                                   mean(absdiff_0.1matches_Lmin20[[2]])),
                    time=c(absdiff_0.005matches_Lmin5[[3]],
                           absdiff_0.005matches_Lmin10[[3]],
                           absdiff_0.005matches_Lmin20[[3]],
                           absdiff_0.01matches_Lmin5[[3]],
                           absdiff_0.01matches_Lmin10[[3]],
                           absdiff_0.01matches_Lmin20[[3]],
                           absdiff_0.025matches_Lmin5[[3]],
                           absdiff_0.025matches_Lmin10[[3]],
                           absdiff_0.025matches_Lmin20[[3]],
                           absdiff_0.05matches_Lmin5[[3]],
                           absdiff_0.05matches_Lmin10[[3]],
                           absdiff_0.05matches_Lmin20[[3]],
                           absdiff_0.1matches_Lmin5[[3]],
                           absdiff_0.1matches_Lmin10[[3]],
                           absdiff_0.1matches_Lmin20[[3]]),
                    entropy_match=c(mean(absdiff_0.005matches_Lmin5[[4]]),
                                     mean(absdiff_0.005matches_Lmin10[[4]]),
                                     mean(absdiff_0.005matches_Lmin20[[4]]),
                                     mean(absdiff_0.01matches_Lmin5[[4]]),
                                     mean(absdiff_0.01matches_Lmin10[[4]]),
                                     mean(absdiff_0.01matches_Lmin20[[4]]),
                                     mean(absdiff_0.025matches_Lmin5[[4]]),
                                     mean(absdiff_0.025matches_Lmin10[[4]]),
                                     mean(absdiff_0.025matches_Lmin20[[4]]),
                                     mean(absdiff_0.05matches_Lmin5[[4]]),
                                     mean(absdiff_0.05matches_Lmin10[[4]]),
                                     mean(absdiff_0.05matches_Lmin20[[4]]),
                                     mean(absdiff_0.1matches_Lmin5[[4]]),
                                     mean(absdiff_0.1matches_Lmin10[[4]]),
                                     mean(absdiff_0.1matches_Lmin20[[4]])),
                    entropy_unmatch=c(mean(absdiff_0.005matches_Lmin5[[5]]),
                                       mean(absdiff_0.005matches_Lmin10[[5]]),
                                       mean(absdiff_0.005matches_Lmin20[[5]]),
                                       mean(absdiff_0.01matches_Lmin5[[5]]),
                                       mean(absdiff_0.01matches_Lmin10[[5]]),
                                       mean(absdiff_0.01matches_Lmin20[[5]]),
                                       mean(absdiff_0.025matches_Lmin5[[5]]),
                                       mean(absdiff_0.025matches_Lmin10[[5]]),
                                       mean(absdiff_0.025matches_Lmin20[[5]]),
                                       mean(absdiff_0.05matches_Lmin5[[5]]),
                                       mean(absdiff_0.05matches_Lmin10[[5]]),
                                       mean(absdiff_0.05matches_Lmin20[[5]]),
                                       mean(absdiff_0.1matches_Lmin5[[5]]),
                                       mean(absdiff_0.1matches_Lmin10[[5]]),
                                       mean(absdiff_0.1matches_Lmin20[[5]])),
                    minmatchfrac=factor(c(0.5,0.5,0.5,1,1,1,2.5,2.5,2.5,5,5,5,10,10,10),levels=c(0.5,1,2.5,5,10)),
                    L_minmatch=factor(rep(c(5,10,20),5),levels=c(5,10,20))
                    )

library(ggplot2)
a=ggplot(data=data_all,aes(x=L_minmatch,y=mean,colour=minmatchfrac))+geom_point(size=4)+theme_bw()+
  ggtitle("(a) Average absolute painting difference")+
  ylab("Mean")+xlab("Minimum length of matches")+
  scale_color_discrete(name="minmatchfrac(%)")+
  theme(
    title=element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    legend.title = element_text(size=14),
    legend.text = element_text(size=14))
b=ggplot(data=data_all,aes(x=L_minmatch,y=median,colour=minmatchfrac))+geom_point(size=4)+theme_bw()+
  ggtitle("(b) Median of absolute painting difference")+
  ylab("Median")+xlab("Minimum length of matches")+
  scale_color_discrete(name="minmatchfrac(%)")+
  theme(
    title=element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    legend.title = element_text(size=14),
    legend.text = element_text(size=14))
c=ggplot(data=data_all,aes(x=L_minmatch,y=quantile99,colour=minmatchfrac))+geom_point(size=4)+theme_bw()+
  ggtitle("(c) 99% quantile of absolute painting difference")+
  ylab("99% quantile")+xlab("Minimum length of matches")+
  scale_color_discrete(name="minmatchfrac(%)")+
  theme(title=element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    legend.title = element_text(size=14),
    legend.text = element_text(size=14))
d=ggplot(data=data_all,aes(x=L_minmatch,y=maxmatchfrac,colour=minmatchfrac))+geom_point(size=4)+theme_bw()+
  ggtitle("(d) Average max-match fraction")+
  ylab("Average max-match fraction")+xlab("Minimum length of matches")+
  scale_color_discrete(name="minmatchfrac(%)")+
  theme(title=element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    legend.title = element_text(size=14),
    legend.text = element_text(size=14))

data_use=data.frame(rbind(data_all[,8:9],data_all[,8:9]),
                    entropy=c(data_all[,6],data_all[,7]),
                    entropy_type=factor(rep(c("max_match","max_unmatch"),each=15),levels=c("max_match","max_unmatch")))

e=ggplot(data=data_use,aes(x=L_minmatch,y=entropy,colour=minmatchfrac,shape=entropy_type))+
  geom_point(size=4)+theme_bw()+
  ggtitle("(e) Entropy of SNPs for d-pbwt painting")+
  ylab("Entropy")+xlab("Minimum length of matches")+
  scale_color_discrete(name="minmatchfrac(%)")+scale_shape_discrete(name="Entropy_type")+
  theme(title=element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14))

f=ggplot(data=data_all,aes(x=L_minmatch,y=time,colour=minmatchfrac))+geom_point(size=4)+theme_bw()+
  ggtitle("(f) Running time (parallel) of d-pbwt painting")+
  ylab("Time (s)")+xlab("Minimum length of matches")+
  scale_color_discrete(name="minmatchfrac(%)")+
  theme(title=element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14))

library(lemon)
grid_arrange_shared_legend(a, b, c, d, e, f, nrow = 2, ncol = 3, position = "bottom")



  setwd("C:/Ubuntu/longmatchquery/bigeg2")
  starttime=proc.time()
  paintings=paintingalldense(refindex,nind=50,targetfrac=1.0,fixrho=FALSE,method="Viterbi",
                             minmatchfrac = 0.01,L_initial=160,L_minmatch=10)
  diff_all<-vector()
  p_chromo_all <- vector()
  p_pbwt_all <- vector()
  reliability <- vector()
  setwd("C:/Ubuntu/longmatchquery/bigeg2/chromopainter")
  for(i in 1:100){
    print(i)
    p_chromo=read.table(paste0("test_",i,".txt"))[,-1]
    p_pbwt=as.data.frame(paintings[[i]])[,-1]
    #reliability=c(reliability,p_pbwt[,1])
    diff=unlist(abs(p_chromo-p_pbwt))
    diff_all=c(diff_all,diff)
    reliability <- c(reliability,as.data.frame(paintings[[i]])[,1])
    p_chromo_all <- c(p_chromo_all,unlist(p_chromo))
    p_pbwt_all <- c(p_pbwt_all,unlist(p_pbwt))
  }

mean_use2<-vector()
lengthall=length(reliability)
for(i in 1:lengthall){
  mean_use2[i]=(diff_all[i]+diff_all[i+lengthall]+diff_all[i+lengthall*2]+diff_all[i+lengthall*3])/4
}
plot(reliability,mean_use2)
library(ggplot2)
ggplot(data=NULL,aes(x=reliability,y=mean_use2))+
  geom_point()+xlab("number of matches")+ylab("average absolute difference")+ggtitle("Sparsity>=2.5%")





p_FLARE=read.table("flare.out.anc.vcf")[,10:59]

paintings_FLARE <- vector("list", length = 100)

for (i in 1:length(paintings_FLARE)) {
  paintings_FLARE[[i]] <- data.frame(matrix(ncol = 4, nrow = 2091))
}

for(i in 1:50){
  print(i)
  for(j in 1:2091){
    probs <- as.numeric(regmatches(p_FLARE[j,i], gregexpr("[0-9.]+", p_FLARE[j,i]))[[1]][5:12])
    paintings_FLARE[[2*i-1]][j,]=probs[1:4]
    paintings_FLARE[[2*i]][j,]=probs[5:8]
  }
}


setwd("C:/Ubuntu/longmatchquery/bigeg2")
starttime=proc.time()
paintings=paintingalldense(refindex,nind=50,targetfrac=1.0,fixrho=FALSE,method="Viterbi",
                           minmatchfrac = 0.01,L_initial=160,L_minmatch=10)
proc.time()-starttime
p_chromo_all <- vector()
p_pbwt_all <- vector()
p_FLARE_all <- vector()
setwd("C:/Ubuntu/longmatchquery/bigeg2/chromopainter")
for(i in 1:100){
  print(i)
  p_chromo=read.table(paste0("test2_",i,".txt"))[,-1]
  p_pbwt=as.data.frame(paintings[[i]])[,-1]
  p_FLARE_temp=as.data.frame(paintings_FLARE[[i]])
  p_chromo_all <- c(p_chromo_all,unlist(p_chromo))
  p_pbwt_all <- c(p_pbwt_all,unlist(p_pbwt))
  p_FLARE_all <- c(p_FLARE_all,unlist(p_FLARE_temp))
}

diff_pbwtandcp=abs(p_chromo_all-p_pbwt_all)
diff_FLAREandcp=abs(p_chromo_all-p_FLARE_all)
summary(diff_pbwtandcp)
summary(diff_FLAREandcp)
cor(p_chromo_all,p_pbwt_all)
cor(p_FLARE_all,p_pbwt_all)
