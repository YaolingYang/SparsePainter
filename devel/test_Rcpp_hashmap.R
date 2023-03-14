
setwd("C:/Ubuntu/longmatchquery")
library("Rcpp")
sourceCpp(file="hashmap.cpp")
setwd("C:/Ubuntu/longmatchquery/bigeg2")
map=read.table('p.map')[,3:4]
gd=map[,2]/100000000
refindex=rep(c(0,1,2,3),each=4000)

starttime=proc.time()
painting=paintingalldense(gd,refindex,1000,fixrho=FALSE,method="Viterbi",minmatchfrac = 0.001,L=80)
proc.time()-starttime

painting=paintingalldense(gd,refindex,50,fixrho=TRUE,indfrac=0.001,minmatchfrac = 0.001,L=80)
painting=paintingalldense(gd,refindex,50,ite_time=10,method="EM",
                          indfrac=0.001,minmatchfrac = 0.001,L=80)
painting=paintingalldense(gd,refindex,50,ite_time=10,method="EM",
                          indfrac=0.001,minsnpEM=100,EMsnpfrac=0.1,minmatchfrac = 0.001,L=80)

coancestry=chunklengthall(gd,refindex,minmatchfrac = 0.001,L=80)


