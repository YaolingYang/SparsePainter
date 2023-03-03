
setwd("C:/Ubuntu/longmatchquery")
library("Rcpp")
sourceCpp(file="hashmap_Rcpp.cpp")
setwd("C:/Ubuntu/longmatchquery/bigeg")
map=read.table('p.map')[,3:4]
gd=map[,2]/100000000
refindex=rep(c(0,1),each=10000)

painting=paintingalldense(gd,refindex,50,fixrho=FALSE)
painting=paintingalldense(gd,refindex,50,prop=0.0005)
painting=paintingalldense(gd,refindex,50,ite_time=10,method="EM",prop=0.0005)

aa=hashmaptest()
bb=hashmaptest2()
q=t(data.frame(bb))


est_rho_Viterbi (c(0,0,1,1,2,2,3,3), c(1,3,5,5,6,5,7,7),8, 0.01)

