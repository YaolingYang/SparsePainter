library(data.table)
constant=fread("chr19_1000G_constant_10inds_prob.txt.gz",fill=TRUE)
nsnp=as.numeric(constant[nrow(constant),2])
npop=ncol(constant)-2
split_constant=which(is.na(constant[,2]))
nind=length(split_constant)/2
split_constant=c(split_constant,nrow(constant)+1)
SNPidx=1:nsnp

painting=lapply(1:nind,function(h){
  print(h)
  painting_temp=matrix(0,ncol=nsnp,nrow=npop)
  snppos=as.numeric(unlist(constant[(split_constant[2*h-1]+1):(split_constant[2*h]-1),2]))
  ii=1 # Index of the end of the current interval
  for(i in 1:length(SNPidx)){
    SNPid = SNPidx[i]
    while(snppos[ii]<SNPid){
      ii = ii + 1
    }
    prob=as.numeric(constant[split_constant[2*h-1]+ii,3:(npop+2)])
    painting_temp[,i] = 0.5*prob/sum(prob)
  }
  
  snppos=as.numeric(unlist(constant[(split_constant[2*h]+1):(split_constant[2*h+1]-1),2]))
  ii=1
  for(i in 1:length(SNPidx)){
    SNPid = SNPidx[i]
    while(snppos[ii]<SNPid){
      ii = ii + 1
    }
    prob=as.numeric(constant[split_constant[2*h]+ii,3:(npop+2)])
    painting_temp[,i] = 0.5*prob/sum(prob)+painting_temp[,i]
  }
  painting_temp
})

big_matrix <- do.call(cbind, painting)
painting_final <- lapply(1:npop, function(x){
  matrix(big_matrix[x, ], ncol = i)
})
