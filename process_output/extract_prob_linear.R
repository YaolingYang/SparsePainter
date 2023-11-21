library(data.table)
linear=fread("chr19_1000G_linear_10inds_prob.txt.gz",fill=TRUE)
nsnp=as.numeric(linear[nrow(linear),SNPidx])
npop=ncol(linear)-1
split_linear=which(is.na(linear[,2]))
nind=length(split_linear)/2
split_linear=c(split_linear,nrow(linear)+1)
SNPidx=1:nsnp


painting=lapply(1:nind,function(h){
  print(h)
  painting_temp=matrix(0,ncol=nsnp,nrow=npop)
  snppos=as.numeric(unlist(linear[(split_linear[2*h-1]+1):(split_linear[2*h]-1),1]))
  ii=2 # Index of the end of the current interval
  for(i in 1:length(SNPidx)){
    SNPid = SNPidx[i]
    while(snppos[ii]<=SNPid & snppos[ii]!=nsnp){
      ii = ii + 1
    }
    if(snppos[ii]==nsnp){
      prob=as.numeric(linear[split_linear[2*h-1]+ii-1,-1])
      painting_temp[,i] = 0.5*prob/sum(prob)
    }else{
      prop=(SNPid - snppos[ii-1])/(snppos[ii]-snppos[ii-1])
      start=as.numeric(linear[split_linear[2*h-1]+ii-1,-1])
      end=as.numeric(linear[split_linear[2*h-1]+ii,-1])
      prob=start + (end-start)*prop
      painting_temp[,i] = 0.5*prob/sum(prob)
    }
  }
  
  snppos=as.numeric(unlist(linear[(split_linear[2*h]+1):(split_linear[2*h+1]-1),1]))
  ii=2
  for(i in 1:length(SNPidx)){
    SNPid = SNPidx[i]
    while(snppos[ii]<=SNPid & snppos[ii]!=nsnp){
      ii = ii + 1
    }
    if(snppos[ii]==nsnp){
      prob=as.numeric(linear[split_linear[2*h]+ii-1,-1])
      painting_temp[,i] = 0.5*prob/sum(prob)+painting_temp[,i]
    }else{
      prop=(SNPid - snppos[ii-1])/(snppos[ii]-snppos[ii-1])
      start=as.numeric(linear[split_linear[2*h]+ii-1,-1])
      end=as.numeric(linear[split_linear[2*h]+ii,-1])
      prob=start + (end-start)*prop
      painting_temp[,i] = 0.5*prob/sum(prob)+painting_temp[,i]
    }
  }
  painting_temp
})

big_matrix <- do.call(cbind, painting)
painting_final <- lapply(1:npop, function(x){
  matrix(big_matrix[x, ], ncol = nsnp)
})
