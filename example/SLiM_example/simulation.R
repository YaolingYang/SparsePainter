combinedata1 <- function(ncore){
  print(ncore)
  name1=paste('p_merged_filtered/stage7/test_p_merged_filtered_stage7_tmp_mainrun.linked_file1_ind',(1+25*(ncore-1)),'.copyprobsperlocus.out.gz',sep='')
  a <- read.table(name1,sep='')
  n_snp <- (nrow(a)-3)/2
  p1 <- c(a[2,3],a[c(3:(n_snp+2)),2])
  p1 <- rbind(p1,c(a[2,3],a[c((n_snp+4):(2*n_snp+3)),2]))
  for(i in 2:25){
    name=paste('p_merged_filtered/stage7/test_p_merged_filtered_stage7_tmp_mainrun.linked_file1_ind',(i+25*(ncore-1)),'.copyprobsperlocus.out.gz',sep='')
    a <- read.table(name,sep='')
    p1 <- rbind(p1,c(a[2,3],a[c(3:(n_snp+2)),2]))
    p1 <- rbind(p1,c(a[2,3],a[c((n_snp+4):(2*n_snp+3)),2]))
  }
  data1 <- as.data.frame(p1)
  
  colnames(data1)[1] <- 'ind'
  colnames(data1)[-1]<-a[3:(n_snp+2),1]
  return(data1)
}

combinedata2 <- function(ncore){
  print(ncore)
  name1=paste('p_merged_filtered/stage7/test_p_merged_filtered_stage7_tmp_mainrun.linked_file1_ind',(1+25*(ncore-1)),'.copyprobsperlocus.out.gz',sep='')
  a <- read.table(name1,sep='')
  n_snp <- (nrow(a)-3)/2
  p2 <- c(a[2,3],a[c(3:(n_snp+2)),3])
  p2 <- rbind(p2,c(a[2,3],a[c((n_snp+4):(2*n_snp+3)),3]))
  for(i in 2:25){
    name=paste('p_merged_filtered/stage7/test_p_merged_filtered_stage7_tmp_mainrun.linked_file1_ind',(i+25*(ncore-1)),'.copyprobsperlocus.out.gz',sep='')
    a <- read.table(name,sep='')
    p2 <- rbind(p2,c(a[2,3],a[c(3:(n_snp+2)),3]))
    p2 <- rbind(p2,c(a[2,3],a[c((n_snp+4):(2*n_snp+3)),3]))
  }
  data2 <- as.data.frame(p2)
  
  colnames(data2)[1] <- 'ind'
  colnames(data2)[-1]<-a[3:(n_snp+2),1]
  return(data2)
}
res <- lapply(1:20, combinedata1)
data1 <- rbind(res[[1]],res[[2]],res[[3]],res[[4]],res[[5]],
               res[[6]],res[[7]],res[[8]],res[[9]],res[[10]],
               res[[11]],res[[12]],res[[13]],res[[14]],res[[15]],
               res[[16]],res[[17]],res[[18]],res[[19]],res[[20]])
res <- lapply(1:20, combinedata2)
data2 <- rbind(res[[1]],res[[2]],res[[3]],res[[4]],res[[5]],
               res[[6]],res[[7]],res[[8]],res[[9]],res[[10]],
               res[[11]],res[[12]],res[[13]],res[[14]],res[[15]],
               res[[16]],res[[17]],res[[18]],res[[19]],res[[20]])
write.table(data1,file='painting_p1.csv',sep=',',row.names=FALSE)
write.table(data2,file='painting_p2.csv',sep=',',row.names=FALSE)