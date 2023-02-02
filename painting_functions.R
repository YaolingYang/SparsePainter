get_gd <- function(map,rate){
  # if the map file doesn't have a recombination distance
  if(all(!map[,1])){
    if(is.null(rate)){
      rate=1e-8
    }
    return(map[,2]*rate)
  }else{
    return(map[,1])
  }
}

## assume we get a single file called matchdata
## It includes 4 columns, which are ref_ind, mod_ind, match_start and match_end
## the ref_ind ranges from 1 to sum(n_ref_ancestry)
data_process <- function(matchdata,gd){
  # matchdata is the original match data
  # gd is the genetic position of each SNP in Morgan
  colnames(matchdata)=c('ref_ind','mod_ind','start','end')
  matchdata[,2]=as.numeric(gsub('q','',matchdata[,2]))+1
  #if(!is.null(i)){
  #  matchdata=matchdata[which(matchdata[,2]==i),]
  #}
  matchdata[,1]=matchdata[,1]+1
  matchdata[,3]=as.numeric(gsub('[','',gsub(',','',matchdata[,3]),fixed=TRUE))+1
  matchdata[,4]=as.numeric(gsub(')','',matchdata[,4]))
  matchdata=data.frame(matchdata,length=gd[matchdata$end]-gd[matchdata$start])
  return(matchdata)
}

#estimate \rho
est_rho <- function(match,gd){
  # match is the match data after data_process
  # i is the index of target individual
  #data=match[which(match$mod_ind==i),c(3,4)]
  data=match[,c(3,4)]
  nsnp=length(gd)
  # recombination position
  rec_pos=vector()
  j=1
  while(j<nsnp){
    data_use=data[which(data$start<=j & data$end>=j),]
    if(nrow(data_use)==0){
      rec_pos=c(rec_pos,j+1)
      j=j+1
    }else{
      if(max(data_use$end)<nsnp){
        rec_pos=c(rec_pos,max(data_use$end)+1)
      }
      j=max(data_use$end)+1
    }
  }
  n_rec=length(rec_pos)
  return(n_rec/(gd[nsnp_use]-gd[1]))
}

# a logical matrix describing whether the ith target individual match
# each reference samples at each SNP
matchmatrix <- function(match,n_ref,nsnp){
  matchmat=matrix(FALSE,nrow=n_ref,ncol=nsnp)
  for(j in 1:nsnp){
    #data=match[which(match$mod_ind==i),]
    data=match[which(match$start<=j & match$end>=j),]
    matchmat[data$ref_ind,j]=TRUE
  }
  return(matchmat)
}

# the probability for transitting to the same state
cal_sameprob <- function(nsnp,rho,gd){
  sameprob=vector()
  for(j in 1:(nsnp-1)){
    gd_use=(gd[j+1]-gd[j])*100
    sameprob[j]=exp(-rho*gd_use)
  }
  return(sameprob)
}

# calculate the forward probability for target individual i
cal_forward <- function(matchmat,nsnp,n_ref,sameprob){
  otherprob = (1-sameprob)/n_ref
  forward_prob <- matrix(0,nrow=n_ref,ncol=nsnp)
  if(sum(matchmat[,1])!=0){
    forward_prob[,1] = as.numeric(matchmat[,1])/sum(matchmat[,1])
  }else{
    forward_prob[,1] = 1/n_ref
  }
  
  for(j in 2:nsnp){
    if(sum(matchmat[,j])==0){
      matchmat[,j]=1
    }
    forward_prob[,j] = matchmat[,j] * (sameprob[j-1] * forward_prob[,j-1] + otherprob[j-1])
    forward_prob[,j] = forward_prob[,j]/sum(forward_prob[,j])
  }
  return(forward_prob)
}

# calculate the backward probability for target individual i
cal_backward <- function(matchmat,nsnp,n_ref,sameprob){
  otherprob=(1-sameprob)/n_ref
  backward_prob <- matrix(0,nrow=n_ref,ncol=nsnp)
  backward_prob[,nsnp]=1
  for(j in (nsnp-1):1){
    if(sum(matchmat[,j+1])==0){
      matchmat[,j+1]=1
    }
    Bjp1=matchmat[,j+1] * backward_prob[,j+1]
    backward_prob[,j] = otherprob[j] * sum(Bjp1) + sameprob[j] * Bjp1
    backward_prob[,j] = backward_prob[,j]/sum(backward_prob[,j])
  } 
  return(backward_prob)
}

# calculate the marginal probability for each reference
cal_marginal <- function(forward_prob,backward_prob){
  marginal_prob <- forward_prob*backward_prob
  for(j in 1:ncol(marginal_prob)){
    marginal_prob[,j] = marginal_prob[,j]/sum(marginal_prob[,j])
  }
  return(marginal_prob)
}

# calculate the probability for each ancestry
cal_painting <- function(marginal_prob,n_ref_each){
  # n_ref_each is a vector representing the number of reference for each ancestry
  painting = list()
  painting[[1]] = colSums(marginal_prob[1:n_ref_each[1],])
  start=n_ref_each[1]+1
  for(k in 2:length(n_ref_each)){
    painting[[k]] = colSums(marginal_prob[start:(start+n_ref_each[k]-1),])
    start=start+n_ref_each[k]
  }
  return(painting)
}

cal_painting_full <- function(match,gd,n_ref_each,
                              fix_rho=NULL,returneachref=FALSE){
  n_ref=sum(n_ref_each)
  
  nsnp=length(gd)
  
  #match=data_process(matchdata,gd,i)
  
  if(is.null(fix_rho)){
    rho=est_rho(match,i,gd)
  }else{
    rho=fix_rho
  }
  
  matchmat=matchmatrix(match,n_ref,nsnp)
  
  sameprob=cal_sameprob(nsnp,rho,gd)
  
  forward_prob = cal_forward(matchmat,nsnp,n_ref,sameprob)
  
  backward_prob = cal_backward(matchmat,nsnp,n_ref,sameprob)
  
  marginal_prob = cal_marginal(forward_prob,backward_prob)
  
  if(returneachref){
    return(marginal_prob)
  }else{
    return(cal_painting(marginal_prob,n_ref_each))
  }
}



