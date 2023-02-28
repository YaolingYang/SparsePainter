library("Rcpp")
sourceCpp(file="test_Rcpp_hashmap.cpp")

set.seed(3)
L=1000
f=runif(L,0.1,0.5)
beta=20
## f = alpha/(alpha + beta)
## alpha f + beta f = alpha
## alpha = beta f/(1-f)
N=50000
n_ref=2*N
rho=1/500 # Recombination rate
min_length=10

popf1=rbeta(L,f/(1-f)*beta,beta)
popf2=rbeta(L,f/(1-f)*beta,beta)
pop1=sapply(popf1,rbinom,n=N,size=1)
pop2=sapply(popf2,rbinom,n=N,size=1)
pop=rbind(pop1,pop2)
breaks=cumsum(rpois(2*L*rho,1/rho))
breaks=c(0,breaks[breaks<L],L) # Breakpoints for the target
nbp=length(breaks)-1
pathway=sample(1:dim(pop)[1],nbp,TRUE)

target=rep(NA,L)
for(i in 1:nbp){
    target[breaks[i]:breaks[i+1]]=pop[pathway[i],breaks[i]:breaks[i+1]]
}


cal_O_matrix <- function(pop,j){
  O=matrix(0,nrow=dim(pop)[1],ncol=dim(pop)[1])
  diag(O)[matchmat[,j]]=1
  return(O)
}

vec_to_match_length<-function(x){
    r=rep(0,length(x))
    r[1]=x[1]
    for(i in 2:length(r)){
        if(x[i]) r[i]=r[i-1]+1
    }
    for(i in (length(r)-1):1){
        if(x[i+1] & x[i]) r[i]=r[i+1]
    }
    r
}
matrix_to_match_length<-function(m){
    t(apply(m,1,vec_to_match_length))
}


matchmat=t(t(pop)==target)
matchlenmat=matrix_to_match_length(matchmat)
matchmat=matchlenmat>min_length
matchsums=apply(matchmat,2,sum)
sapply(1:20,function(x)sum(matchlenmat>=x))/prod(dim(matchlenmat))


sameprob=exp(-rho)
otherprob=(1-sameprob)/n_ref


#image(1:dim(matchmat)[1],1:dim(matchmat)[2],matchlenmat,xlab="Individual",ylab="SNP",main="Match matrix")#,col=c("white","blue"))
#ylocs=as.numeric(rbind(breaks[1:nbp],breaks[1+(1:nbp)]))
#xlocs=as.numeric(rep(pathway,each=2))
#lines(xlocs,ylocs-0.5,lwd=2)

##################3
#trans = matrix(otherprob,nrow=n_ref,ncol=n_ref)
#diag(trans) = diag(trans) + sameprob
############## FULL ALGORITHM
## forward algorithm
#forward_prob <- matrix(0,nrow=n_ref,ncol=L)
#forward_prob[,1]=diag(cal_O_matrix(pop,1))
#forward_prob[,1]=forward_prob[,1]/sum(forward_prob[,1])
#for(j in 2:L){
#    forward_prob[,j]=cal_O_matrix(pop,j) %*% trans %*% forward_prob[,j-1]
#    forward_prob[,j]=forward_prob[,j] /sum(forward_prob[,j])
#}
  
### backward algorithm
#backward_prob <- matrix(0,nrow=n_ref,ncol=L)
#backward_prob[,L]=1
#for(j in (L-1):1){
#    backward_prob[,j]=trans %*% cal_O_matrix(pop,j+1) %*%backward_prob[,j+1]
#    backward_prob[,j]=backward_prob[,j] /sum(backward_prob[,j])
#} 


## ### backward algorithm (broken)
## backward_prob <- matrix(0,nrow=n_ref,ncol=L)
## backward_prob[,L]=1/n_ref
## for(j in (L-1):1){
##     backward_prob[,j]= cal_O_matrix(pop,j) %*% trans %*%backward_prob[,j+1]
##     backward_prob[,j]=backward_prob[,j]/sum(backward_prob[,j])
## }

### marginal probability
#marginal_prob <- forward_prob*backward_prob
#for(j in 1:L){
#    marginal_prob[,j]=marginal_prob[,j] /sum(marginal_prob[,j])
#}

## Plot of matches and probabilities
#par(mfrow=c(1,2))
#image(1:dim(matchmat)[1],1:dim(matchmat)[2],matchlenmat,xlab="Individual",ylab="SNP",main="Match matrix")
#ylocs=as.numeric(rbind(breaks[1:nbp],breaks[1+(1:nbp)]))
#xlocs=as.numeric(rep(pathway,each=2))
#lines(xlocs,ylocs-0.5,lwd=2,col="blue")

#image(1:dim(matchmat)[1],1:dim(matchmat)[2],marginal_prob,xlab="Individual",ylab="SNP",main="Posterior matrix")
#ylocs=as.numeric(rbind(breaks[1:nbp],breaks[1+(1:nbp)]))
#xlocs=as.numeric(rep(pathway,each=2))
#lines(xlocs,ylocs-0.5,lwd=2,col="blue")


marginalProb<-function(f,b){
### marginal probability
    marginal_prob <- f*b
    for(j in 1:L){
        marginal_prob[,j]=marginal_prob[,j]/sum(marginal_prob[,j])
    }
    marginal_prob
}
########## TEST ALGORITHM 2
## Only compute probabilities passing from valid to valid sets
## forward algorithm 2
forwardProb2<-function(matchmat,sameprob,otherprob){
    forward_prob2 <- matrix(0,nrow=n_ref,ncol=L)
    forward_prob2[,1]=matchmat[,1]/sum(matchmat[,1])
    for(j in 2:L){
        forward_prob2[,j]=matchmat[,j] * ( sameprob * forward_prob2[,j-1] + otherprob)
        forward_prob2[,j]=forward_prob2[,j]/sum(forward_prob2[,j])
    }
    forward_prob2
}
backwardProb2<-function(matchmat,sameprob,otherprob){
### backward algorithm 2
    backward_prob2 <- matrix(0,nrow=n_ref,ncol=L)
    backward_prob2[,L]=1
    backward_prob2[,L]=1/sum(backward_prob2[,L])
    for(j in (L-1):1){
        Bjp1=matchmat[,j+1] * backward_prob2[,j+1]
        backward_prob2[,j] = otherprob * sum(Bjp1) + sameprob * Bjp1
        backward_prob2[,j]=backward_prob2[,j] /sum(backward_prob2[,j])
    }
    backward_prob2
}
forwardBackward2<-function(matchmat,sameprob,otherprob){
    f<-forwardProb2(matchmat,sameprob,otherprob)
    b<-backwardProb2(matchmat,sameprob,otherprob)
    marginalProb(f,b)
}
############### END TEST ALGORITHM 2
########## TEST ALGORITHM 3
## Merge O and T matrices
forwardProb3<-function(matchmat,sameprob,otherprob){
## forward algorithm 3
    n_ref=dim(matchmat)[1]
    L=dim(matchmat)[2]
    forward_prob3 <- matrix(0,nrow=n_ref,ncol=L)
    forward_prob3[,1]=matchmat[,1]/sum(matchmat[,1])
    for(j in 2:L){
        tw=which(matchmat[,j])
        forward_prob3[tw,j] = sameprob * forward_prob3[tw,j-1] + otherprob 
        forward_prob3[matchmat[,j],j]= forward_prob3[matchmat[,j],j]/sum(forward_prob3[matchmat[,j],j])
    }
    forward_prob3
}
backwardProb3<-function(matchmat,sameprob,otherprob){
### backward algorithm 3
    n_ref=dim(matchmat)[1]
    L=dim(matchmat)[2]
    backward_prob3 <- matrix(otherprob,nrow=n_ref,ncol=L)
    backward_prob3[,L]=1
    backward_prob3[,L]=1/sum(backward_prob3[,L])
    for(j in (L-1):1){
        tw=which(matchmat[,j+1])
        #nm=length(tw)
        Bjp1=backward_prob3[tw,j+1]
        sumBjp1=sum(Bjp1)
        #missingprob= (dim(backward_prob3)[1]-nm)*otherprob*sumBjp1
        val=otherprob * sumBjp1 + sameprob * Bjp1
        #presentprob=sum(val)
        ### presentprob+missingprob == sumBjp1
        #backward_prob3[tw,j] = val /(presentprob+missingprob)
        backward_prob3[tw,j] = val /sumBjp1
    }
    backward_prob3
}
forwardBackward3<-function(matchmat,sameprob,otherprob){
    f<-forwardProb3(matchmat,sameprob,otherprob)
    b<-backwardProb3(matchmat,sameprob,otherprob)
    marginalProb(f,b)
}
############### END TEST ALGORITHM 3
############### TEST ALGORITHM 4: Sparse implementation
library(hash) # provides hash. std::unordered_map/sets in C++,
hVec<-function(default=0){
    ## Create an empty hash vector
    v=list(h=hash(),default=default)
    class(v)="hVec"
    v
}
hVecIdx<-function(v){
    ## Extract the indexes (keys) stored of the hash vector v
    as.numeric(names(v$h))
}
hVecIn<-function(v,idx){
    ## Return a binary vector of whether the indices idx are stored in v
    idxk=make.keys(idx) # convert values to characters
    has.key(idxk,v$h) # whether the key idxk exists in hash table v$h
}
hVecVal<-function(v,idx){
    ## Get the value of the indices idx of v
    if(length(idx)==0) return(numeric())
    idxk=make.keys(idx)
    hk=has.key(idxk,v$h)
    ret=rep(v$default,length(idxk))
    if(length(hk)>0) ret[hk]<-as.numeric(values(v$h,idxk[hk]))
    ret
    ## when the value of indices doesn't exist, use default
}
hVecSet<-function(v,idx,val){
    ## Set the indices idx of v to the vector val
    idxk=make.keys(idx)
    v$h[idxk]=val
    v
}
hVecCirc<-function(v1,v2){
    ## v1 \circ v2
    ret=hVec(v1$default*v2$default)
    v1idx=hVecIdx(v1)
    v2idx=hVecIdx(v2)
    idx=union(v1idx,v2idx)    
    v1val=hVecVal(v1,idx)
    v2val=hVecVal(v2,idx)
    hVecSet(ret,idx,v1val*v2val)
    #ret
}
hVecDense<-function(x,nrow){
    ## Convert a hash vector x into a dense vector of length nrow (the maximum array index)
    v=rep(x$default,nrow)
    vval=values(x$h)
    v[as.numeric(names(vval))]=as.numeric(vval)
    v

}
hVecSum<-function(v,idx=NULL){
    ## Compute the column sum of a hash vector
    ## NOTE: CURRENTLY ASSUMES default=0!!!
    if(all(is.null(idx))) idx=hVecIdx(v)
    sum(hVecVal(v,idx))
}
hVecScale<-function(v,x){
    ## Scale a hash vector by x
    idx=hVecIdx(v)
    val=hVecVal(v,idx)
    v$default=v$default*x
    hVecSet(v,idx,val*x)
}
matchmat2hMatrix<-function(m,default=0){
    ## Convert a binary matrix into a hash matrix
    r=apply(m,2,function(x){
        idx= which(x) ## positions of matches
        r=hVec(default)
        hVecSet(r,idx,1)
    },simplify=FALSE)
    attr(r,"matdim")=dim(m) # convert to a matrix with the same dimension as m
    class(r)="hMatrix"
    r
}
hMatrix2matrix<-function(m){
    ## Convert a hash matrix into a dense matrix
    nrow=attr(m,"matdim")[1]
    ret=sapply(m,hVecDense,nrow=nrow)
    ret
}
hMatrix<-function(dim,default=0){
    ## Create a hash matrix with the specified dimensions
    r=lapply(1:dim[2],function(x)hVec(default))
    attr(r,"matdim")=dim
    class(r)="hMatrix"
    r
}
forwardProb4<-function(m,sameprob,otherprob){
## Merge O and T matrices
## forward algorithm 4
    forward_prob4 <- hMatrix(attr(m,"matdim"))
    twj=hVecIdx(m[[1]])
    forward_prob4[[1]]=hVecSet(forward_prob4[[1]],twj,1/hVecSum(m[[1]],twj))
    for(j in 2:L){
        twj=hVecIdx(m[[j]])
        val=sameprob * hVecVal(forward_prob4[[j-1]],twj) + otherprob
        cs=sum(val)
        forward_prob4[[j]] = hVecSet(forward_prob4[[j]],twj, val/cs)
    }
    forward_prob4
}
backwardProb4<-function(m,sameprob,otherprob){
### backward algorithm 4
    backward_prob4 <- hMatrix(attr(m,"matdim"))
    ## First entry is not sparse, but all has the same value so no need to specify it explicitly
    backward_prob4[[L]]$default=1/attr(m,"matdim")[1]
    for(j in (L-1):1){
        twj=hVecIdx(m[[j+1]])
        #nm=length(twj)
        Bjp1=hVecVal(backward_prob4[[j+1]],twj)
        sumBjp1=sum(Bjp1)
        #missingelement=otherprob*sumBjp1
        #missingprob= (attr(m,"matdim")[1]-nm)* missingelement
        #val = missingelement + sameprob * Bjp1
        #presentprob=sum(val)
        #missingprob+presentprob==sumBjp1
        #hVecSet(backward_prob4[[j]],twj,val/(presentprob+missingprob))
        val = otherprob*sumBjp1 + sameprob * Bjp1
        backward_prob4[[j]]=hVecSet(backward_prob4[[j]],twj,val/sumBjp1)
        backward_prob4[[j]]$default=otherprob
    }
    backward_prob4
}
marginalProb4<-function(f,b){
### marginal probability
    marginal_prob4 <- hMatrix(attr(f,"matdim"))
    for(j in 1:L){
        marginal_prob4[[j]]<-hVecCirc(f[[j]],
                                      b[[j]])
        marginal_prob4[[j]]<-hVecScale(marginal_prob4[[j]],
                                       1.0/hVecSum(marginal_prob4[[j]]))
    }
    marginal_prob4
}
forwardBackward4<-function(matchmat,sameprob,otherprob){
    m <- matchmat2hMatrix(matchmat)
    forward_prob4<-forwardProb4(m,sameprob,otherprob)
    backward_prob4<-backwardProb4(m,sameprob,otherprob)
    marginal_prob4<-marginalProb4(forward_prob4,backward_prob4)
    hMatrix2matrix(marginal_prob4)
}

matchmatlist <- list()
for(i in 1:nrow(matchmat)){
  matchmatlist[[i]]=matchmat[i,]
}

#minmatchlen=min(apply(matchlenmat,2,max))
#m <- matchmat2hMatrix(matchmat)
## Timing comparison:
system.time(marginal_prob2<-forwardBackward2(matchmat,sameprob,otherprob))
system.time(marginal_prob3<-forwardBackward3(matchmat,sameprob,otherprob))
system.time(marginal_prob4<-forwardBackward4(matchmat,sameprob,otherprob))
system.time(marginal_prob<-forwardBackward(matchmatlist,sameprob,otherprob))



########################
## Plot of matches and probabilities
par(mfrow=c(1,2))
image(1:dim(matchmat)[1],1:dim(matchmat)[2],marginal_prob,xlab="Individual",ylab="SNP",main="Posterior matrix")
ylocs=as.numeric(rbind(breaks[1:nbp],breaks[1+(1:nbp)]))
xlocs=as.numeric(rep(pathway,each=2))
lines(xlocs,ylocs-0.5,lwd=2,col="blue")

image(1:dim(matchmat)[1],1:dim(matchmat)[2],marginal_prob2,xlab="Individual",ylab="SNP",main="Posterior matrix 2")
ylocs=as.numeric(rbind(breaks[1:nbp],breaks[1+(1:nbp)]))
xlocs=as.numeric(rep(pathway,each=2))
lines(xlocs,ylocs-0.5,lwd=2,col="blue")



forward_prob4test=hMatrix2matrix(forward_prob4)
backward_prob4test=hMatrix2matrix(backward_prob4)

marginal_prob4test=hMatrix2matrix(marginal_prob4)

marginal_prob4 <- forward_prob4test*backward_prob4test
for(j in 1:L){
    marginal_prob4[,j]=marginal_prob4[,j]/sum(marginal_prob4[,j])
}
############### END TEST ALGORITHM 4
