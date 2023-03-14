## 15 reference haps, 10 snps, constructed so that 
## e.g. L=2 gives a complete match set but L=5 does not. 
## How far is the L=2 solution from the L=5 one? 
## How far are they both from the L=1 solution?

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

set.seed(49)
n_ref=15
L=10

ori_matchmat=matrix(rbinom(n_ref*L,1,0.5),nrow=n_ref,ncol=L)
matchlenmat=matrix_to_match_length(ori_matchmat)

#min_length=4

matchmat1=matchlenmat>=1
matchmat1[matchmat1]=1

matchmat=matchlenmat>=2
matchmat[matchmat]=1

matchmat2=matchlenmat>=4
matchmat2[matchmat2]=1

matchmat_allmatch=matchmat2
matchmat_allmatch[,6]=1

matchmat_mut=matchmat2
matchmat_mut[,6]=0.001

matchmat_difflength=matchmat2
matchmat_difflength[14,4:6]=1

sameprob=0.95
otherprob=(1-sameprob)/n_ref



marginalProb<-function(f,b){
  ### marginal probability
  marginal_prob <- f*b
  for(j in 1:L){
    marginal_prob[,j]=marginal_prob[,j]/sum(marginal_prob[,j])
  }
  marginal_prob
}
forwardProb2<-function(matchmat,sameprob,otherprob){
  n_ref=nrow(matchmat)
  L=ncol(matchmat)
  forward_prob2 <- matrix(0,nrow=n_ref,ncol=L)
  forward_prob2[,1]=matchmat[,1]/sum(matchmat[,1])
  for(j in 2:L){
    if(sum(matchmat[,j])==0){
      forward_prob2[,j]=matchmat[,j-1] * ( sameprob * forward_prob2[,j-1] + otherprob)
      forward_prob2[,j]=forward_prob2[,j]/sum(forward_prob2[,j])
    }else{
      forward_prob2[,j]=matchmat[,j] * ( sameprob * forward_prob2[,j-1] + otherprob)
      forward_prob2[,j]=forward_prob2[,j]/sum(forward_prob2[,j])
    }
    
  }
  forward_prob2
}
backwardProb2<-function(matchmat,sameprob,otherprob){
  ### backward algorithm 2
  n_ref=dim(matchmat)[1]
  L=dim(matchmat)[2]
  backward_prob2 <- matrix(0,nrow=n_ref,ncol=L)
  backward_prob2[,L]=1
  backward_prob2[,L]=1/sum(backward_prob2[,L])
  for(j in (L-1):1){
    if(sum(matchmat[,j+1])==0){
      Bjp1=matchmat[,j+2] * backward_prob2[,j+1]
      backward_prob2[,j] = otherprob * sum(Bjp1) + sameprob * Bjp1
      backward_prob2[,j]=backward_prob2[,j] /sum(backward_prob2[,j])
    }else{
      Bjp1=matchmat[,j+1] * backward_prob2[,j+1]
      backward_prob2[,j] = otherprob * sum(Bjp1) + sameprob * Bjp1
      backward_prob2[,j]=backward_prob2[,j] /sum(backward_prob2[,j])
    }

    
  }
  backward_prob2
}
forwardBackward2<-function(matchmat,sameprob,otherprob){
  f<-forwardProb2(matchmat,sameprob,otherprob)
  b<-backwardProb2(matchmat,sameprob,otherprob)
  marginalProb(f,b)
}

## the same as chromoPainter, i.e. minlen=1
marginal_prob_chromopainter<-forwardBackward2(matchmat1,sameprob,otherprob)
## matchmat[,j]=matchmat[,j-1] for forward, matchmat[,j]=matchmat[,j+1] for backward
marginal_prob_skip<-forwardBackward2(matchmat2,sameprob,otherprob)
## minlen=2
marginal_prob_original<-forwardBackward2(matchmat,sameprob,otherprob)
## assume all are matches for the no-match SNP
marginal_prob_allmatch<-forwardBackward2(matchmat_allmatch,sameprob,otherprob)
## assume all are mutations for the no-match SNP
marginal_prob_mut<-forwardBackward2(matchmat_mut,sameprob,otherprob)
## use different minlen, i.e. the idea of L_use for dpbwt
marginal_prob_difflength<-forwardBackward2(matchmat_difflength,sameprob,otherprob)

diff=marginal_prob_original-marginal_prob_allmatch

diff2=marginal_prob_original-marginal_prob_difflength

diff3=marginal_prob_original-marginal_prob_chromopainter

diff4=marginal_prob_skip-marginal_prob_chromopainter

diff5=marginal_prob_mut-marginal_prob_chromopainter

