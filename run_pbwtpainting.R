source("painting_functions.R")
library(LDAandLDAS)

map=read.table('p1new.map')[,3:4]

gd=get_gd(map,rate=1e-8) ### gd is the genetic distance in Morgan

for(i in 1:1000){
  print(i)
  matchdata=read.table(paste0('match',i-1,'.txt'))[,c(1,3,5,6)]
  
  match=data_process(matchdata,gd)
  
  #match=match[which(match$mod_ind==1),]
  
  painting=cal_painting_full(match,gd,n_ref_each=c(40,40),
                             fix_rho=90.6536,returneachref=FALSE)
  
  n_ref=length(n_ref_each)
  if(i==1){
    painting_all=painting
  }else{
    for(j in 1:n_ref){
      painting_all[[j]]=rbind(painting_all[[j]],painting[[j]])
    }
    if(i==2){
      for(j in 1:n_ref){
        colnames(painting_all[[j]])=map[,2]
      }
    }
  }
}


#compute LDA and LDAS
lda=LDA(painting_all)
# gd should be in cM in function LDAS
map_use=data.frame(SNP=map[nrow(map):1,2],gd=gd[length(gd):1]*100) 
ldas = LDAS(lda,map_use,window=2)
plot(x=ldas$SNP/1e6,y=ldas$LDAS,xlab='Physical position (cM)',ylab='LDAS')
