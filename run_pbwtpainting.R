source("painting_functions.R")
library(LDAandLDAS)

# N should be extracted from popnames.ids and diploid
# n_ref_each should be extracted from popnames.ids
# the other parameters should be written in a txt file,
# which should be given by user or using the default value

#compute painting
painting_all = cal_painting_all(N=50,n_refind_each=c(20,20),map=read.table('p.map')[,3:4],
                              method='Viterbi',fix_rho=TRUE,rate=1e-8,normalize=FALSE,
                              matchtype='donor',ite_time=20,
                              theta=NULL,theta_EM=FALSE,fix_theta=TRUE,diploid=TRUE)

#compute coancestry_matrix
coa_matrix = cal_coancestry(n_refind_each=c(20,20),map=read.table('p.map')[,3:4],
                            method='Viterbi',fix_rho=TRUE,rate=1e-8,
                            ite_time=20,theta=NULL,theta_EM=FALSE,fix_theta=TRUE,
                            diploid=TRUE)

# setting theta_EM=TRUE greatly reduces the speed...
coa_matrix2 = cal_coancestry(n_refind_each=c(20,20),map=read.table('p.map')[,3:4],
                             method='EM',fix_rho=TRUE,rate=1e-8,
                             ite_time=2,theta=NULL,theta_EM=TRUE,fix_theta=TRUE,
                             diploid=TRUE)

# if we greatly increase iteration times for EM, the coancestry matrix becomes closer
coa_matrix3 = cal_coancestry(n_refind_each=c(20,20),map=read.table('p.map')[,3:4],
                             method='EM',fix_rho=TRUE,rate=1e-8,
                             ite_time=60,theta=NULL,theta_EM=FALSE,fix_theta=TRUE,
                             diploid=TRUE)

#compute LDA and LDAS
lda=LDA(painting_all)
# gd should be in cM in function LDAS
map_use=data.frame(SNP=map[nrow(map):1,2],gd=gd[length(gd):1]*100) 
ldas = LDAS(lda,map_use,window=2)
plot(x=ldas$SNP/1e6,y=ldas$LDAS,xlab='Physical position (cM)',ylab='LDAS')
