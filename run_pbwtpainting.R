source("painting_functions.R")
library(LDAandLDAS)

#later we read popnames.ids to get all information.
#compute painting
painting_all = cal_painting_all(N=1000,n_ref_each=c(40,40),map=read.table('p.map')[,3:4],
                              method='Viterbi',fix_rho=TRUE,rate=1e-8)

#compute coancestry_matrix
coa_matrix = cal_coancestry(n_ref_each=c(40,40),map=read.table('p.map')[,3:4],
                            method='Viterbi',fix_rho=TRUE,rate=1e-8)

# comparison -- significant difference, iteration time by default is 20
coa_matrix2 = cal_coancestry(n_ref_each=c(40,40),map=read.table('p.map')[,3:4],
                            method='EM',fix_rho=TRUE,rate=1e-8)

# if we greatly increase iteration times for EM, the coancestry matrix becomes closer
coa_matrix3 = cal_coancestry(n_ref_each=c(40,40),map=read.table('p.map')[,3:4],
                            method='EM',fix_rho=TRUE,rate=1e-8,ite_time=60)

#compute LDA and LDAS
lda=LDA(painting_all)
# gd should be in cM in function LDAS
map_use=data.frame(SNP=map[nrow(map):1,2],gd=gd[length(gd):1]*100) 
ldas = LDAS(lda,map_use,window=2)
plot(x=ldas$SNP/1e6,y=ldas$LDAS,xlab='Physical position (cM)',ylab='LDAS')
