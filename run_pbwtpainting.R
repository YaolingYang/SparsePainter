source('painting_functions.R')

matchdata=read.table('matchdata.txt')[,c(1,3,5,6)]

map=read.table('p1new.map')[,4]

painting=cal_painting_full(matchdata,i=1,map,n_ref_each=c(40,40),
                           fix_rho=90.6536,returneachref=TRUE)
