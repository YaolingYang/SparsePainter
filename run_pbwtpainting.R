source('painting_functions.R')

p1=read.table('p1match.txt')[,c(1,3,5,6)]
p2=read.table('p2match.txt')[,c(1,3,5,6)]
p2[,1]=p2[,1]+40
matchdata=rbind(p1,p2)
## above should be done in bash line later

map=read.table('p1new.map')[,4]

painting=cal_painting_full(matchdata,i=1,map,n_ref_each=c(40,40),
                           fix_rho=90.6536,returneachref=FALSE)
