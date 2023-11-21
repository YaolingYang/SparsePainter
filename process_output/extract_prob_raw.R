library(data.table)
raw=fread("chr19_1000G_raw_10inds_prob.txt.gz",header = T,skip=1)
SNPs=as.numeric(colnames(raw)[-1])
nsnp=length(SNPs)
nind=nrow(raw)
indnames=raw[,1]

# each SNP is stored in a list
painting <- vector("list", length = nsnp)
for(i in 2:ncol(raw)){
  print(i)
  painting_temp=raw[[i]]
  painting_temp_list <- strsplit(painting_temp, split = " ")
  painting_temp_list <- lapply(painting_temp_list, function(x) as.numeric(unlist(strsplit(gsub("\\|", ",", x), split = ","))))
  painting_temp_list <- lapply(painting_temp_list, function(x) 0.5*(x[1:26]/sum(x[1:26])+x[27:52]/sum(x[27:52])))
  painting[[i-1]] <- do.call(rbind, painting_temp_list)
}

array_data <- array(unlist(painting), dim = c(nind, npop, nsnp))

painting_final <- lapply(1:npop, function(x) {
  array_data[, x, ]
})
