setwd("C:/Ubuntu/longmatchquery")
library("Rcpp")
sourceCpp(file="hashmap_Rcpp.cpp")

sourceCpp(file="hashmap_Rcppsugar.cpp")
