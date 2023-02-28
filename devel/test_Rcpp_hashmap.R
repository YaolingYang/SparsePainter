setwd("C:/Ubuntu/longmatchquery")
library("Rcpp")
sourceCpp(file="test_Rcpp_hashmap.cpp")
hash_test()
hash_test2()
hash_test3()
hash_test4()

