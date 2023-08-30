sim_map <- read.table("sim_map.txt", header = TRUE)
sim_map$chromosome <- 1 
sim_map$snp <- paste0("SNP_", 1:nrow(sim_map)) 
sim_map <- sim_map[, c("chromosome", "snp", "gd", "pd")]
write.table(sim_map, "sim.map", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
