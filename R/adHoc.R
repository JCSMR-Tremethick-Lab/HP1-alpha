toCheck <- c("DGKI", "CADM3", "PADI1", "SCG5")
results$kallisto_table_genes[target_id %in% toCheck]
p1 <- ggplot(data = results$kallisto_table_genes[target_id %in% toCheck & condition %in% c("MCF10A_WT", "MCF10A_shH2AZHP1a")], aes(condition,tpm)) + geom_jitter(width = 0.1, aes(colour = sample)) + facet_wrap(~target_id, scales = "free_y")
