# prepare BED files of DE genes
# double KD
library(deepToolsUtils)
library(GenomicRanges)

dT1 <- ensGenes[external_gene_name %in% results$sleuth_results_genes$conditionMCF10A_shHP1ab[qval < 0.1]$target_id]
gr1 <- GenomicRanges::GRanges(seqnames = dT1$chromosome_name,
                              IRanges(start = dT1$start_position,
                                      end = dT1$end_position,
                                      names = dT1$external_gene_name),
                              strand = c("+", "-")[match(dT1$strand, c(1, -1))],
                              category = rep("shHP1ab", nrow(dT1)))
deepToolsUtils::WriteGRangesToBED(gr = gr1, out_file = "test.bed")
