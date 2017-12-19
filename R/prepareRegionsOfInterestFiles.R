# prepare BED files of DE genes
# double KD
library(deepToolsUtils)
library(GenomicRanges)

annotationVersion <- "GRCh37_hg19_ensembl75"
if (!dir.exists(paste("~/Data/Tremethick/HP1-alpha/AnnotationData", annotationVersion, sep = "/"))) {
  dir.create(paste("~/Data/Tremethick/HP1-alpha/AnnotationData", annotationVersion, sep =  "/"))
  pathPrefix <- paste("~/Data/Tremethick/HP1-alpha/AnnotationData", annotationVersion, sep = "/")
} else {
  pathPrefix <- paste("~/Data/Tremethick/HP1-alpha/AnnotationData", annotationVersion, sep = "/")
}

lapply(names(results$sleuth_results_genes), function(x){
  id <- results$sleuth_results_genes[[x]][qval < 0.1]$target_id
  dT1 <- ensGenes[external_gene_name %in% id]
  gr1 <- GenomicRanges::GRanges(seqnames = dT1$chromosome_name,
                                IRanges(start = dT1$start_position,
                                        end = dT1$end_position,
                                        names = dT1$external_gene_name),
                                strand = c("+", "-")[match(dT1$strand, c(1, -1))],
                                category = rep(x, nrow(dT1)))
  out_file <- paste(pathPrefix, "/", x, ".bed", sep = "")
  deepToolsUtils::WriteGRangesToBED(gr = gr1, out_file = out_file)
})

lapply(names(results$sleuth_results_genes), function(x){
  id <- results$sleuth_results_genes[[x]][qval < 0.1]$target_id
  dT1 <- ensGenes[external_gene_name %in% id]
  gr1 <- GenomicRanges::GRanges(seqnames = dT1$chromosome_name,
                                IRanges(start = dT1$start_position,
                                        end = dT1$end_position,
                                        names = dT1$external_gene_name),
                                strand = c("+", "-")[match(dT1$strand, c(1, -1))],
                                category = rep(x, nrow(dT1)))
  strand(gr1) <- "*"
  out_file <- paste(pathPrefix, "/", x, "_unstranded.bed", sep = "")
  deepToolsUtils::WriteGRangesToBED(gr = gr1, out_file = out_file)
})

dT1 <- ensGenes

# makes GR from ensGenes
gr1 <-  GenomicRanges::GRanges(seqnames = dT1$chromosome_name,
                                 IRanges(start = dT1$start_position,
                                         end = dT1$end_position,
                                         names = dT1$external_gene_name),
                                 strand = c("+", "-")[match(dT1$strand, c(1, -1))],
                                 category = rep("allGenes", nrow(dT1)))
strand(gr1) <- "*"
out_file <- paste(pathPrefix, "/", "allGenes", ".bed", sep = "")
deepToolsUtils::WriteGRangesToBED(gr = gr1, out_file = out_file)

library("rtracklayer")
rm1 <- import("~/Data/References/Annotations/Homo_sapiens/hg19/UCSC/hg19_repeat_regions.bed")
rm1
# map Ensembl to UCSC chromosome IDs
chromMap <- data.table::fread("/home/sebastian/Data/References/Annotations/ChromosomeMappings/GRCh37_UCSC2ensembl.txt", sep = "\t")
seqlevels(rm1) <- chromMap$V2[match(seqlevels(rm1), chromMap$V1)]

deepToolsUtils::WriteGRangesToBED(gr = rm1[grep("GSAT", rm1$name)], out_file = paste(pathPrefix, "hg19_GSAT_repeats.bed", sep = "/"))
deepToolsUtils::WriteGRangesToBED(gr = rm1[grep("LTR", rm1$name)], out_file = paste(pathPrefix, "hg19_LTR_repeats.bed", sep = "/"))

repeatMaskerTab <- data.table::fread("/home/sebastian/Data/References/Annotations/Homo_sapiens/hg19/UCSC/repeatMaskerhg19.txt", sep = "\t")
setkey(repeatMaskerTab, repClass)
for (i in c("LINE", "SINE", "Low_complexity")){ 
  gr1 <- GenomicRanges::GRanges(seqnames = repeatMaskerTab[i]$genoName,
                                IRanges(start = repeatMaskerTab[i]$genoStart,
                                        end = repeatMaskerTab[i]$genoEnd,
                                        names = repeatMaskerTab[i]$repName),
                                strand = repeatMaskerTab[i]$strand,
                                repClass = repeatMaskerTab[i]$repClass)
  seqlevels(gr1) <- chromMap$V2[match(seqlevels(gr1), chromMap$V1)]
  fn <- paste(pathPrefix, "/hg19_", i, "_repeats.bed", sep = "")
  deepToolsUtils::WriteGRangesToBED(gr = gr1, out_file = fn)
}





