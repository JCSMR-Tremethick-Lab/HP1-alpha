# prepare BED files of DE genes
# double KD
library(deepToolsUtils)
library(GenomicRanges)
library(snowfall)

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


# use TxDB for extracting gene structures ---------------------------------
require(GenomicFeatures)
require(RSQLite)

TxDBFile <- "~/Data/References/Annotations/Homo_sapiens/GRCh37_hg19_ensembl75/hsapiens_gene_ensembl_GRCh37_TxDB.sqlite"
file.exists(TxDBFile)
hsapEnsemblTxDB <- loadDb(TxDBFile)
canonicalChr <- c(seq(1,22,1), "X", "MT")
seqlevels(hsapEnsemblTxDB) <- canonicalChr

# get exons ---------------------------------------------------------------
exons <- GenomicFeatures::exonsBy(hsapEnsemblTxDB, by = "gene")
exons <- unlist(exons)
keep <- seqlevels(exons) %in% unique(ChIPClusters$`#chrom`)
seqlevels(exons, pruning.mode = "tidy") <- seqlevels(exons)[keep]
exons <- sort(exons)
exons <- reduce(exons)

expressedGenes <- results$kallisto_table_genes_wide$ensembl_gene_id
exons[expressedGenes]
strand(exons) <- "*"
deepToolsUtils::WriteGRangesToBED(reduce(exons[names(exons) %in% expressedGenes]), out_file = "~/Data/Tremethick/HP1-alpha/AnnotationData/expressedExons.bed")


# get introns -------------------------------------------------------------
canonicalChr <- c(seq(1,22,1), "X", "MT")

introns <- intronicParts(hsapEnsemblTxDB, linked.to.single.gene.only = T)
grl.introns <- intronsByTranscript(hsapEnsemblTxDB, use.name = T)
sfInit(parallel = T, cpus = 16)
sfExport("grl.introns")
n1 <- sfLapply(grl.introns, function(x) {length(x)})
grl.introns <- grl.introns[names(n1[which(n1 > 0)])]

sfExport("grl.introns")
sfLibrary(biovizBase)
sfLibrary(GenomicRanges)
grl <- sfLapply(canonicalChr, function(x) biovizBase::flatGrl(grl.introns[which(GenomicRanges::seqnames(grl.introns) == x)]))
grl <- GenomicRanges::GRangesList(grl)
grl <- GenomicRanges::sort(grl)
grl <- GenomicRanges::reduce(grl)
grl.introns <- grl

sfExport("grl.introns")
grl.exonIntron <- GRangesList(sfLapply(grl.introns, function(x) {
  if (length(x) > 0){ 
    d <- data.frame(chr = seqnames(x), start = start(x) - 49, width = 100, strand = as.character(strand(x)))
    gr <- GRanges(as.character(d[,"chr"]), IRanges(start = as.numeric(d[,"start"]), width = as.numeric(d[,"width"])), strand = as.character(d[,"strand"]))
    return(gr)
  }
})
)

sfExport("grl.exonIntron")
grl <- sfLapply(canonicalChr, function(x) flatGrl(grl.exonIntron[which(seqnames(grl.exonIntron) == x)]))
grl <- GRangesList(grl)
grl <- sort(grl)
grl <- reduce(grl)
grl.exonIntron <- grl

grl.intronExon <- GRangesList(sfLapply(grl.introns, function(x) {
  if (length(x) > 0){ 
    d <- data.frame(chr = seqnames(x), start = end(x) - 49, width = 100, strand = as.character(strand(x)))
    gr <- GRanges(as.character(d[,"chr"]), IRanges(start = as.numeric(d[,"start"]), width = as.numeric(d[,"width"])), strand = as.character(d[,"strand"]))
    return(gr)
  }
}
)
)

sfExport("grl.intronExon")
grl <- sfLapply(canonicalChr, function(x) flatGrl(grl.intronExon[which(seqnames(grl.intronExon) == x)]))
grl <- GRangesList(grl)
grl <- sort(grl)
grl <- reduce(grl)
grl.intronExon <- grl

grl.intronExonFlank25 <- flank(grl.intronExon, width = 25, both = TRUE)
grl.exonIntronFlank25 <- flank(grl.exonIntron, width = 25, both = TRUE)
grl.intronsFlank25 <- flank(grl.introns, width = 25, both = TRUE)

gr.introns <- introns
seqlevels(gr.introns, pruning.mode = "tidy") <- canonicalChr
deepToolsUtils::WriteGRangesToBED(reduce(introns), 
                                  out_file = "~/Data/Tremethick/HP1-alpha/AnnotationData/GRCh37_hg19_ensembl75/allIntrons.bed")

