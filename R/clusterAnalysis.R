library(vegan)
library(ggplot2)
library(data.table)
# prepare bed files of genes which are covered in RNA-Seq data ----------
setkey(ensGenes, "ensembl_gene_id")
dT1 <- ensGenes[results$kallisto_table_genes_wide$ensembl_gene_id]
dT1 <- dT1[!duplicated(dT1$ensembl_gene_id)]
gr1 <-  GenomicRanges::GRanges(seqnames = dT1$chromosome_name,
                               IRanges(start = dT1$start_position,
                                       end = dT1$end_position,
                                       names = dT1$external_gene_name),
                               strand = c("+", "-")[match(dT1$strand, c(1, -1))],
                               category = rep("allGenesRNASeq", nrow(dT1)))
# strand(gr1) <- "*"
out_file <- paste(pathPrefix, "/", "allGenesRNASeq", ".bed", sep = "")
deepToolsUtils::WriteGRangesToBED(gr = gr1, out_file = out_file)

# genes ranked by WT expression
kT <- results$kallisto_table_genes_wide
dT <- data.table("ext_gene" = kT$target_id, subset(kT, select=(grep("MCF10A", colnames(kT)))), "cluster" = "NA")
dTm <- melt.data.table(dT)
dTm$value <- log2(dTm$value + 1)
dTm$group <- unlist(lapply(strsplit(as.character(dTm$variable), "_"), function(x) x[2]))
dTmeanExpression <- dTm[, j=mean(value), by=list(group, ext_gene)]
setkey(dTmeanExpression, "group")
dTmeanExpression["WT"][order(-V1)]
rankedWTgenes <- dTmeanExpression["WT"][order(-V1)]$ext_gene
setkey(ensGenes, "external_gene_name")
dT1 <- ensGenes[rankedWTgenes]
gr1 <-  GenomicRanges::GRanges(seqnames = dT1$chromosome_name,
                               IRanges(start = dT1$start_position,
                                       end = dT1$end_position,
                                       names = dT1$external_gene_name),
                               strand = c("+", "-")[match(dT1$strand, c(1, -1))],
                               category = rep("allGenesWTExpressionRanked", nrow(dT1)))

out_file <- paste(pathPrefix, "/", "allGenesWTExpressionRanked", ".bed", sep = "")
deepToolsUtils::WriteGRangesToBED(gr = gr1[rankedWTgenes], out_file = out_file)

# use k-means to cluster expression data ----------------------------------
kT <- results$kallisto_table_genes_wide
dT <- data.table("ext_gene" = kT$target_id, subset(kT, select=(grep("MCF10A", colnames(kT)))), "cluster" = "NA")
mat1 <- as.matrix(subset(kT, select=(grep("MCF10A", colnames(kT)))))
mat1 <- log2(mat1 + 1)
rownames(mat1) <-  kT$target_id

# run pca for scaling and centering
pca1 <- ade4::dudi.pca(t(mat1), scannf = F, nf = 5)
screeplot(pca1)
head(t(pca1$tab))
adegraphics::s.arrow(pca1$li)
s.class(pca1$li, fac = condition)

mat2 <- t(pca1$tab)
ccKM1 <- cascadeKM(mat2, 2, 10, 100, criterion = "ssi")
plot(ccKM1)
group7 <- rownames(ccKM1$partition[which(ccKM1$partition[,6] == 7),])
group6 <- rownames(ccKM1$partition[which(ccKM1$partition[,6] == 6),])
group5 <- rownames(ccKM1$partition[which(ccKM1$partition[,6] == 5),])
group4 <- rownames(ccKM1$partition[which(ccKM1$partition[,6] == 4),])
group3 <- rownames(ccKM1$partition[which(ccKM1$partition[,6] == 3),])
group2 <- rownames(ccKM1$partition[which(ccKM1$partition[,6] == 2),])
group1 <- rownames(ccKM1$partition[which(ccKM1$partition[,6] == 1),])


setkey(dT, "ext_gene")
dT[group1, "cluster"] <- "cluster1"
dT[group2, "cluster"] <- "cluster2"
dT[group3, "cluster"] <- "cluster3"
dT[group4, "cluster"] <- "cluster4"
dT[group5, "cluster"] <- "cluster5"
dT[group6, "cluster"] <- "cluster6"
dT[group7, "cluster"] <- "cluster7"
table(dT$cluster)
dTm <- melt.data.table(dT)
dTm$value <- log2(dTm$value + 1)
dTm$group <- unlist(lapply(strsplit(as.character(dTm$variable), "_"), function(x) x[2]))
bp1 <- ggplot(data = dTm, aes(x = cluster, y = value)) + geom_boxplot()
bp1
bp2 <- ggplot(data = dTm, aes(x = group, y = value)) + geom_boxplot() + facet_wrap(~cluster)
bp2
vp1 <- ggplot(data = dTm, aes(x = cluster, y = value)) + geom_violin() + facet_wrap(~group)
vp1

setkey(dTm, "ext_gene")
dTm[group1, j=mean(value), by=group]
dTm[group2, j=mean(value), by=group]
dTm[group3, j=mean(value), by=group]
dTm[group4, j=mean(value), by=group]
dTm[group5, j=mean(value), by=group]
dTm[group6, j=mean(value), by=group]
dTm[group7, j=mean(value), by=group]
dTm[group1, j=median(value), by=group]
dTm[group2, j=median(value), by=group]
dTm[group3, j=median(value), by=group]
dTm[group4, j=median(value), by=group]
dTm[group5, j=median(value), by=group]
dTm[group6, j=median(value), by=group]
dTm[group7, j=median(value), by=group]


# do k-means clustering on single contrasts -------------------------------
selCol <- grep("Scramble|shHP1ab_", colnames(kT))
dT <- data.table(subset(kT, select=c(1, selCol)), "cluster" = "NA")
mat1 <- as.matrix(subset(kT, select=(selCol)))
mat1 <- log2(mat1 + 1)
rownames(mat1) <-  kT$ext_gene
ccKM1 <- cascadeKM(mat1, 2, 10, 100, criterion = "calinski")
plot(ccKM1)
group3 <- rownames(ccKM1$partition[which(ccKM1$partition[,2] == 3),])
group2 <- rownames(ccKM1$partition[which(ccKM1$partition[,2] == 2),])
group1 <- rownames(ccKM1$partition[which(ccKM1$partition[,2] == 1),])
setkey(dT, "ext_gene")
dT[group1, "cluster"] <- "cluster1"
dT[group2, "cluster"] <- "cluster2"
dT[group3, "cluster"] <- "cluster3"
table(dT$cluster)
dTm <- melt.data.table(dT)
dTm$value <- log2(dTm$value + 1)
dTm$group <- unlist(lapply(strsplit(as.character(dTm$variable), "_"), function(x) x[2]))
setkey(dTm, "ext_gene")
dTm[sample(group1, size = 20, replace = F), j = mean(value), by = list(group, ext_gene)]
# get ChIP based k-means clustering info -----------------------------------
ChIPClusters <- data.table::fread("/home/sebastian/Data/Tremethick/HP1-alpha/Figures/bigwigCompare.allGenesRNASeq.normal.RPKM.kmeans6.bed")
table(ChIPClusters$deepTools_group)
cluster_1 <- ChIPClusters[deepTools_group == "cluster_1"]$name
cluster_2 <- ChIPClusters[deepTools_group == "cluster_2"]$name
cluster_3 <- ChIPClusters[deepTools_group == "cluster_3"]$name
cluster_4 <- ChIPClusters[deepTools_group == "cluster_4"]$name
cluster_5 <- ChIPClusters[deepTools_group == "cluster_5"]$name
cluster_6 <- ChIPClusters[deepTools_group == "cluster_6"]$name

kT <- results$kallisto_table_genes_wide
dT <- data.table("ext_gene" = kT$target_id, subset(kT, select=(grep("MCF10A", colnames(kT)))), "cluster" = "NA")
dT <- dT[ext_gene %in% ChIPClusters$name]
setkey(dT, "ext_gene")
dT[cluster_1, "cluster"] <- "cluster_1"
dT[cluster_2, "cluster"] <- "cluster_2"
dT[cluster_3, "cluster"] <- "cluster_3"
dT[cluster_4, "cluster"] <- "cluster_4"
dT[cluster_5, "cluster"] <- "cluster_5"
dT[cluster_6, "cluster"] <- "cluster_6"
table(dT$cluster,useNA = "ifany")
dTm <- melt.data.table(dT)
dTm$value <- log2(dTm$value + 1)
dTm$group <- unlist(lapply(strsplit(as.character(dTm$variable), "_"), function(x) x[2]))
dTm$group <- (factor(dTm$group, levels = c("WT", "shHP1a", "shHP1b", "shHP1ab", "shH2AZHP1a")))
library(ggplot2)
bp1 <- ggplot(data = dTm, aes(x = cluster, y = value)) + geom_boxplot() + facet_wrap(~group)
bp1
vp1 <- ggplot(data = dTm, aes(x = cluster, y = value)) + geom_violin() + facet_wrap(~group)
vp1
setkey(dTm, "ext_gene")
dTm[cluster_1, j=mean(value), by=group]
dTm[cluster_2, j=mean(value), by=group]
dTm[cluster_3, j=mean(value), by=group]
dTm[cluster_4, j=mean(value), by=group]
dTm[cluster_5, j=mean(value), by=group]
dTm[cluster_6, j=median(value), by=group]
dTm[cluster_1, j=median(value), by=group]
dTm[cluster_2, j=median(value), by=group]
dTm[cluster_3, j=median(value), by=group]
dTm[cluster_4, j=median(value), by=group]
dTm[cluster_5, j=median(value), by=group]
dTm[cluster_6, j=median(value), by=group]

vp2 <- ggplot(data = dTm, aes(x = group, y = value)) + geom_violin() + facet_wrap(~cluster) + theme(axis.text.x = element_text(angle = 45))
vp2 +geom_boxplot(width=.1)
bp2 <- ggplot(data = dTm, aes(x = group, y = value)) + geom_boxplot() + facet_wrap(~cluster) + theme(axis.text.x = element_text(angle = 45))
bp2



pca2 <- dudi.pca(t(log2(txi$abundance[sd1 > 1,] + 1)), scannf = F, nf = 6)
scatter(pca2, density.plot = T, box = F)

ade4::s.arrow(pca2$li, boxes = F)
ade4::s.class(pca2$li, fac = as.factor(condition))
dim(pca2$tab)


# organising data into quartiles, based on WT -----------------------------
l1 <- lapply(unique(dTmeanExpression$group), function(x){
  dTsort <- dTmeanExpression[x][order(-V1)]
  quartiles <- quantile(dTsort$V1)
  sd1 <- sd(dTsort$V1)
  dTsort$quartile <- "NA"
  dTsort[V1 <= quartiles[2]]$quartile <- "Q1"
  dTsort[V1 > quartiles[2] & V1 <= quartiles[3]]$quartile <- "Q2"
  dTsort[V1 > quartiles[3] & V1 <= quartiles[4]]$quartile <- "Q3"
  dTsort[V1 > quartiles[4]]$quartile <- "Q4"
  return(list(quartiles = quartiles, sd = sd1, dTsort = dTsort))
})
names(l1) <- unique(dTmeanExpression$group)

dTq <- l1$WT$dTsort
colnames(dTq)[4] <- "WT_quartile"
dTq <- merge(dTq, subset(l1$shHP1a$dTsort, select = c("ext_gene", "quartile")))
colnames(dTq)[5] <- "shHP1a_quartile"
dTq <- merge(dTq, subset(l1$shHP1b$dTsort, select = c("ext_gene", "quartile")))
colnames(dTq)[6] <- "shHP1b_quartile"
dTq <- merge(dTq, subset(l1$shHP1ab$dTsort, select = c("ext_gene", "quartile")))
colnames(dTq)[7] <- "shHP1ab_quartile"
dTq <- merge(dTq, subset(l1$shH2AZHP1a$dTsort, select = c("ext_gene", "quartile")))
colnames(dTq)[8] <- "shH2AZHP1a_quartile"
table(dTq[,c(4:8)])
library(GGally)
ggparcoord(data = dTq, columns = c(4:8))
library(ggparallel)
ggparallel(list("WT_quartile", "shHP1a_quartile", "shHP1b_quartile", "shHP1ab_quartile", "shH2AZHP1a_quartile"), data=dTq)
ggparallel(list("WT_quartile", "shHP1a_quartile"), data=dTq, method="hammock", ratio=0.25)
ggparallel(list("WT_quartile", "shH2AZHP1a_quartile"), data=dTq, method="hammock", ratio=0.25)
ggparallel(list("WT_quartile", "shHP1ab_quartile"), data=dTq, method="hammock", ratio=0.25)

dTmeanExpression$group <- (factor(dTmeanExpression$group, levels = c("WT", "shHP1a", "shHP1b", "shHP1ab", "shH2AZHP1a")))
setkey(dTmeanExpression, "ext_gene")
ggplot(dTmeanExpression[dTq[WT_quartile == "Q2" & shHP1ab_quartile == "Q4"]$ext_gene], aes(x = group, y = V1)) + geom_boxplot()

dTmq <- merge(dTm, subset(dTsort, select = c("ext_gene", "quartile")), by.x = "ext_gene", by.y = "ext_gene", all.x = T)
dTmq$group <- (factor(dTmq$group, levels = c("WT", "shHP1a", "shHP1b", "shHP1ab", "shH2AZHP1a")))
ggplot(dTmq, aes(x = group, y = value)) + geom_violin() + facet_wrap(~quartile) + geom_boxplot(width = .1)
ggplot(dTmq, aes(x = group, y = value)) + geom_violin()

setkey(dTsort, "quartile")
lapply(as.character(unique(dTsort$quartile)), function (x) {
  out_file <- paste(pathPrefix, "/", "allGenesWTExpression", x,  ".bed", sep = "")
  g1 <- dTsort[x]$ext_gene
  dT1 <- ensGenes[g1]
  gr1 <-  GenomicRanges::GRanges(seqnames = dT1$chromosome_name,
                                 IRanges(start = dT1$start_position,
                                         end = dT1$end_position,
                                         names = dT1$external_gene_name),
                                 strand = c("+", "-")[match(dT1$strand, c(1, -1))],
                                 category = rep("allGenesQuartiles", nrow(dT1)))
  deepToolsUtils::WriteGRangesToBED(gr = gr1, out_file = out_file)
})

# clusterProfiler analysis of genes changing quartiles
library(clusterProfiler)
library(org.Hs.eg.db)

gene <- c(dTq[which(dTq$WT_quartile == "Q3" & (dTq$shHP1a_quartile == "Q4" & dTq$shHP1b_quartile == "Q4" & dTq$shHP1ab_quartile == "Q4"))]$ext_gene,
          dTq[which(dTq$WT_quartile == "Q2" & (dTq$shHP1a_quartile == "Q3" & dTq$shHP1b_quartile == "Q3" & dTq$shHP1ab_quartile == "Q3"))]$ext_gene,
          dTq[which(dTq$WT_quartile == "Q3" & (dTq$shHP1a_quartile == "Q2" & dTq$shHP1b_quartile == "Q2" & dTq$shHP1ab_quartile == "Q2"))]$ext_gene,
          dTq[which(dTq$WT_quartile == "Q4" & (dTq$shHP1a_quartile == "Q3" & dTq$shHP1b_quartile == "Q3" & dTq$shHP1ab_quartile == "Q3"))]$ext_gene)

gene.df <- bitr(gene, fromType = "SYMBOL", 
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
universe <- bitr(dTq$ext_gene, fromType = "SYMBOL", 
                 toType = c("ENSEMBL", "ENTREZID"),
                 OrgDb = org.Hs.eg.db)

pvalueCutoff <- 0.01

egoMF <- enrichGO(gene         = gene.df$ENTREZID,
                  OrgDb         = 'org.Hs.eg.db',
                  ont           = "MF",
                  pvalueCutoff  = pvalueCutoff)
head(egoMF)

egoBP <- enrichGO(gene         = gene.df$ENTREZID,
                  OrgDb         = "org.Hs.eg.db",
                  ont           = "BP",
                  pvalueCutoff  = pvalueCutoff)
egoBP

egoCC <- enrichGO(gene         = gene.df$ENTREZID,
                  OrgDb         = "org.Hs.eg.db",
                  ont           = "CC",
                  pvalueCutoff  = pvalueCutoff)
egoCC
dotplot(egoBP)
dotplot(egoCC)


# genes rankes by WT expression -------------------------------------------
out_file
dT1 <- ensGenes[rankedWTgenes]
gr1 <-  GenomicRanges::GRanges(seqnames = dT1$chromosome_name,
                               IRanges(start = dT1$start_position,
                                       end = dT1$end_position,
                                       names = dT1$external_gene_name),
                               strand = c("+", "-")[match(dT1$strand, c(1, -1))],
                               category = rep("allGenesWTExpressionRanked", nrow(dT1)))
strand(gr1) <- "*"
out_file <- paste(pathPrefix, "/", "allGenesWTExpressionRanked", ".bed", sep = "")
deepToolsUtils::WriteGRangesToBED(gr = gr1[rankedWTgenes], out_file = out_file)



# look at exons
exonClusters <- data.table::fread("/home/sebastian/Data/Tremethick/HP1-alpha/Figures/bigwigCompare.expressedExons.normal.RPKM.kmeans4.bed")

#
