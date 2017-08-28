library(clusterProfiler)
library(org.Hs.eg.db)
gmtfiles <- list.files("~/Data/References/Annotations/MSigDB/msigdb_v6.0_GMTs", pattern = "symbols", full.names = T)
gmtfiles <- gmtfiles[grep("all", gmtfiles)]
gmtfilesNames <- unlist(lapply(strsplit(gmtfiles, "/"), function(x) x[9]))

library(snowfall)
sfInit(parallel = T, cpus = 5)
sfExport(list = c("results", "gmtfiles", "gmtfilesNames"))
sfLibrary(clusterProfiler)
sfLibrary(org.Hs.eg.db)

EnrichmentData <- sfLapply(results$sleuth_results_genes, function(x) {
  geneList <- x$b
  names(geneList) <- x$target_id
  geneList <- sort(geneList, decreasing = T)
  gene <- names(geneList)[abs(geneList) > 1]
  # for the GO analysis we include only differentially expressed genes
  # effect size is > |1|, but no p/q-value cutoff
  egoMF <- enrichGO(gene         = gene,
                    universe = names(geneList),
                    OrgDb         = org.Hs.eg.db,
                    keytype       = 'SYMBOL',
                    ont           = "MF",
                    pAdjustMethod = "fdr",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
  egoBP <- enrichGO(gene         = gene,
                    universe = names(geneList),
                    OrgDb         = org.Hs.eg.db,
                    keytype       = 'SYMBOL',
                    ont           = "BP",
                    pAdjustMethod = "fdr",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
  egoCC <- enrichGO(gene         = gene,
                    universe = names(geneList),
                    OrgDb         = org.Hs.eg.db,
                    keytype       = 'SYMBOL',
                    ont           = "CC",
                    pAdjustMethod = "fdr",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
  
  # for the MSigDB analysis we use the whole list of ranked genes
  gmtList <- lapply(gmtfiles, function(y){
    gmt <- read.gmt(y)
    egmt <- GSEA(geneList = geneList, 
                 TERM2GENE=gmt, 
                 nPerm = 1000, 
                 exponent = 1, 
                 seed = T, 
                 pAdjustMethod = "fdr")
  })
  names(gmtList) <- gmtfilesNames
  return(list(egoBP = egoBP, 
              egoMF = egoMF, 
              egoCC = egoCC,
              gmtList = gmtList))
})

lapply(names(EnrichmentData), function(x){
  sapply(names(EnrichmentData[[x]]$gmtList), function(y){
    df <- data.frame(EnrichmentData[[x]]$gmtList[[y]])
    if(nrow(df) > 0) {
      fileName <- paste(x, y, ".csv", sep = ".")
      write.csv(df, file = fileName)
    }
  })
})



