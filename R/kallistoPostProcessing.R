library(sleuth)
library(tximport)
library(biomaRt)
library(rtracklayer)
library(GenomicFeatures)
library(readr)
library(data.table)
library(rentrez)
library(dplyr)
library(gtools)
library(edgeR)
# external functions ------------------------------------------------------
source("~/Development/GeneralPurpose/R/amILocal.R")
source("~/Development/GeneralPurpose/R/heatmap.3.R")
source("~/Development/GeneralPurpose/R/lsos.R")
source("~/Development/GeneralPurpose/R/loadConfigFile.R")

# load snakemake configuration --------------------------------------------
runConfig <- loadConfigFile("~/Development/JCSMR-Tremethick-Lab/HP1-alpha/snakemake/configs/config_RNA-Seq.json", file_type = "JSON")

# local functions ---------------------------------------------------------
lDir <- function(x, y){
  paste(x, y, sep = "/")
}

# global variables --------------------------------------------------------
colors <- RColorBrewer::brewer.pal(3, "Set2")
#refVersion <- runConfig$references$active
# manually set to hg19
refVersion <- "hg19"
dataset <- runConfig$references$database$dataset
biomart <- runConfig$references$database$biomart

#ensemblHost <- runConfig$references[[refVersion]]$ensemblHost
ensemblHost <- "http://grch37.ensembl.org"
annotationVersion <- runConfig$references[[refVersion]]$version
# reset to Ensembl75 ERCC
annotationVersion <- annotationVersion[2]
# manually set to "RNA-Seq"
assay <- names(runConfig$samples)[1]
processedDir <- runConfig$processed_dir


# environment setup -------------------------------------------------------
if (amILocal("JCSMR027564ML")){
  pathPrefix <- "~/mount/gduserv"
  mount <- system("mount", intern = T)
  if (length(grep("gduserv", mount)) == 0) {system("sshfs skurscheid@gduserv.anu.edu.au: ~/mount/gduserv/")}
  cpus <- 8
} else {
  pathPrefix <- "~"
  cpus <- 8
  options(width = 137)
}
options(mc.cores = cpus)


# change to WD -----------------------------------------------------------
devPath <- "~/Development"
runID <- runConfig$samples$`RNA-Seq`$runID
if (length(runID) == 1){
  wd <- lDir(pathPrefix, 
             paste("Data/Tremethick/HP1-alpha", 
                   assay, 
                   runID, 
                   processedDir,
                   annotationVersion, 
                   "R_Analysis",
                   sep = "/"))
  annotationDataPath <- lDir(pathPrefix, 
                             paste("Data/Tremethick/HP1-alpha", 
                                   assay, 
                                   runID, 
                                   processedDir,
                                   annotationVersion, 
                                   "R_Analysis", "",
                                   sep = "/"))
  base_dir <- paste(pathPrefix, 
                    "Data/Tremethick/HP1-alpha",
                    assay, 
                    runID, 
                    processedDir,
                    annotationVersion, 
                    "kallisto",
                    sep = "/")
  if (!dir.exists(wd)){
    dir.create(wd, recursive = T)
    setwd(wd)
  } else {
    setwd(wd)
  }
} else {
  wd <- lDir(pathPrefix, 
             paste("Data/Tremethick/HP1-alpha", 
                   assay, 
                   processedDir,
                   annotationVersion, 
                   "R_Analysis",
                   sep = "/"))
  annotationDataPath <- lDir(pathPrefix, 
                             paste("Data/Tremethick/HP1-alpha", 
                                   assay, 
                                   processedDir,
                                   annotationVersion, 
                                   "R_Analysis", "",
                                   sep = "/"))
  base_dir <-unlist(lapply(runID, function(x){
    paste(pathPrefix, 
          "Data/Tremethick/HP1-alpha",
          assay, 
          x, 
          processedDir,
          annotationVersion, 
          "kallisto",
          sep = "/")
  }))
  if (!dir.exists(wd)){
    dir.create(wd, recursive = T)
    setwd(wd)
  } else {
    setwd(wd)
  }
}

# file names for data output  ---------------------------------------------
analysis_version <- "1_wo_scrambled"
sleuth_results_output <- paste("sleuthResults_", annotationVersion, "_V", analysis_version, ".rda", sep = "")
sleuth_resultsCompressed_file <- paste("sleuthResultsCompressed_", annotationVersion, "_V", analysis_version, ".rda", sep = "")
clusterProfiler_results_output <- paste("clusterProfilerResults_", annotationVersion, "_V", analysis_version, ".rda", sep = "")

# read in additional parameters from JSON file ----------------------------
transcriptBiotype <- runConfig$samples$`RNA-Seq`$sleuth_parameters$transcript_biotype

# preparing annotation data from Ensembl ----------------------------------
# ToDo: refactor annotation data preparation 
# so that exsiting files will be re-used within the same project

# create biomaRt object
mart <- biomaRt::useEnsembl(biomart = biomart, dataset = dataset, host = ensemblHost)
attribs <- biomaRt::listAttributes(mart)

if (length(grep("ensembl", annotationVersion)) >= 1){
  ensGenes_file <- paste(annotationDataPath, "ensGenes_", annotationVersion, ".rda", sep = "")
  ensTranscripts_file <- paste(annotationDataPath, "ensTranscripts_", annotationVersion, ".rda", sep = "")
  t2g_file <- paste(annotationDataPath, "t2g_", annotationVersion, ".rda", sep = "")
  myLength_file <- paste(annotationDataPath, "mylength_", annotationVersion, ".rda", sep = "")
  myGC_file <- paste(annotationDataPath, "myGC_", annotationVersion, ".rda", sep = "")
  myBiotypes_file <- paste(annotationDataPath, "myBiotypes_", annotationVersion, ".rda", sep = "")
  myChroms_file <- paste(annotationDataPath, "myChroms_", annotationVersion, ".rda", sep = "")
  annotationFileList <- list(ensGenes_file, ensTranscripts_file, t2g_file, myLength_file, myGC_file, myBiotypes_file, myChroms_file) 
  if (!all(sapply(annotationFileList, file.exists))){
    ensGenes <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                              "external_gene_name",
                                              "chromosome_name",
                                              "start_position",
                                              "end_position",
                                              "strand",
                                              "band",
                                              "description",
                                              "percentage_gc_content", # renames in current ensembl release
                                              "gene_biotype",
                                              "entrezgene"),
                               mart = mart)
    ensGenes <- data.table::as.data.table(ensGenes)
    save(ensGenes, file = ensGenes_file)
    
    # get Ensembl transcripts
    ensTranscripts <- biomaRt::getBM(attributes = c("ensembl_transcript_id",
                                                    "ensembl_gene_id",
                                                    "transcript_length",
                                                    "version", 
                                                    "transcript_version",
                                                    "external_gene_name",
                                                    "description",
                                                    "transcript_biotype",
                                                    "transcript_version"),
                                     mart = mart,
                                     filter = "ensembl_gene_id",
                                     values = ensGenes$ensembl_gene_id)
    ensTranscripts <- data.table::data.table(ensTranscripts)
    save(ensTranscripts, file = ensTranscripts_file)
    # create t2g object
    t2g <- ensTranscripts[, c("ensembl_transcript_id", 
                              "ensembl_gene_id", 
                              "external_gene_name",
                              "transcript_biotype",
                              "description",
                              "transcript_version")]
    if(refVersion == "hg38"){
      t2g$ensembl_transcript_id <- paste(t2g$ensembl_transcript_id, t2g$transcript_version, sep = ".")
    }
    t2g <- data.table::data.table(t2g[, c("ensembl_transcript_id", 
                                          "ensembl_gene_id", 
                                          "external_gene_name",
                                          "transcript_biotype",
                                          "description")])
    t2g <- subset(t2g, transcript_biotype %in% transcriptBiotype)
    t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
    if(length(grep("ERCC", ensGenes_file)) > 0){
      erccGenes <- import("~/Data/References/Transcriptomes/ERCC/ERCC92.gtf")
      erccGenes <- data.table::data.table(target_id = erccGenes$gene_id, 
                                          ens_gene = erccGenes$gene_id, 
                                          ext_gene = erccGenes$gene_id,
                                          transcript_biotype = "errc_control")
      erccGenes$target_id <- paste(erccGenes$target_id, "1", sep = ".")
      t2g <- rbind(t2g, erccGenes)
    }
    save(t2g, file = t2g_file)
    
    mylength <- sapply(ensGenes$ensembl_gene_id, function(x){
      y <- ensTranscripts[which(ensTranscripts$ensembl_gene_id == x), ]
      y <- y[which.max(y$transcript_length), ]$transcript_length})
    save(mylength, file = myLength_file)
    mygc <- ensGenes$percentage_gc_content
    names(mygc) <- ensGenes$ensembl_gene_id
    save(mygc, file = myGC_file)
    mybiotypes <- ensGenes$gene_biotype
    names(mybiotypes) <- ensGenes$ensembl_gene_id
    save(mybiotypes, file = myBiotypes_file)
    mychroms <- data.frame(Chr = ensGenes$chromosome_name, GeneStart = ensGenes$start_position, GeneEnd = ensGenes$end_position)
    save(mychroms, file = myChroms_file)
    geneIdCol <- "ens_gene"
    targetIdCol <- "target_id"
  } else {
    load(ensGenes_file)
    load(ensTranscripts_file)
    load(myLength_file)
    load(myGC_file)
    load(myBiotypes_file)
    load(myChroms_file)
    load(t2g_file)
    geneIdCol <- "ext_gene"
    aggregation_column = "ext_gene"
    targetIdCol <- "target_id"
  }
} else if (length(grep("RefSeq", annotationVersion)) >= 1) {
  # load RefSeq NM IDs prior to accessing biomaRt in order to incorporate transcript versions
  t2g_file <- mRNArefSeq_IDs_file <- paste(unlist(strsplit(runConfig$references$hg38$kallisto$GRCh38_RefSeq_NM, "\\."))[1], "_t2g.", "rds", sep = "")
  if (file.exists(t2g_file)) {
    load(t2g_file)
  } else {
    mRNArefSeq_IDs_file <- paste(unlist(strsplit(runConfig$references$hg38$kallisto$GRCh38_RefSeq_NM, "\\."))[1], "txt", sep = ".")
    mRNArefSeq_IDs <- read_lines(mRNArefSeq_IDs_file)
    seq1 <- seq(1, length(mRNArefSeq_IDs), 200)
    library(snowfall)
    sfInit(parallel = T, cpus = 16)
    sfExport("seq1")
    sfExport("mRNArefSeq_IDs")
    sfLibrary(rentrez)
    sfLibrary(data.table)
    l1 <- sfLapply(seq_along(seq1), function(x) {
      if (x < length(seq1)){ 
        IDs <- mRNArefSeq_IDs[seq1[x] : seq1[x+1]]
      } else {
        IDs <-  mRNArefSeq_IDs[seq1[x] : length(mRNArefSeq_IDs)]
      }
      es <- entrez_summary(db = "nuccore", id = IDs)
      el <- entrez_link(dbfrom="nuccore", id = names(es), db = "gene", by_id = T)
      nuccore_gene <- unlist(lapply(1:length(el), function(x) el[[x]]$links$nuccore_gene))
      es1 <- entrez_summary(db="gene", id=nuccore_gene)
      dT <- data.table("refseq_mRNA" = unlist(lapply(es, function(x) x$accessionversion)),
                       "entrez_id" = names(es1), 
                       "gene_symbol" = unlist(lapply(es1, function(x) x$name)))
     return(dT)
    })
    t2g <- do.call(rbind, l1)
    t2g <- t2g[,c("refseq_mRNA", "gene_symbol")]
    data.table::setnames(t2g, "refseq_mRNA", "target_id")
    save(t2g, file = t2g_file)
  }
  geneIdCol = "gene_symbol"
  aggregation_column = "gene_symbol"
  targetIdCol = "target_id"
}


# load kallisto data with tximport and inspect via PCA -------------------------
if (length(base_dir) == 1){
  if (dir.exists(base_dir)){
    sample_id <- dir(base_dir)
    kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
    condition <- unlist(lapply(strsplit(sample_id, "_"), function(x) paste(x[1:2], collapse = "_")))
    names(condition) <- sample_id
    files <- paste(kal_dirs, "abundance.h5", sep = "/")
    names(files) <- sample_id
    } else { #end of directory IF
      stop("Directory with kallisto input data is missing!")
    }
  } else {
    kal_dirs <- unlist(lapply(runID, function(x) {
      paste(base_dir[grep(x, base_dir)], names(runConfig$samples$`RNA-Seq`[[x]]), sep = "/")
    }))
    sample_id <- unlist(lapply(runID, function(x) {
      names(runConfig$samples$`RNA-Seq`[[x]])
    }))
    condition <- unlist(lapply(strsplit(sample_id, "_"), function(x) paste(x[1:2], collapse = "_")))
    names(kal_dirs) <- sample_id
    names(condition) <- sample_id
    files <- paste(kal_dirs, "abundance.h5", sep = "/")
    names(files) <- sample_id
}

txi <- tximport::tximport(files, 
                          type = "kallisto",
                          tx2gene = subset(t2g, select = c(targetIdCol, geneIdCol)), geneIdCol = geneIdCol, txIdCol = targetIdCol,
                          ignoreTxVersion = F)

# perform PCA for first inspection of data --------------------------------
sd1 <- apply(log2(txi$abundance + 1), 1, sd)
summary(sd1)
pca1 <- ade4::dudi.pca(t(log2(txi$abundance[sd1 > 0, ] + 1)), scannf = F, nf = 6)
pdf(paste("PCA_MCF10A_HP1-alpha_", annotationVersion, ".pdf", sep = ""))
  ade4::s.arrow(pca1$li, boxes = F)
  ade4::s.class(pca1$li, fac = as.factor(condition))
dev.off()
# based on visual inspection:
# MCF10A_WT_[123] removed
# MCF10A_Scramble_2 removed 
# MCF10A_shHP1b_3 removed due to low knock-down efficiency 
# MCF10A_shHP1ab_3 removed due to low knock-down efficiency

toRemove <- c("MCF10A_WT_1|MCF10A_WT_2|MCF10A_WT_3|MCF10A_Scramble_")
files <- files[-grep(toRemove, names(files))]
length(files)
condition <- condition[-grep(toRemove, names(condition))]
length(condition)
sample_id <- sample_id[-grep(toRemove, sample_id)]
kal_dirs <- kal_dirs[-grep(toRemove, kal_dirs)]

txi <- tximport::tximport(files, 
                          type = "kallisto",
                          tx2gene = subset(t2g, select = c(targetIdCol, geneIdCol)), geneIdCol = geneIdCol, txIdCol = targetIdCol,
                          ignoreTxVersion = F)
sd1 <- apply(log2(txi$abundance + 1), 1, sd)
summary(sd1)
pca1 <- ade4::dudi.pca(t(log2(txi$abundance[sd1 > 0.25, ] + 1)), scannf = F, nf = 6)
pdf(paste("PCA_MCF10A_HP1-alpha_post_sample_removal_", annotationVersion, ".pdf", sep = ""))
  ade4::s.arrow(pca1$li, boxes = F)
  ade4::s.class(pca1$li, fac = as.factor(condition))
dev.off()


# prepare data.frame for importing kallisto data into sleuth --------------
s2c <- data.frame(sample = sample_id, condition = condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
s2c$sample <- as.character(s2c$sample)
s2c$condition <- relevel(s2c$condition, ref = "MCF10A_WT")
# build a list for all contrasts to test ----------------------------------
contrastsList <- runConfig$samples$`RNA-Seq`$sleuth_parameters$contrasts
conditionsList <- runConfig$samples$`RNA-Seq`$sleuth_parameters$conditions
design <- model.matrix(~ condition, data = s2c)


# run sleuth --------------------------------------------------------------
if(!file.exists(sleuth_results_output)){
    #----------------------------------------------------------------------
    # only doing gene level DE
    so.gene <- sleuth::sleuth_prep(s2c, 
                                   ~ condition, 
                                   target_mapping = t2g, 
                                   aggregation_column = aggregation_column)
    so.gene <- sleuth::sleuth_fit(so.gene, formula = design)
    so.gene <- sleuth::sleuth_fit(so.gene, ~1, "reduced")
    so.gene <- sleuth::sleuth_lrt(so.gene, "reduced", "full")
    for (i in colnames(design)[grep("Intercept", colnames(design), invert = T)]){
      so.gene <- sleuth::sleuth_wt(so.gene, i)  
    }
    rt.gene.list <- lapply(colnames(design)[grep("Intercept", colnames(design), invert = T)], function(x){
      rt.gene <- sleuth::sleuth_results(so.gene, x, show_all = F)
      rt.gene <- data.table::data.table(rt.gene[order(rt.gene$qval),])
    })
    names(rt.gene.list) <- colnames(design)[grep("Intercept", colnames(design), invert = T)]
    kt.gene <- data.table::data.table(sleuth::kallisto_table(so.gene, 
                                                             use_filtered = T, 
                                                             normalized = T, 
                                                             include_covariates = T))
    kt_wide.gene <- tidyr::spread(kt.gene[, c("target_id", "sample", "tpm")], 
                                  sample, 
                                  tpm)
    kt_wide.gene <- data.table::as.data.table(merge(kt_wide.gene, 
                                                    unique(subset(ensGenes, 
                                                                  select =  c("ensembl_gene_id", 
                                                                              "external_gene_name",
                                                                              "description")),
                                                                  by = "ensembl_gene_id"), 
                                                    by.x = "target_id", 
                                                    by.y = "external_gene_name",
                                                    all.x = TRUE, 
                                                    all.y = FALSE))
    kt_wide.gene <- kt_wide.gene[!duplicated(kt_wide.gene$target_id),]
    results <- list(sleuth_object_genes = so.gene,
                    sleuth_results_genes = rt.gene.list,
                    kallisto_table_genes = kt.gene,
                    kallisto_table_genes_wide = kt_wide.gene)
  save(results, file = sleuth_results_output)
} else {
  load(sleuth_results_output)
}

# export sleuth_results_genes as CSV --------------------------------------
lapply(names(results$sleuth_results_genes), function(x){
  tab1 <- data.table::as.data.table(results$sleuth_results_genes[[x]])
  tab1 <- data.table::as.data.table(merge(tab1, 
                                          unique(subset(ensGenes, 
                                                        select =  c("ensembl_gene_id", 
                                                                    "external_gene_name", 
                                                                    "description")), 
                                                 by = "ensembl_gene_id"), 
                                          by.x = "target_id", 
                                          by.y = "external_gene_name", 
                                          all.x = TRUE, 
                                          all.y = FALSE))
  tab1 <- tab1[order(tab1$qval, decreasing = F),]
  tab1 <- tab1[!duplicated(tab1$target_id),]
  # added filtering for genes
  tab1 <- subset(tab1, qval < 0.2 | abs(b) > 1) 
  fn <- paste(x, "_differential_gene_expression_results_table.csv", sep = "")
  write.csv(x = tab1, file = fn, row.names = FALSE)
})


# perform gene set enrichment analysis --------------------------------------------------
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

if (!file.exists(clusterProfiler_results_output)){
clusterProfilerResults <- sfLapply(results$sleuth_results_genes, function(x) {
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
  save(clusterProfilerResults, file = clusterProfiler_results_output)
} else {
  load(clusterProfiler_results_output)
}


# write enrichment analysis results to separate CSV files -----------------
lapply(names(clusterProfilerResults), function(x){
  sapply(names(clusterProfilerResults[[x]]$gmtList), function(y){
    df <- data.frame(clusterProfilerResults[[x]]$gmtList[[y]])
    if(nrow(df) > 0) {
      df$http_link <- paste("http://software.broadinstitute.org/gsea/msigdb/cards", df$ID, sep = "/")
      fileName <- paste(x, y, "csv", sep = ".")
      write.csv(df, file = fileName)
    }
  })
})


# print dotplots for GO enrichment ----------------------------------------
pdf(file = "Dotplots.pdf", paper = "a4r")
lapply(names(clusterProfilerResults), function(x){
  egos <- names(clusterProfilerResults[[x]])[grep("ego", names(clusterProfilerResults[[x]]))]
  obj <- clusterProfilerResults[[x]]
  lapply(egos, function(y){
    dotplot(obj[[y]], title = paste(x, y, sep ="_"))
    write.csv(data.frame(obj[[y]]), file = paste(x, y, "table.csv", sep = ""))
  })
})
dev.off()



