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
refVersion <- runConfig$references$active
dataset <- runConfig$references$database$dataset
biomart <- runConfig$references$database$biomart
ensemblHost <- runConfig$references$hg38$ensemblHost
annotationVersion <- runConfig$references[[refVersion]]$version
# reset to Ensembl84 ERCC
annotationVersion <- annotationVersion[1]
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
analysis_version <- "2"
sleuth_results_output <- paste("sleuthResults_", annotationVersion, "_V", analysis_version, ".rda", sep = "")
sleuth_resultsCompressed_file <- paste("sleuthResultsCompressed_", annotationVersion, "_V", analysis_version, ".rda", sep = "")

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
                                                    "external_gene_name"),
                                     mart = mart,
                                     filter = "ensembl_gene_id",
                                     values = ensGenes$ensembl_gene_id)
    save(ensTranscripts, file = ensTranscripts_file)
    # create t2g object
    t2g <- ensTranscripts[, c("ensembl_transcript_id", 
                              "ensembl_gene_id", 
                              "external_gene_name", 
                              "version", 
                              "transcript_version")]
    if(refVersion == "hg38"){
      t2g$ensembl_transcript_id <- paste(t2g$ensembl_transcript_id, t2g$transcript_version, sep = ".")
    }
    t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
    if(length(grep("ERCC", ensGenes_file)) > 0){
      erccGenes <- import("~/Data/References/Transcriptomes/ERCC/ERCC92.gtf")
      erccGenes <- data.frame(target_id = erccGenes$gene_id, 
                              ens_gene = erccGenes$gene_id, 
                              ext_gene = erccGenes$gene_id)
      erccGenes$version <- 1
      erccGenes$transcript_version <- 1
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
    geneIdCol <- "ext_gene"
    targetIdCol <- "target_id"
  } else {
    load(ensGenes_file)
    load(ensTranscripts_file)
    load(myLength_file)
    load(myGC_file)
    load(myBiotypes_file)
    load(myChroms_file)
    load(t2g_file)
    geneIdCol <- "external_gene_name"
    aggregation_column = "external_gene_name"
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
                          tx2gene = t2g, geneIdCol = geneIdCol, txIdCol = targetIdCol,
                          ignoreTxVersion = F)

# perform PCA for first inspection of data --------------------------------
sd1 <- apply(txi$abundance, 1, sd)
summary(sd1)
pca1 <- ade4::dudi.pca(t(txi$abundance[sd1 > 0, ]), scannf = F, nf = 6)
pdf(paste("PCA_MCF10A_HP1-alpha_", annotationVersion, ".pdf", sep = ""))?
  ade4::s.arrow(pca1$li, boxes = F)
  ade4::s.class(pca1$li, fac = as.factor(condition))
dev.off()
# based on visual inspection:
# MCF10A_WT_[123] removed
# MCF10A_Scramble_3 removed 
# MCF10A_shHP1b_3 removed due to low knock-down efficiency 
# MCF10A_shHP1ab_3 removed due to low knock-down efficiency

toRemove <- c("MCF10A_WT_1|MCF10A_WT_2|MCF10A_WT_3|MCF10A_shHP1ab_3|MCF10A_shHP1b_3|MCF10A_Scramble_3|MCF10A_shH2AZHP1a_3|MCF10A_shHP1a_3")
files <- files[-grep(toRemove, names(files))]
length(files)
condition <- condition[-grep(toRemove, names(condition))]
length(condition)
sample_id <- sample_id[-grep(toRemove, sample_id)]
kal_dirs <- kal_dirs[-grep(toRemove, kal_dirs)]

txi <- tximport::tximport(files, 
                          type = "kallisto",
                          tx2gene = t2g, geneIdCol = geneIdCol, txIdCol = targetIdCol,
                          ignoreTxVersion = F)
sd1 <- apply(txi$abundance, 1, sd)
summary(sd1)
pca1 <- ade4::dudi.pca(t(txi$abundance[sd1 > 0, ]), scannf = F, nf = 6)
pdf(paste("PCA_MCF10A_HP1-alpha_post_sample_removal_", annotationVersion, ".pdf", sep = ""))
  ade4::s.arrow(pca1$li, boxes = F)
  ade4::s.class(pca1$li, fac = as.factor(condition))
dev.off()

# prepare data.frame for importing kallisto data into sleuth --------------
s2c <- data.frame(sample = sample_id, condition = condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
s2c$sample <- as.character(s2c$sample)


# build a list for all contrasts to test ----------------------------------
contrastsList <- runConfig$samples$`RNA-Seq`$sleuth_parameters$contrasts
conditionsList <- runConfig$samples$`RNA-Seq`$sleuth_parameters$conditions
s2c.list <- lapply(names(contrastsList), function(x){
  # working under the assumption that these are pairwise comparisons only!
  return(s2c[grep(paste(unlist(conditionsList[contrastsList[[x]][c(1,2)]]), collapse = "|"), s2c$sample),]) 
})
names(s2c.list) <- names(contrastsList)
if(!file.exists(sleuth_results_output)){
  results <- lapply(names(s2c.list), function(x){
    s2c.list[[x]]$condition <- droplevels(s2c.list[[x]]$condition)
    design <- model.matrix(~ condition, data = s2c.list[[x]])
    print(paste("Processing ", x, " transcript-level analysis",sep = ""))
    #-----------------------------------------------------------------------------
    # transcript-level DE
    so <- sleuth::sleuth_prep(s2c.list[[x]],
                              ~ condition, 
                              target_mapping = t2g, 
                              max_bootstrap = 30,
                              read_bootstrap_tpm = T,
                              extra_bootstrap_summary = T)
    so <- sleuth::sleuth_fit(so, formula = design)
    so <- sleuth::sleuth_fit(so, ~1, "reduced")
    so <- sleuth::sleuth_lrt(so, "reduced", "full")
    for (i in colnames(design)[grep("Intercept", colnames(design), invert = T)]){
      so <- sleuth::sleuth_wt(so, i)  
    }
    rt.list <- lapply(colnames(design)[grep("Intercept", colnames(design), invert = T)], function(x){
      rt <- sleuth::sleuth_results(so, x)
      rt <- data.table::data.table(rt[order(rt$qval),])
    })
    names(rt.list) <- colnames(design)[grep("Intercept", colnames(design), invert = T)]
    kt <- sleuth::kallisto_table(so, normalized = T, include_covariates = T)
    kt_wide <- tidyr::spread(kt[, c("target_id", "sample", "tpm")], sample, tpm)
    kt_wide <- data.table::as.data.table(kt_wide)
    data.table::setkey(kt_wide, "target_id")
    # summarise gene level expression by adding all transcripts of a gene - the data.table way!
    so.gene <- sleuth::sleuth_prep(s2c.list[[x]], ~ condition, target_mapping = t2g, aggregation_column = aggregation_column)
    so.gene <- sleuth::sleuth_fit(so.gene, formula = design)
    so.gene <- sleuth::sleuth_fit(so.gene, ~1, "reduced")
    so.gene <- sleuth::sleuth_lrt(so.gene, "reduced", "full")
    for (i in colnames(design)[grep("Intercept", colnames(design), invert = T)]){
      so.gene <- sleuth::sleuth_wt(so.gene, i)  
    }
    rt.gene.list <- lapply(colnames(design)[grep("Intercept", colnames(design), invert = T)], function(x){
      rt.gene <- sleuth::sleuth_results(so.gene, x)
      rt.gene <- data.table::data.table(rt.gene[order(rt.gene$qval),])
      rt.gene <- rt.gene[!duplicated(rt.gene[,-c("transcript_version")])]
    })
    names(rt.gene.list) <- colnames(design)[grep("Intercept", colnames(design), invert = T)]
    # gene-level expression is summed from transcript level data (sum(TPM))
    target_mapping <- data.table::as.data.table(so$target_mapping)
    kt.gene <- sleuth::kallisto_table(so.gene, use_filtered = T)
    kt_wide.gene <- tidyr::spread(kt.gene[, c("target_id", "sample", "tpm")], sample, tpm)
    kt_wide.gene <- data.table::as.data.table(merge(kt_wide.gene, subset(t2g, select = c("target_id", aggregation_column)), all.x = TRUE, all.y = FALSE))
    cols.chosen <- as.character(s2c.list[[x]]$sample)
    kt_wide.gene <- kt_wide.gene[,lapply(.SD,sum),by=ens_gene, .SDcols = cols.chosen]
    kt.gene <- data.table::as.data.table(reshape::melt(kt_wide.gene))
    kt.gene$groups <- unlist(lapply(strsplit(as.character(kt.gene$variable), "_"), function(x) paste(x[1:2], collapse = "_")))
    kt_wide.gene <- data.table::as.data.table(merge(kt_wide.gene, unique(subset(ensGenes, select =  c("ensembl_gene_id", "external_gene_name", "description")), by = "ensembl_gene_id"), by.x = "ens_gene", by.y = "ensembl_gene_id", all.x = TRUE, all.y = FALSE))
    return(list(sleuth_object = so,
                sleuth_object_genes = so.gene,
                sleuth_results = rt.list,
                sleuth_results_genes = rt.gene.list,
                kallisto_table = kt,
                kallisto_table_wide = kt_wide,
                kallisto_table_genes = kt.gene,
                kallisto_table_genes_wide = kt_wide.gene))
  }) 
  names(results) <- names(s2c.list)
  save(results, file = sleuth_results_output)
} else {
  load(sleuth_results_output)
}


# export sleuth_results_genes as CSV --------------------------------------
lapply(names(results), function(x){
  tab1 <- data.table::as.data.table(results[[x]]$sleuth_results_genes)
  colnames(tab1) <- unlist(lapply(strsplit(colnames(tab1), "\\."), function(x) x[2]))
  fn <- paste(x, "_coding_sequences_results_table.csv", sep = "")
  write.csv(x = tab1, file = fn, row.names = FALSE)
})







