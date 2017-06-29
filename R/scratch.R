tb1 <- data.table::data.table(biomaRt::getBM(attributes = c("ensembl_transcript_id",
                                                                    "refseq_mrna"),
                                                    filter = "refseq_mrna",
                                                    value = mRNArefSeq_IDs,
                                                    mart = mart))

ncRNARefSeq <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                             "external_gene_name",
                                             "chromosome_name",
                                             "start_position",
                                             "end_position",
                                             "strand",
                                             "band",
                                             "description",
                                             "percentage_gc_content", # renames in current ensembl release
                                             "gene_biotype",
                                             "entrezgene",
                                             "refseq_ncrna_predicted",
                                             "refseq_ncrna"),
                              mart = mart)


data.table::data.table(biomaRt::getBM(attributes = c("ensembl_transcript_id",
                                                     "version",
                                                     "external_gene_name",
                                                     "entrezgene",
                                                     "refseq_mrna"),
                                      filter = "refseq_mrna",
                                      values = "NM_000014",
                                      mart = mart))

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
  dT <- data.table("refseq_mRNA" = IDs,
                   "entrez_id" = names(es1), 
                   "gene_symbol" = unlist(lapply(es1, function(x) x$name)),
                   "refseq_mRNA_version" = unlist(lapply(es, function(x) x$accessionversion)))
  return(dT)
})

library(data.table)
library(rentrez)
seq1 <- seq(1, length(mRNArefSeq_IDs), 499)
l1 <- lapply(seq_along(seq1), function(x) {
  print(x)
  if (x + 1 <= length(seq1)){ 
    IDs <- mRNArefSeq_IDs[seq1[x] : seq1[x+1]]
    es <- entrez_summary(db = "nuccore", id = IDs)
    el <- entrez_link(dbfrom="nuccore", id = names(es), db = "gene", by_id = T)
    nuccore_gene <- unlist(lapply(1:length(el), function(x) el[[x]]$links$nuccore_gene))
    es1 <- entrez_summary(db="gene", id=nuccore_gene)
    dT <- data.table("refseq_mRNA" = IDs,
                     "entrez_id" = names(es1), 
                     "gene_symbol" = unlist(lapply(es1, function(x) x$name)),
                     "refseq_mRNA_version" = unlist(lapply(es, function(x) x$accessionversion)))
    return(dT)
  }
})


