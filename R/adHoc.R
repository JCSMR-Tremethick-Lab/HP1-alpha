toCheck <- c("DGKI", "CADM3", "PADI1", "SCG5")
results$kallisto_table_genes[target_id %in% toCheck]
p1 <- ggplot(data = results$kallisto_table_genes[target_id %in% toCheck & condition %in% c("MCF10A_WT", "MCF10A_shH2AZHP1a")], aes(condition,tpm)) + geom_jitter(width = 0.1, aes(colour = sample)) + facet_wrap(~target_id, scales = "free_y")


function (gr = NULL, out_file = NULL) 
{
  if (is.null(gr)) {
    stop("GRanges object missing")
  }
  if (is.null(out_file)) {
    stop("Output file is missing")
  }
  st <- GenomicRanges::strand(gr)
  levels(st)[3] <- "."
  if (is.null(names(gr))) {
    df <- data.frame(seqnames = seqnames(gr), starts = as(start(gr) - 
                                                            1, "integer"), ends = as(end(gr), "integer"), scores = ".", 
                     strand = st)
  }
  else {
    df <- data.frame(seqnames = seqnames(gr), starts = as(start(gr) - 
                                                            1, "integer"), ends = as(end(gr), "integer"), names = names(gr), 
                     scores = ".", strand = st)
  }
  write.table(df, file = out_file, quote = F, sep = "\t", row.names = F, 
              col.names = F)
}



library(httr)
library(jsonlite)
library(xml2)

server <- "http://grch37.rest.ensembl.org"
ext <- "/regulatory/species/homo_sapiens/microarray/HumanWG_6_V2/vendor/illumina?"

r <- GET(paste(server, ext, sep = ""), content_type("application/json"))

stop_for_status(r)

# use this if you get a simple nested list back, otherwise inspect its structure
# head(data.frame(t(sapply(content(r),c))))
head(fromJSON(toJSON(content(r))))
