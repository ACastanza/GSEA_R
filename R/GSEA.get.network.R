#' Constructs a gene weight network from HumanNet v2
#'
#' `GSEA.get.network` constructs a network weight table for all genes
#'
#' Internal `GSEA` function.
#'
#' @keywords internal
#'

GSEA.get.network <- function(msigdbversion, gene.labels) {
 library("igraph")
 net <- read.table(url("https://www.inetbio.org/humannet/networks/HumanNet-XN.tsv"), 
  header = FALSE, stringsAsFactors = FALSE, sep = "\t")
 chip <- read.table(url(paste0("https://data.broadinstitute.org/gsea-msigdb/msigdb/annotations_versioned/Human_NCBI_Entrez_Gene_ID_MSigDB.v", 
  msigdbversion, ".chip")), header = TRUE, stringsAsFactors = FALSE, sep = "\t", 
  quote = "", fill = TRUE)
 chip = chip[, -c(3)]
 labels.chip <- merge(x = gene.labels, y = chip, by.x = 1, by.y = 2)
 net_map = merge(x = net, y = labels.chip, by.x = 1, by.y = 2)
 net_map = net_map[, -c(1)]
 colnames(net_map)[3] <- "Gene.Symbol.1"
 net_map = merge(x = net_map, y = labels.chip, by.x = 1, by.y = 2)
 net_map = net_map[, -c(1)]
 colnames(net_map)[3] <- "Gene.Symbol.2"
 net_map = net_map[, c(2, 3, 1)]
 colnames(net_map)[3] <- "LLS"
 net_graph <- graph_from_data_frame(net_map, directed = FALSE, )
 net_weight <- strength(net_graph)
 net_weight <- log(net_weight)/median(log(na.omit(net_weight)))
 net_weight <- 1 + net_weight
 return(list(weights = net_weight, map = net_map))
}
