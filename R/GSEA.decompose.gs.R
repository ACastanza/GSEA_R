#' Constructs a gene weight network from HumanNet v2 for a gene set
#'
#' `GSEA.decompose.gs` extracts subnetworks from gene sets to perform more 
#' granular subgraph level enrichment analysis
#'
#' Internal `GSEA` function.
#'
#' @keywords internal
#'

GSEA.decompose.gs <- function(gsdb, gs.names, gs.desc, gene.labels, gs.size.threshold.min, 
 msigdbversion) {
 
 gs.subgraphs <- list()
 gs.subgraph.names <- list()
 gs.subgraph.desc <- list()
 
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
 global.network = net_map
 global_graph <- graph_from_data_frame(global.network, directed = FALSE)
 
 for (i in 1:length(gs.names)) {
  set <- gsdb[i, gsdb[i, ] != "null"]
  gene.set.graph <- induced.subgraph(graph = global_graph, vids = set[set %in% 
   V(global_graph)$name])
  subsets <- decompose(gene.set.graph, min.vertices = gs.size.threshold.min)
  if (length(subsets) > 0) {
   for (n in 1:length(subsets)) {
    gs.subgraphs[[length(gs.subgraphs) + 1]] <- V(subsets[[n]])
    gs.subgraph.names[length(gs.subgraph.names) + 1] <- paste0(gs.names[i], 
      "_subgraph_", names(which.max(eigen_centrality(subsets[[n]])$vector)), 
      "_centered")
    gs.subgraph.desc[length(gs.subgraph.desc) + 1] <- paste0("Subgraph of: ", 
      gs.desc[i])
   }
  }
 }
 
 gs.subgraph.names <- make.names(unlist(gs.subgraph.names))
 gs.subgraph.desc <- unlist(gs.subgraph.desc)
 return(list(names = gs.subgraph.names, desc = gs.subgraph.desc, gs = gs.subgraphs))
}
