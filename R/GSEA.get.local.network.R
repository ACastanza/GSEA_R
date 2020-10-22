#' Constructs a gene weight network from HumanNet v2 for a gene set
#'
#' `GSEA.get.local.network` constructs a network weight table for the genes in a gene set
#'
#' Internal `GSEA` function.
#'
#' @keywords internal
#'

GSEA.get.local.network <- function(global.netowrk, set, score.type = "strength") {
 
 gene.set.frame <- as.data.frame(set, stringsAsFactors = FALSE)
 gene.set.map <- merge(x = global.netowrk, y = gene.set.frame, by.x = 2, by.y = 1)
 gene.set.map <- merge(x = gene.set.map, y = gene.set.frame, by.x = 2, by.y = 1)
 gene.set.graph <- graph_from_data_frame(gene.set.map, directed = FALSE)
 if (score.type == "strength") 
  {
   gene.set.weight <- strength(gene.set.graph)
   gene.set.weight <- gene.set.weight[set]
   names(gene.set.weight) <- set
   gene.set.weight <- 1 + (log(1 + (gene.set.weight/median(na.omit(gene.set.weight)))))
  # gene.set.weight <- 1 + ((gene.set.weight-min(gene.set.weight))/(max(gene.set.weight)-min(gene.set.weight)))
  }  #else if (score.type == 'centrality') {
 # gene.set.weight <- eigen_centrality(gene.set.graph)$vector gene.set.weight <-
 # gene.set.weight[set] names(gene.set.weight) <- set gene.set.weight <- 1 +
 # (log(1 + (gene.set.weight/median(na.omit(gene.set.weight))))) }
 gene.set.weight[is.na(gene.set.weight)] <- 1
 
 return(unname(as.character(gene.set.weight)))
}
