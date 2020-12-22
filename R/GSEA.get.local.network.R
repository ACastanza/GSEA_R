#' Constructs a gene weight network from HumanNet v2 for a gene set
#'
#' `GSEA.get.local.network` constructs a network weight table for the genes in a gene set
#'
#' Internal `GSEA` function.
#'
#' @keywords internal
#'

GSEA.get.local.network <- function(global.network, set, score.type = "strength", 
 weighted.score.type, expression.matrix = NULL) {
 
 gene.set.frame <- as.data.frame(set, stringsAsFactors = FALSE)
 gene.set.map <- merge(x = global.network, y = gene.set.frame, by.x = 2, by.y = 1)
 gene.set.map <- merge(x = gene.set.map, y = gene.set.frame, by.x = 2, by.y = 1)
 gene.set.graph <- graph_from_data_frame(gene.set.map, directed = FALSE)
 if (score.type == "strength") {
  gene.set.weight <- strength(gene.set.graph)
  gene.set.weight <- as.numeric(weighted.score.type) + (log(1 + (gene.set.weight/median(na.omit(gene.set.weight)))))
  # gene.set.weight <- 1 +
  # ((gene.set.weight-min(gene.set.weight))/(max(gene.set.weight)-min(gene.set.weight)))
 } else if (score.type == "centrality") {
  gene.set.weight <- eigen_centrality(gene.set.graph)$vector
  gene.set.weight <- as.numeric(weighted.score.type) + (log(1 + (gene.set.weight/median(na.omit(gene.set.weight)))))
 } else if (score.type == "tif") {
  # Topology Influence Factor adapted from Hung et al. PMID:20187943
  d <- distances(gene.set.graph)
  set.expression.matrix <- expression.matrix[rownames(d), ]
  if (length(set.expression.matrix) > 0) {
   pcc <- stats::cor(t(set.expression.matrix), method = "pearson")#, use = "pairwise.complete.obs")
   d <- d[rownames(pcc), rownames(pcc)]
   f <- d/abs(pcc)
   diag(f) <- NA
   valid <- f <= -log(0.05)
   f <- as.data.frame(t(f))
   valid <- as.data.frame(t(valid))
   tif <- unlist(Map(function(x, y) if (sum(y, na.rm = TRUE) > 0) 
    mean(x[y], na.rm = TRUE) else NA, f, valid))
   tif <- as.numeric(weighted.score.type) + exp(-tif)
   gene.set.weight <- tif
  } else {
   gene.set.weight <- rep(as.numeric(weighted.score.type), length(set))
  }
 } else if (score.type == "stif") {
  # Topology Influence Factor adapted from Hung et al. PMID:20187943
  d <- distances(gene.set.graph)
  set.expression.matrix <- expression.matrix[rownames(d), ]
  if (length(set.expression.matrix) > 0) {
   pcc <- stats::cor(t(set.expression.matrix), method = "spearman")#, use = "pairwise.complete.obs")
   d <- d[rownames(pcc), rownames(pcc)]
   f <- d/abs(pcc)
   diag(f) <- NA
   valid <- f <= -log(0.05)
   f <- as.data.frame(t(f))
   valid <- as.data.frame(t(valid))
   tif <- unlist(Map(function(x, y) if (sum(y, na.rm = TRUE) > 0) 
    mean(x[y], na.rm = TRUE) else NA, f, valid))
   tif <- as.numeric(weighted.score.type) + exp(-tif)
   gene.set.weight <- tif
  } else {
   gene.set.weight <- rep(as.numeric(weighted.score.type), length(set))
  }
 }
 gene.set.weight <- gene.set.weight[set]
 names(gene.set.weight) <- set
 gene.set.weight[is.na(gene.set.weight)] <- as.numeric(weighted.score.type)
 
 return(unname(as.character(gene.set.weight)))
}
