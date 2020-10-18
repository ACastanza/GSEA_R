#' Computes the enrichment score of a gene set
#'
#' `GSEA.Network.EnrichmentScore` computes the weighted GSEA score of gene.set in gene.list using the local set specific network to calcualte the weights for genes in the gene set and the global netowrk to clculate weights for the genes not in the gene set
#'
#' Internal `GSEA` function.
#' Computes the weighted GSEA score of gene.set in gene.list.  The weighted score
#' type is the exponent of the correlation weight, where the weight is prodived by the 
#' strength score of each gene in the local network for genes in the gene set or the 
#' strength score of genes in the global network for genes not in the gene set.  
#' Inputs: gene.list: The ordered gene list (e.g. integers indicating the original 
#' position in the input dataset) gene.set: A gene set (e.g. integers indicating the 
#' location of those genes in the input dataset). correl.vector: A vector with the 
#' coorelations (e.g. signal to noise scores) corresponding to the genes in the gene 
#' list. net.set: the netowrk scores for each gene in the gene set. correl.weight: the 
#' global network scores for all genes in the dataset. Outputs: ES: Enrichment score 
#' (real number between -1 and +1) arg.ES: Location in gene.list where the peak 
#' running enrichment occurs (peak of the 'mountain') RES: Numerical vector containing 
#' the running enrichment score for all locations in the gene list tag.indicator: 
#' Binary vector indicating the location of the gene sets (1's) in the gene list
#'
#' @keywords internal
#'

GSEA.Network.EnrichmentScore <- function(gene.list, gene.set, correl.vector = NULL, 
 net.set, correl.weight) {
 
 tag.indicator <- sign(match(gene.list, gene.set, nomatch = 0))  # notice that the sign is 0 (no tag) or 1 (tag)
 no.tag.indicator <- 1 - tag.indicator
 N <- length(gene.list)
 Nh <- length(gene.set)
 Nm <- N - Nh
 
 weighted.score.matrix <- correl.weight
 weighted.score.matrix[gene.set] <- net.set
 weighted.score.matrix <- weighted.score.matrix[gene.list]
 
 correl.vector2 <- abs(as.complex(correl.vector)^unname(weighted.score.matrix))
 sum.correl.tag <- sum(correl.vector2[tag.indicator == 1])
 norm.tag <- 1/sum.correl.tag
 norm.no.tag <- 1/Nm
 RES <- cumsum(tag.indicator * correl.vector2 * norm.tag - no.tag.indicator * 
  norm.no.tag)
 max.ES <- max(RES)
 min.ES <- min(RES)
 if (max.ES > -min.ES) {
  # ES <- max.ES
  ES <- signif(max.ES, digits = 5)
  arg.ES <- which.max(RES)
 } else {
  # ES <- min.ES
  ES <- signif(min.ES, digits = 5)
  arg.ES <- which.min(RES)
 }
 return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
}
