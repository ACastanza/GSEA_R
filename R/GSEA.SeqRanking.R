#' Ranks genes according to the specified ranking metric
#'
#' `GSEA.SeqRanking` computes the GSEA ranking metric for each gene in the gene raw counts matrix using DESeq2
#'
#' Compute the GSEA ranking metric for each gene in the gene list. 
#' This function ranks the genes for the actual phenotype and also random 
#' permutations and bootstrap subsamples of both the observed and random 
#' phenotypes. It uses matrix operations to implement the rank calculations 
#' in stages and achieves fast execution speed. It supports two types of 
#' permutations: random (unbalanced) and balanced.  It also supports 
#' subsampling and bootstrap by using masking and multiple-count variables.
#' When 'fraction' is set to 1 (default) the there is no subsampling or 
#' boostrapping and the matrix of observed rank metrics will have the
#' same value for all permutations. This is wasteful but allows to support all the
#' multiple options with the same code. Notice that the second matrix for the null
#' distribution will still have the values for the random permutations (null
#' distribution). This mode (fraction = 1.0) is the default, the recommended one
#' and the one used in the examples. It is also the one that has be tested more
#' thoroughly. The resampling and boostrapping options are intersting to obtain
#' smooth estimates of the observed distribution but its is left for the expert
#' user who may want to perform some sanity checks before trusting the code.
#' Inputs: A: Matrix of gene expression values (rows are genes, columns are
#' samples) class.labels: Phenotype of class disticntion of interest. A vector of
#' binary labels having first the 1's and then the 0's gene.labels: gene labels.
#' Vector of probe ids or accession numbers for the rows of the expression matrix
#' nperm: Number of random permutations/bootstraps to perform permutation.type:
#' Permutation type: 0 = unbalanced, 1 = balanced. For experts only (default: 0)
#' sigma.correction: Correction to the signal to noise ratio (Default =
#' GeneCluster, a choice to support the way it was handled in a previous package)
#' fraction: Subsampling fraction. Set to 1.0 (no resampling). For experts only
#' (default: 1.0) replace: Resampling mode (replacement or not replacement). For
#' experts only (default: F) reverse.sign: Reverse direction of gene list (default
#' = F) rank.metric: metric to use for ranking genes, supports 'change' which ranks 
#' by the DESeq2 computed Log2(FC), 'signedsig' which ranks by the  * the Sign of 
#' the Log2(FC), and 'scaledchange' which ranks by the DESeq2 computed 
#' Log2(FC)*-log10(DESeq2 pValue), Outputs: rnk.matrix: Matrix with 
#' random permuted or bootstraps rank metrics signal to noise ratios by default 
#' (rows are genes, columns are permutations or bootstrap subsamplings 
#' obs.rnk.matrix: Matrix with observed rank metrics (rows are genes, columns 
#' are boostraps subsamplings. If fraction is set to 1.0 then all the columns 
#' have the same values order.matrix: Matrix with the orderings that will sort 
#' the columns of the obs.rnk.matrix in decreasing rnk order obs.order.matrix: 
#' Matrix with the orderings that will sort the columns of the rnk.matrix in 
#' decreasing rnk order.
#'
#' @keywords internal
#'

GSEA.SeqRanking <- function(A, class.labels, gene.labels, nperm, permutation.type = 0, 
 fraction = 1, replace = F, reverse.sign = F, rank.metric, progress, total, stage) {
 
 mode(A) <- "integer"
 
 N <- length(A[, 1])
 Ns <- length(A[1, ])
 
 subset.mask <- matrix(0, nrow = Ns, ncol = nperm)
 reshuffled.class.labels1 <- matrix(0, nrow = Ns, ncol = nperm)
 reshuffled.class.labels2 <- matrix(0, nrow = Ns, ncol = nperm)
 class.labels1 <- matrix(0, nrow = Ns, ncol = nperm)
 class.labels2 <- matrix(0, nrow = Ns, ncol = nperm)
 
 order.matrix <- matrix(0, nrow = N, ncol = nperm)
 obs.order.matrix <- matrix(0, nrow = N, ncol = nperm)
 rnk.matrix <- matrix(0, nrow = N, ncol = nperm)
 obs.rnk.matrix <- matrix(0, nrow = N, ncol = nperm)
 
 obs.gene.labels <- vector(length = N, mode = "character")
 obs.gene.descs <- vector(length = N, mode = "character")
 obs.gene.symbols <- vector(length = N, mode = "character")
 
 M1 <- matrix(0, nrow = N, ncol = nperm)
 M2 <- matrix(0, nrow = N, ncol = nperm)
 S1 <- matrix(0, nrow = N, ncol = nperm)
 S2 <- matrix(0, nrow = N, ncol = nperm)
 
 gc()
 
 C <- split(class.labels, class.labels)
 class1.size <- length(C[[1]])
 class2.size <- length(C[[2]])
 class1.index <- seq(1, class1.size, 1)
 class2.index <- seq(class1.size + 1, class1.size + class2.size, 1)
 
 for (r in 1:nperm) {
  class1.subset <- sample(class1.index, size = ceiling(class1.size * fraction), 
   replace = replace)
  class2.subset <- sample(class2.index, size = ceiling(class2.size * fraction), 
   replace = replace)
  class1.subset.size <- length(class1.subset)
  class2.subset.size <- length(class2.subset)
  subset.class1 <- rep(0, class1.size)
  for (i in 1:class1.size) {
   if (is.element(class1.index[i], class1.subset)) {
    subset.class1[i] <- 1
   }
  }
  subset.class2 <- rep(0, class2.size)
  for (i in 1:class2.size) {
   if (is.element(class2.index[i], class2.subset)) {
    subset.class2[i] <- 1
   }
  }
  subset.mask[, r] <- as.numeric(c(subset.class1, subset.class2))
  fraction.class1 <- class1.size/Ns
  fraction.class2 <- class2.size/Ns
  
  if (permutation.type == 0) {
   # random (unbalanced) permutation
   full.subset <- c(class1.subset, class2.subset)
   label1.subset <- sample(full.subset, size = Ns * fraction.class1)
   reshuffled.class.labels1[, r] <- rep(0, Ns)
   reshuffled.class.labels2[, r] <- rep(0, Ns)
   class.labels1[, r] <- rep(0, Ns)
   class.labels2[, r] <- rep(0, Ns)
   for (i in 1:Ns) {
    m1 <- sum(!is.na(match(label1.subset, i)))
    m2 <- sum(!is.na(match(full.subset, i)))
    reshuffled.class.labels1[i, r] <- m1
    reshuffled.class.labels2[i, r] <- m2 - m1
    if (i <= class1.size) {
      class.labels1[i, r] <- m2
      class.labels2[i, r] <- 0
    } else {
      class.labels1[i, r] <- 0
      class.labels2[i, r] <- m2
    }
   }
  } else if (permutation.type == 1) {
   # proportional (balanced) permutation
   
   class1.label1.subset <- sample(class1.subset, size = ceiling(class1.subset.size * 
    fraction.class1))
   class2.label1.subset <- sample(class2.subset, size = floor(class2.subset.size * 
    fraction.class1))
   reshuffled.class.labels1[, r] <- rep(0, Ns)
   reshuffled.class.labels2[, r] <- rep(0, Ns)
   class.labels1[, r] <- rep(0, Ns)
   class.labels2[, r] <- rep(0, Ns)
   for (i in 1:Ns) {
    if (i <= class1.size) {
      m1 <- sum(!is.na(match(class1.label1.subset, i)))
      m2 <- sum(!is.na(match(class1.subset, i)))
      reshuffled.class.labels1[i, r] <- m1
      reshuffled.class.labels2[i, r] <- m2 - m1
      class.labels1[i, r] <- m2
      class.labels2[i, r] <- 0
    } else {
      m1 <- sum(!is.na(match(class2.label1.subset, i)))
      m2 <- sum(!is.na(match(class2.subset, i)))
      reshuffled.class.labels1[i, r] <- m1
      reshuffled.class.labels2[i, r] <- m2 - m1
      class.labels1[i, r] <- 0
      class.labels2[i, r] <- m2
    }
   }
  }
 }
 
 if (rank.metric == "change" | rank.metric == "signedsig" | rank.metric == "scaledchange") {
  library(DESeq2)
  coldata <- as.data.frame(colnames(A), stringsAsFactors = FALSE)
  rownames(coldata) <- coldata[, 1]
  colnames(coldata) <- "condition"
  coldata.rand <- coldata
  rownames(rnk.matrix) <- rownames(A)
  rownames(obs.rnk.matrix) <- rownames(A)
  
  if (stage == "permute") {
   for (d in 1:nperm) {
    
    message("Computing permutation ", progress + d - 1, " of ", total, 
      "...")
    
    coldata.rand[, 1] <- reshuffled.class.labels1[, d]
    coldata.rand$condition <- as.factor(coldata.rand$condition)
    dds <- DESeqDataSetFromMatrix(countData = A, colData = coldata.rand, 
      design = ~condition)
    dds <- DESeq(dds)
    res <- results(dds)
    if (rank.metric == "change") {
      rnk.matrix[, d] <- res[, 2]  #rank by Log2(FC)
    }
    if (rank.metric == "scaledchange") {
      rnk.matrix[, d] <- res[, 2] * -log10(res[, 5])  #rank by Log2(FC)*-log10(pValue)
    }
    if (rank.metric == "signedsig") {
      rnk.matrix[, d] <- sign(res[, 2]) * -log10(res[, 5])  #rank by the -log10(pValue) signed by the Log2(FC)
    }
    gc()
    
    if (length(rnk.matrix[is.na(rnk.matrix)]) > 0) {
      warning(print(length(rnk.matrix[is.na(rnk.matrix)])), " N/A values were found in the permuted rank matrix. Setting N/As to Zero because these cause GSEA to fail.")
      rnk.matrix[is.na(rnk.matrix)] <- 0
    }
    if (fraction < 1) {
      
      coldata.obs <- coldata
      coldata.obs[, 1] <- class.labels1[, 1]
      coldata.obs$condition <- as.factor(coldata.obs$condition)
      dds <- DESeqDataSetFromMatrix(countData = A, colData = coldata.obs, 
     design = ~condition)
      dds <- DESeq(dds, quiet = TRUE)
      res <- results(dds)
      if (rank.metric == "change") {
     obs.rnk.matrix[, d] <- res[, 2]  #rank by Log2(FC)
      }
      if (rank.metric == "scaledchange") {
     obs.rnk.matrix[, d] <- res[, 2] * -log10(res[, 5])  #rank by Log2(FC)*-log10(pValue)
      }
      if (rank.metric == "signedsig") {
     obs.rnk.matrix[, d] <- sign(res[, 2]) * -log10(res[, 5])  #rank by the -log10(pValue) signed by the Log2(FC)
      }
      gc()
      
      if (length(obs.rnk.matrix[is.na(obs.rnk.matrix)]) > 0) {
     warning(print(length(obs.rnk.matrix[is.na(obs.rnk.matrix)])), 
       " N/A values were found in the observed ranked list. Setting N/As to Zero because these cause GSEA to fail.")
     obs.rnk.matrix[is.na(obs.rnk.matrix)] <- 0
      }
    }
   }
  }
  if (stage == "rank") {
   if (fraction == 1) {
    coldata.obs <- coldata
    coldata.obs[, 1] <- class.labels1[, 1]
    coldata.obs$condition <- as.factor(coldata.obs$condition)
    dds <- DESeqDataSetFromMatrix(countData = A, colData = coldata.obs, 
      design = ~condition)
    message("Computing actual gene rankings...")
    dds <- DESeq(dds)
    res <- results(dds)
    if (rank.metric == "change") {
      obs.rnk.matrix[, 1:nperm] <- res[, 2]  #rank by Log2(FC)
    }
    if (rank.metric == "scaledchange") {
      obs.rnk.matrix[, 1:nperm] <- res[, 2] * -log10(res[, 5])  #rank by Log2(FC)*-log10(pValue)
    }
    if (rank.metric == "signedsig") {
      obs.rnk.matrix[, 1:nperm] <- sign(res[, 2]) * -log10(res[, 5])  #rank by the -log10(pValue) signed by the Log2(FC)
    }
    gc()
    
    if (length(obs.rnk.matrix[is.na(obs.rnk.matrix)]) > 0) {
      warning(print(length(obs.rnk.matrix[is.na(obs.rnk.matrix)])), " N/A values were found in the observed ranked list. Setting N/As to Zero because these cause GSEA to fail.")
      obs.rnk.matrix[is.na(obs.rnk.matrix)] <- 0
    }
   }
  }
 }
 
 
 if (reverse.sign == T) {
  obs.rnk.matrix <- -obs.rnk.matrix
 }
 
 if (fraction < 1) {
  for (r in 1:nperm) {
   obs.order.matrix[, r] <- order(obs.rnk.matrix[, r], decreasing = T)
  }
 } else if (fraction == 1) {
  obs.order.matrix[, 1:nperm] <- order(obs.rnk.matrix[, 1], decreasing = T)
 }
 
 if (reverse.sign == T) {
  rnk.matrix <- -rnk.matrix
 }
 for (r in 1:nperm) {
  order.matrix[, r] <- order(rnk.matrix[, r], decreasing = T)
 }
 
 return(list(rnk.matrix = rnk.matrix, obs.rnk.matrix = obs.rnk.matrix, order.matrix = order.matrix, 
  obs.order.matrix = obs.order.matrix))
}
