#' Read gene symbol mappings from .CHIP file
#'
#' `GSEA.ReadCHIPFile` is a wrapper for `read.table` with specified parameters for reading in a .CHIP file
#'
#' Internal `GSEA` function.
#'
#' @keywords internal
#'

GSEA.ReadCHIPFile <-
function(file = "NULL") {
chipframe <- read.table(file, header = T, 
stringsAsFactors = FALSE, sep = "\t", fill = TRUE, quote = "", na = "", comment.char="")
return(chipframe)
}
