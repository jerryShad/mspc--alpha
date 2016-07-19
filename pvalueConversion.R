#' pvalueConversion
#' @param peakFile input bed file
#' @param pvalueBase cofficient for different type of pvalue standard
#' @return GRanges objects
#' @export
#' @importFrom rtracklayer score
#' @author Julaiti Shayiding

pvalueConversion <- function(x, pvalueBase = 1L) {
  stopifnot(class(x) == "GRanges")
  # explore score of all features
  if(is.null(x$pvalue)){
    x$score <- 10^(score(x)/(- pvalueBase))
    colnames(mcols(x))[2] <- "p.value"
  } else {
    x
  }
  return(x)
}
