# MSPC Package development - Bioconductor Package of Multiple Sample Peak Calling
#
#' @title findoverlapsFromMultipleSample 
#' @param chosen.gr chosen.replicates
#' @param support.gr list of support replicates
#' @param option.multiple.intersection 
#' @return IntegerList hit-index vector
#' @export
#' @importFrom GenomicRanges findOverlaps
#' @importFrom XVector extractList
#' @author Julaiti Shayiding
#' @example 

findoverlapsFromMultipleSample <- function(chosen.gr, support.gr, option=c("MostStringent", "LeastStringent")) {
  # parameter check
  stopifnot(class(chosen.gr)=="GRanges")
  stopifnot(inherits(support.gr, c("GRanges", "GRangesList")))
  option = match.arg(option)
  ov <- lapply(support.gr, function(ele_) {
    res <- as(findOverlaps(chosen.gr, ele_), "List")
    if(option=="MostStringent") {
      idx <- as(which.min(extractList(ele_$pvalue, res)), "List")
    } else {
      idx <- as(which.max(extractList(ele_$pvalue, res)), "List")
    }
    idx <- idx[!is.na(idx)]
    res[idx]
  })
  return(ov)
}

