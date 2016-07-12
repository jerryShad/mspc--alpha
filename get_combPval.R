# MSPC Project -  Bioconductor Package For Multiple Sample Peak Calling
#'
#' set of utility functions
#' 
#' @title get.CombPval
#' @param chosen.gr
#' @param support.gr
#' @return vector
#' @export
#' @importFrom metap sumlog
#' @author Julaiti Shayiding

get.combPval <- function(chosen.gr, support.gr, ov, output=NULL, ...) {
  # input parameter check
  stopifnot(class(chosen.gr)=="GRanges")
  stopifnot(inherits(support.gr, c("GRanges", "GRangesList")))
  if(!missing(chosen.gr$p.value) & !missing(support.gr$p.value)) {
    q.pval <- chosen.gr$pvalue
    q.pval <- unlist(q.pval)
    if(!length(support.gr)=0) {
      s.pval <- lapply(support.gr, function(ele_) {
        res <- ele_$p.value
        res <- sapply(res, function(ro) {
          tmp <- ifelse(length(x)>0,x,0)
        })
        unlist(res)
      })
      s.pval
    }
    dat <- data.frame(q.pval, s.pval)
    require(metap)
    comb.p <- suppressWarnings(
      res <- apply(dat[, 1:ncol(dat)], 1, function(ro) sumlog(ro)$p)
    )
    comb.p <- unlist(comb.p)
  }
  res <- comb.p
  return(res)
}

