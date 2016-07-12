# MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
#
# utility function
#' 
#' 
#' @title expand.ovlp
#' @param ov.index
#' @param support.gr
#' @return GRanges
#' @export
#' @importFrom  XVector extractList
#' @author Julaiti Shayiding

expand.ovlp <- function(ov.index, chosen.gr, support.gr, ...) {
  # input parameter check
  stopifnot(class(chosen.gr)=="GRanges")
  stopifnot(inherits(support.gr, c("GRanges", "GRangesList")))
  stopifnot(missing(ov.index))
  stopifnot(class(ov.index)=="IntegerList")
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
  s.pval <- sort(s.pval)
  q.pval <- unlist(chosen.gr$p.value)
  comb.p <- get_combPal(chosen.gr, support.gr, ov.index)
  
  dat <- data.frame('p'= s.pval, 'cmb.p'=comb.p)
  dat <- dat[dat$p!=0e+00,]
  
  res <- lapply(support.gr, function(ele_) {
    for(i in 1:length(ov.index)) {
      gr <- ele_[unlist(extractList(seq_along(ele_), ov.index[[i]]))]
      gr$comb.p <- dat[[2]]
      gr
    }
    ele_
  })
  return(res)
}
