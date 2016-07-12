# MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
#
# list of helper function
#'
#'
##======================================================================
## ovelapHit function
#' @title overlapHit
#' @param chosen.gr
#' @param support.gr
#' @param option
#' @return IntegerList objects
#' @export
#' @importFrom GenomicRanges findOverlaps
#' @importFrom XVector extractList
#' @author Julaiti Shayiding
#' @example 

overlapHit <- function(chosen.gr, suspport.gr, 
                       option=c("most-stringnet", "Least-stringent"), ...) {
  # input parameter check
  stopifnot(class(chosen.gr)=="GRanges")
  stopifnot(inherits(support.gr, c("GRanges", "GRangesList")))
  option <- match.arg(option)
  res <- lapply(support.gr, function(ele_) {
    ov <- as(findOverlaps(chosen.gr, ele_), "List")
    if(option=="most-stringent") {
      index <- as(which.min(extractList(ele_$p.value, ov)), "List")
    } else {
      index <- as(which.max(extractList(ele_$p.value, ov)), "List")
    }
    index <- index[!is.na(index)]
    ov <- ov[index]
  })
  return(res)
}

#==================================================================================
# cnt.overlap
#'
#' @title cnt.ov.stats
#' @param ov ovelap hit-index
#' @param chosen.gr chosen query replicates
#' @return IntegerList vector
#' @export
#' @importFrom S4Vectors Reduce
#' @author Julaiti Shayiding
#' @example 

cnt.ov.stats <- function(ov, chosen.gr, verbose=FALSE) {
  # input parameter check
  tmp <- Reduce('+', function(ele_) {
    len <- length(ele_)
  })
  len <- length(chosen.gr)
  res <- tmp+len
  return(res)
}

##===========================================================================
#' @title init.filter
#' @param replicate.type
#' @param chosen.gr
#' @param support.gr
#' @export
#' @return list
#' @example 

init.filter <- function(chosen.gr, support.gr, replicate.type=c("Biological", "Technical"), ...) {
  # input param check
  overlap <- overlapHit(chosen.gr = chosen.gr,
                   suspport.gr = support.gr, ...)
  novlp <- cnt.ov.stats(ov = overlap, chosen.gr = chosen.gr)
  novlp <- abs(novlp)
  
  num.input <- function(x, y) {
    w <- length(list(x))
    v <- lapply(y, length)
    res <- sum(w,v)
    return(res)
  }
  nn <- num.input(chosen.gr, support.gr)
  if(replicate.type=="Biological") {
    min.c <- nn -1
    keep <- novlp >= min.c
    chosen.gr <- chosen.gr[keep]
    overlap.idx <- lapply(overlap, function(ele_) {
      res <- ele_[keep]
    })
    message("extracting pvalue for getting its combined pvalue")
    chosen.pval <- as(chosen.gr$p.value, "List")
    for(i in 1:length(support.gr)) {
      support.pval <- extractList(support.gr[[i]]$p.value, overlap.idx[[i]])
      support.pval <- lapply(support.pval, function(ele_) {
        res <- ifelse(length(ele_)>0, ele_, 0)
      })
      return(support.pval)
    }
    pval.df <- as.data.frame(DataFrame(chosen.pval, support.pval))
    require(metap)
    pval.df$comb.pval <- apply(pval.df,1, function(ro) sumlog(ro))
  }
  if(replicate.type=="Technical") {
    min.c <- nn
    keep <- novlp >= min.c
    chosen.gr <- chosen.gr[keep]
    overlap.idx <- lapply(overlap, function(ele_) {
      res <- ele_[keep]
    })
    
  }
}
