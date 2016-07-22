# MSPC Project - Multiple Sample Peak Calling
#
#' @description 
#' @title 
#' @param 
#' @param 
#' @return 
#' @export
#' @importFrom 
#' @author 
#' @example 

runMSPC <- function(peakFolder, replicate.type=c("Biological", "Technical"), 
                    keep.Onlyone.option=c("most_stringent_peaks", "least_stringent_peaks"), denoise_threshold=1e-4, ...) {
  # check input parameter
  stopifnot(length(peakFolder)>=2L)
  stopifnot(is.numeric(denoise_threshold))
  replicate.type <- match.arg(replicate.type)
  message("reading peakFile foler")
  peak.grs <- readBedAsGRanges(peakFolder = peakFolder)
  
  if(!"p.value" %in% colnames(mcols(peak.grs[[1L]]))) {
    ans <- lapply(peak.grs, function(ele_) {
      res <- pvalueConversion(ele_)
    })
    return(ans)
  }
  chosenPeak <- peak.grs[[1L]]
  ov <- lapply(peak.grs[-1L], function(ele_) {
    query_denoise <- subset(chosenPeak, chosenPeak$p.value >= denoise_threshold)
    hit <- as(findOverlaps(query_denoise, ele_), "List")
    if(keep.Onlyone.option=="most_stringent_peaks") {
      idx0 <- as(which.min(extractList(ele_$p.value, hit)), "List")
      idx0 <- idx0[!is.na(idx0)]
    } else {
      idx0 <- as(which.max(extractList(ele_$p.value, hit)), "List")
      idx0 <- idx0[!is.na(idx0)]
    }
    res <- hit[idx0]
  })
  tot.ovnum <- function(ov, chosenPeak) {
    ans <- Reduce('+', lapply(ov, function(ele_) { 
      lengths(ele_)})) + countOverlaps(chosenPeak)
  }
  min.c <- ifelse(replicate.type=="Biological",
                  bio <- length(peakFolder)-1,
                  tech <- length(peakFolder))
  keep <- tot.ovnum >= min.c
  init.filter <- function(chosenPeak, ov) {
    chosenPeak <- chosenPeak[keep]
    hit <- lapply(ov, function(ele_) {
      ans <- ele_[keep]
    })
    res <- list('chosenPeak'=chosenPeak,
                'ov_hit'=hit)
    return(res)
  }
  # proceed rest of the task 
  
}
