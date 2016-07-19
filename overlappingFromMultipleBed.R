# MSPC Project -  Bioconductor Package of Multiple Sample Peak Calling
#
#' @title multipleBed_intersection
#' @param peakFolder
#' @param replicate.type
#' @param keepOnlyOne.option
#' @param denoise.threshold
#' @param FUN
#' @return IntegerList vector
#' @export
#' @importFrom GenomicRanges findOverlaps
#' @author Julaiti Shayiding

overlappingFromMultipleBeds <- function(peakFolder, denoise.threshold=1e-4, keepOnlyOne.option=c("mostStringentPeak", "LeastStringentPeaks"), verbose=FALSE) {
  # input param check
  stopifnot(length(peakFolder)>=2)
  stopifnot(is.numeric(denoise.threshold))
  if(!is(peakFolder[1], "GRanges")) {
    peaks.gr <- readBedAsGRanges(peakFolder)
    if(!"pvalue" %in% colnames(mcols(peaks.gr[[1]]))) {
      peaks.gr <- lapply(peaks.gr, function(ele_) {
        ans <- pvalueConversion(ele_)
      })
    }
  }
  peaks.gr <- sort(unique(peaks.gr))
  chosen.sample <- peaks.gr[[1]]
  support.sample <- peaks.gr[-1]
  message("we need to filter out background noise features from chosen.sample before proceed overlapping")
  chosen.sample_denoise <- chosen.sample[-(chosen.sample$pvalue >=denoise.threshold)]
  overlap.hit <- lapply(supportSample, function(ele_) {
    ans <- as(findOverlaps(chosen.sample_denoise, ele_), "List")
    if(keepOnlyOne.option=="mostStringentPeaks") {
      hit.idx <- as(which.min(extractList(ele_$p.value, ans)), "List")
      hit.idx <- hit.idx[!is.na(hit.idx)] 
    } else {
      hit.idx <- as(which.max(extractList(ele_$p.value, ans)), "List")
      hit.idx <- hit.idx[!is.na(hit.idx)]
    }
    ans[hit.idx]
  })
  return(overlap.hit)
}
