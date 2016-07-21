## MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
##
##' @title unit test for MSPC Packages
##' @author Julaiti Shayiding

##==========================================================================
## Mission: read peakFile folder as GRanges and do manipulation for score
## Status: perfectly works. SUCCEED
peakset <- readBedAsGRanges("data/TMP/")
peakFiles <- lapply(peakset, function(ele_) {
  ans <- pvalueConversion(ele_)
})

##===========================================================================
## Mission : finding overlapped regions across multiple GRanges simlatanously
## Status : Status: perfectly works. SUCCEED
ov <- lapply(peakFiles[-1L], function(ele_) {
  # remove noiseMyc transcription factor ChIP-seq dataset in K562 cell available in ENCODE project
  query_denoise <- peakFiles[[1L]] %>%  subset(p.value >= 1e-4) %>% unique
  ans <- as(findOverlaps(peakFiles[[1L]], ele_), "List")
  hit.idx <- as(which.max(extractList(ele_$p.value, ans)), "List")
  hit.idx <- hit.idx[!is.na(hit.idx)]
  ans[hit.idx]
})

##============================================================================
## Mission: parallel countning kardinality of total overlapped regions including chosen.gr
## Status: perfectly works. SUCCEED
tot.ovnum <- Reduce('+', lapply(ov, function(ele_) { 
  lengths(ele_)})) + countOverlaps(peakFiles[[1]])

##============================================================================
## Mission : initial filtering process
## Attension: Here I do not expand overlapped supprt.grs (a.k.a, peakFiles[-1L])
## Status: perfectly works. SUCCEED

keep <- tot.ovnum >= 2L
chosen.ok <- peakFiles[[1]][keep]
chosen.no <- peakFiles[[1]][!keep]

ov.ok <- lapply(ov, function(ele_) {
  ele_[keep]
})

##=============================================================================
## Mission: extracting pvalue of overlapped supprt.gr
## Status: perfectly works. SUCCEED
get_pvalue <- function(ov.hit, obj) {
  ans <- extractList(obj$p.value, ov.hit)
  return(ans)
}
pv.ovHit <- mapply(get_pvalue, ov.ok, peakFiles[-1L])

##================================================================================
## Mission : getting correct dimension of pvalue vector for overlapped support.grs
## Status: perfectly works. SUCCEED
tts <- function(li) {
  ans <- sapply(li, function(x) {
    res <- ifelse(length(x)>0,x,0)
  })
  return(ans)
}
pv.list <- mapply(tts, pv.ovHit)

##=================================================================================
## Mission : performing Fisher' method to obtain global combined p.value
##Status: perfectly works. SUCCEED
pval.df <- data.frame(cbind(chosen.ok$p.value, pv.list))
library(metap)
comb.pval <- suppressWarnings(
  res <- apply(pval.df[,],1, function(ro) sumlog(ro)$p)
)

##=================================================================================
## Mission : Already have correct dimension of combined pvalue vector, 
## now can remove extra comb.p for non-overlapped regions
## Status: perfectly works. SUCCEED

sup.pvlist <- lapply(pval.df[-1L], function(ele_) {
  ans <- data.frame(cbind("p"=ele_, "comb.p"=comb.pval))
  dat <- ans[ans$p != 0.000000e+00, ]
})

##===========================================================================================
## Mission : cancel out redundance pvalue vector, because it already exisit in peakFiles[-1L]
## Status: perfectly works. SUCCED
sup_comb.pList <- lapply(sup.pvlist, function(ele_) {
  ele_$p <- NULL
  ele_
})
##================================================================================
## Mission : expanding overlapped support.gr that passed from initial filtering
## Status: perfectly works, SUCCEED
exp_fun <- function(ov.hit, gr) {
  ans <- gr[unlist(extractList(seq_along(gr), ov.hit))]
  return(ans)
}
sup.expand <- mapply(exp_fun, ov.ok, peakFiles[-1L])

##==================================================================================
## Mission: assigning sup_comb.plist to sup.expand as new attributes
## Status : perfectly works, SUCCEED

sup.expand <- base::Map(cbind, sup.expand, sup_comb.pList)
chosen.ok$comb.p <- comb.pval

##==================================================================================
## Mission : working with filitering process based on stringency threshold parameter
## Status: make it better format 

all.dfs <- c(list(chosen.ok), sup.expand)

out <- lapply(all.dfs, function(x) {
  confirmed <- subset(x, x$comb.p <= 1e-08)
})



