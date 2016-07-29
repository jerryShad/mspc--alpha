## MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
#
#  imporoved version of custom function list that contributed to MSPC Packages
#

##========================================================================
readBedAsGRanges<- function(peakFolder) {
  # input param checking
  if(missing(peakFolder)) {
    stop("input param is missing!")
  }
  files <- list.files(peakFolder, full.names = TRUE, "\\.bed$")
  f.read <- lapply(1:length(files), function(ele_) {
    res <- as(import.bed(files[ele_]), "GRanges")
  })
  return(f.read)
}

##===================================================================
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

##====================================================================
denoise_peakFile <- function(peakFolder, denoise_threshold=1E-4, verbose=FALSE) {
  if (verbose) {
    cat(">> filter out all background noise peaks from all replicates simultanously...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  stopifnot(is.numeric(denoise_threshold))
  if(!inherits(peakFolder[[1]], "GRanges")) {
    stop("file entry was not GRanges objects, invalid input")
  }
  filt <- lapply(peakFolder, function(ele_) {
    out <- subset(ele_, ele_$p.value < denoise_threshold)
  })
  res <- filt
  return(res)
}

##=====================================================================
overlap <- function(set, idx=1L, FUN=which.min) {
  chosen <- set[[idx]]
  que.hit <- as(findOverlaps(chosen), "List")
  sup.hit <- lapply(set[- idx], function(ele_) {
    ans <- as(findOverlaps(chosen, ele_), "List")
    out.idx0 <- as(FUN(extractList(ele_$p.value, ans)), "List")
    out.idx0 <- out.idx0[!is.na(out.idx0)]
    ans <- ans[out.idx0]
  })
  res <- c(list(que.hit), sup.hit)
  return(res)
}

##=====================================================================
tot.ovnum <- function(ov.hit, verbose=FALSE) {
  if (verbose) {
    cat(">> getting kardinality for all overlapped peaks from replicates...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  tot.num <- Reduce('+', lapply(ov.hit, lengths))
  return(tot.num)
}

##=====================================================================

init.filt <- function(ov.hit, allPeak_denosie, replicate.type=c("Biological", "Technical"), verbose=FALSE) {
  # input param checking
  replicate.type = match.arg(replicate.type)
  min.c <- ifelse(replicate.type=="Biological",
                  parLen <- length(allPeak_denosie) - 1L,
                  parlen <- length(allPeak_denosie))
  keep <- abs(tot.ovnum) >= abs(min.c)
  peaks_Pass <- lapply(ov.hit, function(ele_) {
    ans <- ele_[keep]
  })
  peaks_drop <- lapply(ov.hit, function(x) {
    message("all peaks regions are dropped because insufficient overlap")
    out <- x[!keep]
  })
  res <- list("peak_pass"=peaks_Pass,
              "peak_fail"=peaks_Pass)
  return(res)
}

##=============================================================================
.allPeak.pval.func <- function(ov.hit, ov.sample, idx, verbose=FALSE) {
  if (verbose) {
    cat(">> getting kardinality for all overlapped peaks from replicates...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  .get.pval_helper <- function(ov.hit, ov.sample, verbose) {
    res <- extractList(ov.sample$p.value, ov.hit)
    return(res)
  }
  allPeaks.pval <- mapply(.get.pval_helper, ov.hit, ov.sample[c(idx, seq_len(length(ov.sample))[-idx])])
  .pval.helper <- function(pv.li) {
    ans <- sapply(pv.li, function(x) {
      out <- ifelse(length(x)>0, x, 0)
    })
    return(ans)
  }
  res <- data.frame(mapply(.pval.helper, allPeaks.pval))
  return(res)
}

##============================================================================

get.fisherScore <- function(allPeask.pval, ...) {
  Fisher.score <- suppressWarnings(
    out <- apply(allPeask.pval[,], 1, function(ro) {
      ans <- sumlog(ro)$p
    })
  )
  return(Fisher.score)
}

##=============================================================================

allhit.expand <- function(ov.hit, ov.sample, idx, allPeak.val) {
  ov.expand <- ov.sample[unlist(extractList(seq_along(ov.sample), ov.hit))]
  all.expand <- mapply(ov.expand, ov.hit, ov.sample[c(idx, seq_len(length(ov.sample))[-idx])])
  corr.VecDim <- lapply(allPeak.val, function(ele_) {
    tmp_pval.df <- data.frame(cbind("p"=ele_, "comb.pvalue"=Fisher.score))
    out <- tmp_pval.df[tmp_pval.df$p != 0.000000e+00, ]
    out$p <- NULL
    out
  })
  res <- Map(cbind, all.expand, corr.VecDim)
  return(res)
}

##==============================================================================











