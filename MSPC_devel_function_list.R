## MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
#
#  imporoved version of custom function list that contributed to MSPC Packages
#  developed version: still under developing and optimize the code make it clean and no error happen

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

get.minOv.param <- function(allPeak.files, replicate.type=c("Biological", "Technical"), verbose=FALSE) {
  # input param check
  stopifnot(length(allPeak.files)>=2)
  replicate.type = match.arg(replicate.type)
  min.c <- ifelse(replicate.type=="Biological",
                  param <- length(allPeak.files)-1,
                  param <- length(allPeak.files))
  min.c <- as.integer(min.c)
  return(min.c)
}
##=====================================================================
.keep.peaks <- function(ov, all.peakFiles, replicate.type=c("Biological", "Technical")) {
  nn <- tot.ovnum(ov)
  min.c <- get.minOv.param(all.peakFiles, "Biological")
  keep <- nn >= min.c
  res <- lapply(ov, function(ele_) {
    ans <- ele_[keep]
  })
  message("keep peaks that sufficiently overlapped")
  return(res)
}

#=====================================================================

.discard.peaks_insufficientOverlap <- function(ov, all.peakFiles, replicate.type=c("Biological", "Technical"), idx=1L) {
  #check input param
  nn <- tot.ovnum(ov)
  min.c <- get.minOv.param(all.peakFiles, "Biological")
  disc.idx <- nn < min.c
  ans <- lapply(ov, function(ele_) {
    out <- ele_[disc.idx]
    out <- Filter(length, out)
  })
  res <- mapply(exp_fun, ans, all.peakFiles[c(idx, seq_len(length(allFiles_denoise))[-idx])])
  return(res)
}

##=============================================================================

.get.pvalue <- function(ov.hit, ov.sample, idx, verbose=FALSE) {
  if (verbose) {
    cat(">> getting kardinality for all overlapped peaks from replicates...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  .get.pval_helper <- function(ov.hit, ov.sample, verbose) {
    res <- extractList(ov.sample$p.value, ov.hit)
    return(res)
  }
  ov.pass <- .keep.peaks(ov.hit, ov.sample, "Biological")
  allPeaks.pval <- mapply(.get.pval_helper, ov.pass, ov.sample[c(idx, seq_len(length(ov.sample))[-idx])])
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

.confirmed.peaks <- function(expanded.peaks, comb.stringThreshold=1E-8, verbose=TRUE) {
  # check input param
  stopifnot(is.numeric(comb.stringThreshold))
  stopifnot(inherits(expanded.peaks[[1]], c("GRanges", "data.frame")))
  if(!"comb.pvalue" %in% colnames(mcols(expanded.peaks[[1]]))) {
    stop("combined pvalue of peaks are not found, unable to proceed")
  }
  res <- lapply(expanded.peaks, function(ele_) {
    ans <- subset(ele_, ele_$comb.pvalue <= comb.stringThreshold)
    ans <- ans[!duplicated(ans),]
    ans
  })
  return(res)
}

##=======================================================================================
.discard.peaks_failedFisherTest <- function(expanded.peaks, comb.stringentThreshold=1E-8, verbose=TRUE) {
  # check input param
  stopifnot(inherits(expanded.peaks, c("GRanges", "data.frame")))
  stopifnot(is.numeric(comb.stringentThreshold))
  if(!"comb.pvalue" %in% colnames(mcols(expanded.peaks[[1]]))) {
    stop("combined pvalue of peaks are not found, unble to proceed")
  }
  res <- lapply(expanded.peaks, function(ele_) {
    disc <- subset(ele_, ele_$comb.pvalue > comb.stringentThreshold)
    disc <- disc[!duplicated(disc),]
    disc$comb.pvalue <- NULL
  })
  return(res)
}

##=======================================================================================

## Waiting problem : integrate discardedPeak from insufficient overlap and discardedPeaks from failed fisher test

# very sketch solution (very first sketch idea)

.discarded.peaks.summary <- mapply(rbind.data.frame, .discard.peaks_insufficientOverlap, .discard.peaks_failedFisherTest)

re_duplicate <- lapply(.discarded.peaks.summary, function(ele_) {
  ele_ <- ele_[!duplicated(ele_),]
  ele_
})

##=======================================================================================
