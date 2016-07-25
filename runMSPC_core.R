peakset <- readBedAsGRanges("data/TMP/")
peakFiles <- lapply(peakset, function(ele_) {
  ans <- pvalueConversion(ele_)
})

denoise_threshold <- 1E-8
overlap <- function(peakFiles, idx=1L, denoise_threshold, FUN=which.min) {
  # check input param
  stopifnot(is.numeric(idx))
  stopifnot(is.numeric(denoise_threshold))
  chosen <- peakFiles[[idx]]
  chosen.denoise <- subset(chosen, chosen$p.value < denoise_threshold)
  chosen.hit <- as(findOverlaps(chosen.denoise), "List")
  support.hit <- lapply(peakFiles[- idx], function(ele_) {
    out <- as(findOverlaps(chosen.denoise, ele_), "List")
    out.idx0 <- as(FUN(extractList(ele_$p.value, out)), "List")
    out.idx0 <- out.idx0[!is.na(out.idx0)]
    res <- out[out.idx0]
  })
  res <- c(list(chosen.hit), support.hit)
  return(res)
}

tot.ovnum <- function(all.hit, verbose=FALSE) {
  if(verbose) {
    cat(">> obtaining total number of overlapped regions...\t\t",
        format(Sys.time(), "Y-%MM-%DD %X"), "\n")
  }
  kardinality <- Reduce('+', lapply(all.hit, lengths))
  return(kardinality)
}

hit.all_pass <- function(ov.hit, replicate.type=c("Biological", "Technical"), peakFolder, verbose=FALSE) {
  # input param checking
  if(verbose) {
    cat(">> obtaining param for minimum required overlap number...\t\t",
        format(Sys.time(), "Y-%MM-%DD %X"), "\n")
  }
  min.c <- ifelse(replicate.type == "Biological",
                  len <- length(peakFolder) - 1L,
                  len <- length(peakFolder))
  keep <- tot.ovnum >= min.c
  .all.pass <- lapply(ov.hit, function(ele_) {
    out <- ele_[keep]
  })
  res <- .all.pass
  return(res)
}

.get_hitAll.pval <- function(ov.hit, ov.obj, verbose) {
  if(verbose) {
    cat(">> obtaining total number of overlapped regions...\t\t",
        format(Sys.time(), "Y-%MM-%DD %X"), "\n")
  }
  res <- extractList(ov.obj$p.value, ov.hit)
  return(res)
}

.hitAll.pval <- mapply(.get_hitAll.pval, hit.all_pass, peakFiles)

.pval.helper <- function(pv.li, verbose=FALSE) {
  ans <- sapply(pv.li, function(x) {
    out <- ifelse(length(x)>0, x, 0)
  })
  return(ans)
}


pval.all <- data.frame(mapply(.pval.helper, pval.list))
library(metap)
comb.p <- suppressWarnings(
  res <- apply(pval.all[,], 1, function(ro) sumlog(ro)$p)
)

func <- lapply(pval.all, function(ele_) {
  ans <- data.frame(cbind("p"=ele_, "comb.p"=comb.p))
  dat <- ans[ans$p != 0.000000e+00, ]
  dat$p <- NULL
  dat
})

exp_fun <- function(ov.hit, gr) {
  ans <- gr[unlist(extractList(seq_along(gr), ov.hit))]
  return(ans)
}

sup.expand <- mapply(exp_fun, all.ok, peakFiles)
sup.expand <- base::Map(cbind, sup.expand, func)


out.conf <- lapply(sup.expand, function(ele_) {
  ans <- subset(ele_, ele_$comb.p <= 1E-8)
})

out.disc <- lapply(sup.expand, function(ele_) {
  ans <- subset(ele_, ele_$comb.p > 1E-8)
})
