# MSPC Project
#'
#' @title pkg.core
#' @param chosenSample
#' @param targetSamples
#' @param replicate.type
#' @param tau.s
#' @param tau.w
#' @param option.multiple.intersection
#' @param comb.threshold
#' @param gamma
#' @return GRanges
#' @export
#' @importFrom GenomicRanges findOverlaps
#' @author Julaiti Shayiding
#' @example

pkg.core <- function(chosenSample, supportSample, replicate.type=c("Biological", "Technical"), tau.s=1e-10, tau.w=1e-4,
                     option.multiple.intersection=c("UseMoestStringentPeak", "UseLeastStringentPeak"),
                     comb.threshold=1e-8, gamma.threshold=0.05, pAdjustMethod="BH", outDir=NULL, ...) {
  # check input parameter
  flagGR <- NULL
  if(!inherits(chosenSample, "GRanges")) {
    message("input format must comply with GRanges objects")
    chosen.gr <- readAsGRanges(chosenSample)
  }
  chosen.gr <- unlist(sort(unique(chosen.gr)))
  stopifnot(inherits(supportSample, "GRangesList"))
  stopifnot(length(supportSample)>0)
  
  if(!is(supportSample[1], "GRangesList")) {
    File <- parse.supportSample.param(supportSample)
    support.gr <- lapply(File, readAsGRanges)
  } else {
    support.gr <- supportSample
    flagGR <- TRUE
  }
  if(!"pvalueLog" %in% colnames(mcols(support.gr[[1]]))) {
    support.gr$pvalueLog <- lapply(support.gr, pvalueConversion)
  }
  support.gr <- sort(unique(support.gr))
  
  replicate.type <- match.arg(replicate.type)
  stopifnot(is.numeric(tau.s))
  stopifnot(is.numeric(tau.w))
  
  # alias function
  num.input <- function(chosenSample, supportSample) {
    w <- length(list(chosenSample))
    v <- lapply(supportSample, length)
    res <- sum(w,v)
    return(res)
  }
  peak.classification <- function(chosen.gr, verbose=FALSE, ) {
    if(verbose) {
      cat(">> loading chosen GRanges object with metadata",
          format(Sys.time(), "%YYYY-%MM-DD%"), "\n")
    }
    if(!"pvalueLog" %in% colnames(mcols(chosen.gr))) {
      chosen.gr$pvalueLog <- pvalueConversion(chosen.gr)
      mcols(chosen.gr)[["pvalueLog"]] <- TRUE
    }
    chosen.stringent.ERs <- chosen.gr %>% subset(pvalue < tau.s) %>% unique %>% sort.GenomicRanges
    chosen.weak.ERs <- chosen.gr %>% subset(pvalue < tau.w & pvalue >= tau.s) %>% unique %>% sort.GenomicRanges
    chosen.background.noise <- chosen.gr %>% subset(pvalue >= tau.w) %>% unique %>% as.data.frame
    res <- GRangesList("AllStringentPeaks"=chosen.stringent.ERs, "AllWeakPeaks"=chosen.weak.ERs)
    return(res)
  }
  source.gr <- peak.classification
  message("perform on all stringent peaks to retrive overlapped regions from set of support sample in pair-wise")
  if(names(source.gr)=="AllStringentPeaks") {
    message("retrieving overlapping regions across set of support sample that performed on chosen sample")
    stringent <- S4Vectors::getListElement(source.gr, 1, exact = TRUE)
    #=============================================================================
    overlap.stringent <- lapply(support.gr, function(ele_) {
      ov <- as(findOverlaps(stringent, ele_), "List")
      hit.idx <- as(which.min(extractList(ele_$pvalue, ov)), "List")
      hit.idx <- hit.idx[!is.na(hit.idx)]
      ov[hit.idx]
    })
    
    #============================================================
    # get lengths of hitIndex
    hitIdx.len <- lapply(ov.hitIndex, lengths)   # must be vector
    hit.sum <- Reduce('+', hitIdx.len)
    kard <- sum(lengths(strigent), hit.sum)
    Kard <- Reduce('+', lapply(ov.s, function(ele_) { lengths(ele_)})) + lengths(foo.s)
    #============================================================
    if(replicate.type %in% c("Biological", "Technical")) {
      if(replicate.type=="Biological") {
        min.c <- num.input - 1
        keep <- kard >= min.c
        
        #=====================================================
        stringent_ <- stringent[keep]
        ov.hitIndex <- lapply(ov.hitIndex, function(ele_) {
          res <- ele_[keep]
        })
        #=====================================================
        # ready to expand ov.hitIndex for performing Fisher.test
        support.ov.gr <- lapply(ov.hitIndex, function(ele_) {
          su <- lapply(support.gr, function(gr) {
            idx <- unlist(extractList(seq_along(gr), ele_))
            gr[idx]
          })
        })
        #=====================================================
      } else {
        min.c <- num.input
        keep <- kard >= min.c
      }
    }

  }
  if(names(source.gr) == "AllWeakPeaks") {
    message("perform on all weak peaks to retrive overlapped regions from set of support sample in pair-wise")
    weak <- S4Vectors::getListElement(source.gr, 2, exact = TRUE)
    overlap.weak <- lapply(support.gr, function(ele_) {
      res <- as(findOverlaps(weak, ele_), "List")
      res
    })
    ov.hitIndex <- lapply(overlap.weak, function(ele_) {
      sup <- lapply(support.gr, function(obj) {
        obj.hitIdx <- as(which.min(extractList(obj$pvalue, ele_)), "List")
        obj.hitIdx <- id0[!is.na(obj.hitIdx)]
      })
      sup <- sup[!duplicated(sup)]
      res <- sup
      return(res)
    })
    #
  }
}




##========================================================================
## most updated code (efficient)
##========================================================================
foo <- GRanges(
  seqnames=Rle(c("chr1", "chr2", "chr3", "chr4"), c(5, 6, 4, 3)),
  ranges=IRanges(seq(1, by=9, len=18), seq(6, by=9, len=18)),
  rangeName=letters[seq(1:18)], score=sample(1:25, 18, replace = FALSE))

foo$pvalue <- 10^(score(foo)/-1);

bar <- GRanges(
  seqnames=Rle(c("chr1", "chr2", "chr3","chr4"), c(4, 7, 5, 4)),
  ranges=IRanges(seq(2, by=11, len=20), seq(8, by=11, len=20)),
  rangeName=letters[seq(1:20)], score=sample(1:25, 20, replace = FALSE))
bar$pvalue <- 10^(score(bar)/-1);

moo <- GRanges(
  seqnames=Rle(c("chr1", "chr2", "chr3","chr4"), c(8, 7, 6, 4)),
  ranges=IRanges(seq(4, by=11, len=25), seq(9, by=11, len=25)),
  rangeName=letters[seq(1:25)], score=sample(1:25, 25, replace = FALSE))

moo$pvalue <- 10^(score(moo)/-1);

# classify a by given threshold value:
foo.s <- foo %>% subset(pvalue < 1e-12)
foo.w <- foo %>% subset(pvalue < 1e-4 & pvalue >= 1e-12)

##################################################################
foo <- GRangesList(stringent=foo.s, weak=foo.w)
stringent <- S4Vectors::getListElement(foo, 1, exact = TRUE)
weak <- S4Vectors::getListElement(foo, 2, exact = TRUE)

##================================================================
## subject list
support.gr <- GRangesList(bar,moo) 

##================================================================
## for all stringent peaks, find overlapped regions from bar, moo
ov.s <- lapply(support.gr, function(ele_) {
  res <- as(findOverlaps(foo.s, ele_), "List")
  idx <- as(which.min(extractList(ele_$pvalue, res)), "List")
  idx <- idx[!is.na(idx)]
  res[idx]
})

##=================================================================
## get kardinality
ov.s.len <- Reduce('+', lapply(ov.s, function(ele_) { lengths(ele_)}))
Kard.s <- ov.s.len + lengths(foo.s)

##================================================================
## do initial filtering
keep <- Kard.s >= 2L
foo.s.ok <- foo.s[keep]

# filter out (correct)
ov.s.idx.ok <- lapply(ov.s, function(ele_) {
  ele_[keep]
})

##=====================================================================================
## Solution: TODO fisher.test, first I don't expand bar.ov, moo.ov,
## instead I get its pvalue by ov.s.idx.ok, then put them in matrix, and do fisher.test
v1 <- as(foo.s.ok$pvalue, "List")
v2 <- extractList(bar$pvalue, ov.s.idx.ok[[1]])
v3 <- extractList(moo$pvalue, ov.s.idx.ok[[2]])

# the solution for replacing numeric(0) with 0
v3 <- lapply(v3, function(x) {
  res <- ifelse(length(x)>0, x, 0)
})

##=============================================
## Result is perfect (very nice)
v1 <- unlist(v1)
v2 <- unlist(v2)
v3 <- unlist(v3)

dat <- DataFrame(v1,v2,v3)

##=============================================
## then, TODO Fisher exact test for v1, v2, v3 by element-wise
library(metap)
comb.pval <- suppressWarnings(
  res <- apply(dat[,1:3],1, function(ro) sumlog(ro)$p)
)

## then, based on the combined pvalue, filter out bar.ov.ok, moo.ov.ok, foo.ok
foo.s.ok$comb.p <- comb.pval
bar.expand.s <- bar[unlist(extractList(seq_along(bar), ov.s.idx.ok[[1]]))]
bar.expand.s$comb.pval <- comb.pval

tmp <- data.frame('p'=unlist(v3), 'comb.p'=unlist(comb.pval))
tmp <- tmp[tmp$p != 0e+00,]

moo.expand.s <- moo[unlist(extractList(seq_along(moo), ov.s.idx.ok[[2]]))]
moo.expand.s$comb.p <- tmp[[2]]

##=====================================================================================
# expand overlapped regions from bar.ov, moo.ov respectively (test is successful)
conf <- comb.pval 

##=============================================
## Just experiment
pval.s.list <- List(foo.s.ok$pvalue, bar.expand.s$pvalue, moo.expand.s$pvalue)

##=================================================================
ov.w <- lapply(support.gr, function(ele_) {
  res <- as(findOverlaps(foo.w, ele_), "List")
  idx <- as(which.min(extractList(ele_$pvalue, res)), "List")
  idx <- idx[!is.na(idx)]
  res[idx]
})

ov.w.len <- Reduce('+', lapply(ov.w, function(ele_) { lengths(ele_)}))
Kard.w <- ov.w.len + lengths(foo.w)

keep <- Kard.w >=2L
foo.w.ok <- foo.w[keep]

ov.w.idx.ok <- lapply(ov.w, function(ele_) {
  ele_[keep]
})

#==========================================================
# expand overlapped regions from bar.ov, moo.ov respectively (test is successful)

bar.expand.w <- bar[unlist(extractList(seq_along(bar), ov.w.idx.ok[[1]]))]
moo.expand.w <- moo[unlist(extractList(seq_along(moo), ov.w.idx.ok[[2]]))]

##=========================================================
# chisq.test
pval.w.list <- List(foo.w.ok$pvalue, bar.expand.w$pvalue, moo.expand.w$pvalue)
