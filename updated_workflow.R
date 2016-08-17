## updated complete work flow for integrating hit list of each independent test before procedding initial filtering

aa <- GRanges( seqnames=Rle("chr1", 3),ranges=IRanges(c(2,7,16), c(5,14,20)),
              rangeName=c("a1", "a2", "a3"), score=c(4, 6,9))

bb <- GRanges(seqnames=Rle("chr1", 3),ranges=IRanges(c(4,13,26), c(11,17,28)),
             rangeName=c("b1", "b2", "b3"), score=c(11, 7, 8))

cc <- GRanges(seqnames=Rle("chr1", 4),ranges=IRanges(c(1,4,10, 23), c(3,8,14, 29)),
             rangeName=c("c1", "c2", "c3", "c4"), score= c(4, 6, 3, 8))


overlap <- function(set, idx=1:3, FUN=which.min) {
  chosen <- set[[idx]]
  que.hit <- as(findOverlaps(chosen), "List")
  sup.hit <- lapply(set[-idx], function(ele_) {
    ans <- as(findOverlaps(chosen, ele_), "List")
    out.idx0 <- as(FUN(extractList(ele_$score, ans)), "List")
    out.idx0 <- out.idx0[!is.na(out.idx0)]
    ans <- ans[out.idx0]
  })
  res <- c(list(que.hit),sup.hit)
  names(res) <- c(names(set[idx]),names(set[-idx]))
  return(res)
}

ss <- list(a=aa,b=bb,c=cc) # each replicate must be assigned to named variable.
v1 <- overlap(set = ss, idx = 1L, FUN = which.min)
v2 <- overlap(set = ss, idx = 2L, FUN = which.min)
v3 <- overlap(set = ss, idx = 3L, FUN = which.min)

##========================================================
## optimized
func <- function(myList) {
  # input param checking
  res <- lapply(myList, function(ele_) {
    ans <- drop(ele_)
    ans[is.na(ans)] <- 0L
    ans
  })
  return(res)
}

hit.Intvec <- Map(func, list(v1, v2, v3))
##========================================================
res_all <- lapply(hit.Intvec, as.data.frame)

library(dplyr)
final.hit <- bind_rows(res_all) %>% distinct %>% arrange(a)
res <- as.matrix(final.hit)

final_hit <- apply(res, 2, function(ele_) {
  ele_ <- as(ele_, "IntegerList")
  ele_[all(ele_==0L)] <- IntegerList(integer(0))
  ele_
})

## compute total overlapped number
Reduce('+', lapply(final_hit, lengths))

## now read to do initial filtering
