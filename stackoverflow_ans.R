## facilitate the output represeantation of MSPC Packages generated:

# Question:  How to group set of data.frame objects in nested list with different order?
# example
res_1 <- list(con=list(a.con_1=airquality[1:4,], b.con_1=iris[2:5,], c.con_1=ChickWeight[3:7,]),
              dis=list(a.dis_1=airquality[5:7,], b.dis_1=iris[8:11,], c.dis_1=ChickWeight[12:17,]))

res_2 <- list(con=list(b.con_2=iris[7:11,], a.con_2=airquality[4:9,], c.con_2=ChickWeight[2:8,]),
              dis=list(b.dis_2=iris[2:5,], a.dis_2=airquality[1:3,], c.dis_2=ChickWeight[12:15,]))

res_3 <- list(con=list(c.con_3=ChickWeight[10:15,], a.con_3=airquality[2:9,], b.con_3=iris[12:19,]),
              dis=list(c.dis_3=ChickWeight[2:7,], a.dis_3=airquality[13:16,], b.dis_3=iris[2:7,]))
              
# first solution (very efficient):
myList <- list(res_1, res_2, res_3)
# make a copy of the list to preserve the structure for the new list
myList2 <- myList

for(i in seq_len(length(myList))) {
  # get ordering of inner list names
  myOrder <- rank(names(myList[[c(i,2)]]))

  for(j in seq_len(length(myList[[i]]))) {
    for(k in seq_len(length(myList[[c(i, j)]]))) {
      # reorder content
      myList2[[c(myOrder[k], j, i)]] <- myList[[c(i, j, k)]]
      # rename element
      names(myList2[[c(myOrder[k], j)]])[i] <- names(myList[[c(i, j)]])[k]
    }
  }
}

# second version of solution from different author:
newOrder <- c(2,1,3,4)
exchangeNestedDF <- function(ListOfListsOfDF, newOrder){
  NewList <- lapply(ListOfListsOfDF, function(x) x[newOrder])
  names(NewList) <- names(ListOfListsOfDF)
  return(NewList)
}
list2_NeW <- exchangeNestedDF(list2, newOrder)

NewOrders <- list(c(2,1,3,4),  c(2,3,1,4),  c(2,3,4,1))
LIST <- lapply(1:3, function(i) exchangeNestedDF(get(paste0("list", i)), NewOrders[[i]]))



###===================================================================================================================
# Question: sorting set of data.frame by paired-row?

# example:
foo <- data.frame( start=seq(1, by=4, len=6), stop=seq(3, by=4, len=6))
bar <- data.frame(start=seq(5, by=2, len=7), stop=seq(7, by=2, len=7))
bleh <- data.frame(start=seq(1, by=5, len=5), stop=seq(3, by=5, len=5))

# solution (very efficient) :
out <- rbind.data.frame(foo, bar, bleh)
out <- out[!duplicated(out),]
out1 <- out[do.call(order, out),]
rownames(out1) <- NULL


###=====================================================================================================================
### Question : How to get relative complement of one data.frame in another?

# example: 
foo <- data.frame( start=seq(1, by=4, len=6), stop=seq(3, by=4, len=6))
bleh <- data.frame(start=seq(1, by=5, len=5), stop=seq(3, by=5, len=5))

# solution:
library(dplyr)
dplyr:: setdiff(foo,bleh)


###======================================================================================================================
# Question: How to flatten out nested list into one list more efficiently instead of using unlist method?

# example: 
mylist <- list(pass=list(Alpha.df1_yes=airquality[2:4,], Alpha.df2_yes=airquality[3:6,],Alpha.df3_yes=airquality[2:5,],Alpha.df4_yes=airquality[7:9,]),
             fail=list(Alpha.df1_no=airquality[5:7,], Alpha.df2_no=airquality[8:10,],  Alpha.df3_no=airquality[13:16,],Alpha.df4_no=airquality[11:13,]))

# 1st solution (very efficient, perfect):
res <- lapply(mylist, function(i){
  x <- do.call(rbind, i)
  x[ !duplicated(x), ]
  rownames(x) <- NULL
  x
})

# 2nd solution (it is OKAY, but not desired): 
res <- do.call(rbind, unlist(mylist, recursive = FALSE))
res <- res[!duplicated(res), ]
rownames(res) <- NULL

###===================================================================================================================================================
## Question : perform set purification test (I come up solution)

# example

mylist <- list(pass=list(Alpha.df1_yes=airquality[2:4,], Alpha.df2_yes=airquality[3:6,],Alpha.df3_yes=airquality[2:5,],Alpha.df4_yes=airquality[7:9,]),
               fail=list(Alpha.df1_no=airquality[5:7,], Alpha.df2_no=airquality[8:10,],  Alpha.df3_no=airquality[13:16,],Alpha.df4_no=airquality[11:13,]))


res <- lapply(mylist, function(i){
  x <- do.call(rbind, i)
  x[ !duplicated(x), ]
})

# solution:
library(dplyr)
dplyr::setdiff(res[[1]], res[[2]])

###===============================================================================================================================================

## remove repeated overlap hit-index and only getting upper triangle of grid

myList <- list(foo, bar, bleh)
nn <- length(myList)
index <- expand.grid("query"=1:nn, "target"=1:nn)
index <- index[index[,1] <= index[,2],]

res <- Map(function(i, j) findOverlaps(mylist[[i]], mylist[[j]]), idx[,1], idx[,2] )

## combination for group of index for set of GRanges objects
i1 <- seq_along(myList)
l1 <- lapply(i1, function(i) seq(i))
l2 <- lapply(seq_along(i1), function(i) i1[!i1 %in% sequence(i-1)] )
Map(function(x,y) {x1 <- expand.grid(x,y)
                     x1[c(TRUE, diff(x1[,1])>=0),]}, l1, l2)



###==============================================================================================================================
## Idea is okay but FIXME
myList <- list(
  a=IntegerList(1,2,3,4),
  b=IntegerList(3,4,integer(0),integer(0)),
  c=IntegerList(1,4,6,7)
)

desired_output <- list(
  a=c(1,2,3,4),
  b=c(3,4,0, 0),
  c(1,4,6,7)
)


tmp <- lapply(myList, drop)
tmp <- tmp[!is.na(tmp)]  # fix me

lapply(desired_output, function(ele_) {
  ans <- as(ele_, "IntegerList")
})
