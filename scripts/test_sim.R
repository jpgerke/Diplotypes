rm(list=ls())

b1 <- readRDS("../data/biparental_41_1_.txt.rds")
b4 <- readRDS("../data/biparental_41_4_.txt.rds")

c1 <- readRDS("../data/fourway_41_1_.txt.rds")
c4 <- readRDS("../data/fourway_41_4_.txt.rds")

AAfreq <- function(m) {
    n <- nrow(m)
    k <- ncol(m)/2
    dip = c()
    for(i in 1:(k-1)) {
      a1 <- i
      a2 <- i+1
      a3 <- k+i
      a4 <- k+i+1
      AA1 <- apply(m[,a1:a2] == c(0,0), 1, function(x){sum(x)==2})
      AA2 <- apply(m[,a3:a4] == c(0,0), 1, function(x){sum(x)==2})
      dip <- c(dip, sum(AA1 & AA2)/n)
    }
    return(dip)
}

AAfreq(b1)
AAfreq(b4)
AAfreq(c1)
AAfreq(c4)
