library(hypred)
library(doParallel)

#Map in morgans
snp_pos <- cumsum( seq(0,20,0.5)/100 )
snp_num <- length(snp_pos)
snp_intervals <- diff(snp_pos)
len <- tail(snp_pos, n=1)

randmap <- hypredGenome(num.chr = 1,
                        len.chr = len,
                        num.snp.chr = snp_num)
regmap <- hypredNewMap(randmap, snp_pos)

#Founders
A <- rep(0, snp_num)
B <- rep(1, snp_num)
C <- rep(2, snp_num)
D <- rep(3, snp_num)

####
# Execute crosses with selfing for 1-4 generations
# and independence across generations.

selfgen <- function(gamete1, gamete2) {
  gam1 <- hypredRecombine(regmap, genomeA = gamete1, genomeB = gamete2,
                          mutate = F, block = F)
  gam2 <- hypredRecombine(regmap, genomeA = gamete1, genomeB = gamete2,
                          mutate = F, block = F)
  return(rbind(gam1, gam2))
}


# After lots of playing around, it seems that for loops
# are as fast as other approaches for hypred.

twoway <- function(gens=4){
  gamete_1 <- A
  gamete_2 <- B
  
  for(i in 1:gens){
    newgen <- selfgen(gamete_1, gamete_2)
    gamete_1 <- newgen[1,]
    gamete_2 <- newgen[2,]
  }
  return(c(gamete_1, gamete_2))
}

fourway <- function(gens=4){
  gamete_1 <- hypredRecombine(regmap, genomeA = A, genomeB = B,
                        mutate = F, block = F)
  gamete_2 <- hypredRecombine(regmap, genomeA = C, genomeB = D,
                        mutate = F, block = F)
  
  for(i in 1:gens){
    newgen <- selfgen(gamete_1, gamete_2)
    gamete_1 <- newgen[1,]
    gamete_2 <- newgen[2,]
  }
  return(c(gamete_1, gamete_2))
}
  
nsims <- 100000
generations <- 4
nproc <- 4
cl = makeCluster(nproc)
registerDoParallel(cl)

#Two-way
biparental <- foreach(i=1:generations, 
                      .packages="hypred") %dopar% {
                res <- matrix(nrow=nsims, ncol=2*snp_num)
                for(j in 1:nsims) {
                  res[j,] <- twoway(gens = i)
                }
                res
              }

quad <- foreach(i=1:generations, 
                .packages="hypred") %dopar% {
                  res <- matrix(nrow=nsims, ncol=2*snp_num)
                  for(j in 1:nsims) {
                    res[j,] <- fourway(gens = i)
                  }
                  res
                }
for(i in 1:generations) {
  biname <- paste("../data/biparental", snp_num, i, ".txt", sep="_")
  fourname <- paste("../data/fourway", snp_num, i, ".txt", sep="_")
  saveRDS(biparental[[i]], file=paste(biname, '.rds', sep=''))
  saveRDS(quad[[i]], file=paste(fourname, '.rds', sep=''))
  write(t(biparental[[i]]), file=biname, ncolumns = 2*snp_num)
  write(t(quad[[i]]), file=fourname, ncolumns = 2*snp_num)
}
