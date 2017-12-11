source('bnedge.r')
simuimag <- load('simuimag.Rdata')
simutemp <- nLimgFun1(as.cimg(simuimag), diag(1, 2, 2), -1, 1)
temp <- postprocess(simutemp$img, simutemp$mag1, simutemp$gr, 1)[[1]]
mstd <- c(0.2, 0.5, 1, 1.2, 1.5)
lmstd <- length(mstd)
nsim <- 100
simuind <- which(temp == 1, arr.ind = T)
lpsbind <- lpsind <- psbind <- psind <- vector('list')
sm <- 1
sm1 <- 2
for(vind in 1:lmstd){
    std <- mstd[vind]
    for(itr in 1:nsim){
        set.seed(itr)
        newimag <- simuimag
        newimag[simuimag == 1] <- 1
        newimag <- newimag + matrix(rnorm(length(simuimag), 0, std), nrow(simuimag), ncol(simuimag))
        temp <- nLimgFun1(as.cimg(newimag), diag(sm, 2, 2), 2, sm1, s = 1, gaussian = TRUE)
        
        psb <- postprocess((temp$rimg),temp$mag1 + abs(min(temp$mag1)),temp$gr,  cutoff = 2)[[1]][, , 1, 1]
        ps <- postprocess(temp$rimg, (temp$mag),temp$gr,  cutoff = 2, bayes = 0)[[1]][, , 1, 1]
        psbind[[itr]] <- which(psb == 1, arr.ind = T)
        psind[[itr]] <- which(ps == 1, arr.ind = T)

    }
    lpsbind[[vind]] <- psbind
    lpsind[[vind]] <- psind
}
save(lpsbind, lpsind, file = 'lpres')
source(boxplot.r)
## plot(as.cimg(newimag),axes=FALSE)
## plot(as.cimg((psb)),axes=FALSE)
## plot(as.cimg(ps),axes=FALSE)
## plot(density(psb), main = '', col = 1, lwd = 3)
## lines(density(ps), main = '', col = 2, lwd = 3, lty = 2)





