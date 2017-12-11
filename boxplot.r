mLmaxdis <- Lmaxdis <- vector('list')
for(vind in 1:lmstd){
    psbind <- lpsbind[[vind]]
    for(itr in 1:nsim){
      nrp <- nrow(psbind[[itr]])
      if(nrp >0 ){
        maxdis <- rep(NA, nrp)

        for(i in 1:nrp){
        
            mindis <- sqrt((psbind[[itr]][i, 1] - simuind[, 1])^2 +  (psbind[[itr]][i, 2] - simuind[, 2])^2) 
            maxdis[i] <- sum(mindis <= sqrt(3^2 + 3^2)) >0# mindis[which(mindis == min(mindis))][1]
        }
        Lmaxdis[[itr]]<- mean(maxdis)#max(maxdis, 0.0)
    }
  }
    mLmaxdis[[vind]] <- unlist(Lmaxdis)
 }

mL1maxdis <- L1maxdis <- vector('list')
for(vind in 1:lmstd){
    psbind <- lpsbind[[vind]]
for(itr in 1:nsim){
nrp <- nrow(simuind)
maxdis <- rep(NA, nrp)
  # ix <- which(psbind[[itr]][, 1] > 2 &psbind[[itr]][, 1] < 52 &psbind[[itr]][, 2] > 18 &  psbind[[itr]][, 2] < 65 )
   #     psbind[[itr]] <- psbind[[itr]][ix, ]
      if(nrow(psbind[[itr]]) >0 ){

for(i in 1:nrp){
    mindis <- sqrt((psbind[[itr]][, 1] - simuind[i, 1])^2 +  (psbind[[itr]][, 2] - simuind[i, 2])^2)
   maxdis[i] <- sum(mindis <= sqrt(3^2 + 3^2)) >0# mindis[which(mindis == min(mindis))][1]
}
L1maxdis[[itr]]<- mean(maxdis)#max(maxdis)
}
}
    mL1maxdis[[vind]] <- unlist(L1maxdis)
}



mLpsmaxdis <- Lpsmaxdis <- vector('list')
for(vind in 1:lmstd){
    psind <- lpsind[[vind]]
    for(itr in 1:nsim){
        nrp <- nrow(psind[[itr]])
        maxdis <- rep(NA, nrp)

for(i in 1:nrp){
    mindis <- sqrt((psind[[itr]][i, 1] - simuind[, 1])^2 +  (psind[[itr]][i, 2] - simuind[, 2])^2)
   maxdis[i] <- sum(mindis <= sqrt(3^2 + 3^2)) >0# mindis[which(mindis == min(mindis))][1]
}
Lpsmaxdis[[itr]]<-mean(maxdis) #max(maxdis)
}
    mLpsmaxdis[[vind]] <- unlist(Lpsmaxdis)
}


mL1psmaxdis <- L1psmaxdis <- vector('list')
for(vind in 1:lmstd){
    psind <- lpsind[[vind]]
    for(itr in 1:nsim){
    #     ix <- which(psind[[itr]][, 1] > 2 &psind[[itr]][, 1] < 52 &psind[[itr]][, 2] > 18 &  psind[[itr]][, 2] < 65 )
    #    psind[[itr]] <- psind[[itr]][ix, ]
         
nrp <- nrow(simuind)
maxdis <- rep(NA, nrp)
for(i in 1:nrp){
    mindis <- sqrt((psind[[itr]][, 1] - simuind[i, 1])^2 +  (psind[[itr]][, 2] - simuind[i, 2])^2)
   maxdis[i] <- sum(mindis <= sqrt(3^2 + 3^2)) >0# mindis[which(mindis == min(mindis))][1]
}
L1psmaxdis[[itr]]<- mean(maxdis)#max(maxdis)
}
    mL1psmaxdis[[vind]] <- unlist(L1psmaxdis)
}

#mLmaxdis <- do.call(rbind, mLmaxdis)
#mL1maxdis <- do.call(rbind, mL1maxdis)
#mLpsmaxdis <- do.call(rbind, mLpsmaxdis)
#mL1psmaxdis <- do.call(rbind, mL1psmaxdis)
hist(sqrt(mLmaxdis[1, ]))


plot(mLmaxdis[1, ] )
for(vind in 1:lmstd){
    filename <- paste('fdr', vind, '.pdf', sep = "")
    pdf(filename)
    boxplot(mLmaxdis[[vind ]], mLpsmaxdis[[vind]], main = paste('SD', '=', mstd[vind], sep = ""),  outline = FALSE, cex = 5, cex.lab = 2, names = c("Bayes", "Canny"), cex.axis = 1.5, lwd = 2, ylim = c(0, 1.2))
    dev.off()
}

for(vind in 1:lmstd){
    filename <- paste('fndr', vind, '.pdf', sep = "")
    pdf(filename)
    boxplot(mL1maxdis[[vind ]], mL1psmaxdis[[vind]], main = paste('SD', '=', mstd[vind], sep = ""),  outline = FALSE, cex = 5, cex.lab = 2, names = c("Bayes", "Canny"), cex.axis = 1.5, lwd = 2, ylim = c(0, 1.2))
    dev.off()
}
