
####generate plot for figure 7 and 8
nohail <-load.image("tryimage/real72.JPG")
realimage  <- nLimgFun1(nohail, diag(c(1,1), 2, 2), -2, blur = 6,  s= 1, TRUE)
psres <- postprocess(realimage$rimg,realimage$mag1 + abs(min(realimage$mag1)), realimage$gr ,cutoff= 2,  bayes = 1)
psrescanny <- postprocess(realimage$rimg,realimage$mag, realimage$gr, cutoff =  2, bayes = 0)
plot(psres[[1]] + realimage$img, interp=TRUE, main = expression(beta == 2 ~~ sigma[s] == 6), axes=FALSE)
plot(psrescanny[[1]] + realimage$img, interp=TRUE, main = expression(sigma[s] == 6), axes=FALSE)
plot(as.cimg(rescanny))

####generate edge result in Figure 9 ##
nohail <-load.image("tryimage/roof-hail-big1.JPG")
realimage  <- nLimgFun1(nohail, diag(c(3,3), 2, 2), 2, blur = 2,  s= 1, TRUE)
psres <- postprocess(realimage$rimg,realimage$mag1 + abs(min(realimage$mag1)), realimage$gr ,cutoff= 2,  bayes = 1)
psrescanny <- postprocess(realimage$rimg,realimage$mag, realimage$gr, cutoff =  2, bayes = 0)
plot(psres[[1]] + realimage$img, interp=TRUE, main = expression(beta == 2 ~~ sigma[s] == 2), axes=FALSE)
plot(psrescanny[[1]] + realimage$img, interp=TRUE, main = expression(sigma[s] == 2), axes=FALSE)
plot(as.cimg(rescanny))

####generate edge result in Figure 10 ##
nohail <-load.image("tryimage/roof-hail-big2.JPG")
realimage  <- nLimgFun1(nohail, diag(c(1,1), 2, 2), 2, blur = 2,  s= 1, TRUE)
psres <- postprocess(realimage$rimg,realimage$mag1 + abs(min(realimage$mag1)), realimage$gr ,cutoff= 2,  bayes = 1)
psrescanny <- postprocess(realimage$rimg,realimage$mag, realimage$gr, cutoff =  2, bayes = 0)
plot(psres[[1]] + realimage$img, interp=TRUE, main = expression(beta == 2 ~~ sigma[s] == 2), axes=FALSE)
plot(psrescanny[[1]] + realimage$img, interp=TRUE, main = expression(sigma[s] == 2), axes=FALSE)
plot(as.cimg(rescanny))
