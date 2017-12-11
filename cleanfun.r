########################################################################################################################
####The Bayesian edge detector
####Function nLimgFun1 is the main function for the edge detection, which include the Bayesian detector and Canny detector. The Canny de####tector is wroten by the cannyEdges function from imager package
####Funcion postprocess is the postprocess routine including the maximum suppression and thresholoding to thin the edges. 
########################################################################################################################
library(spatialfil)
library(imager)
library(rmutil)
library(R2Cuba)
library(mvtnorm)
library(R.matlab)
phimd <- function(theta, theta0, sigmaphi, beta){
 diag(exp(-( (theta - theta0 ) %*% solve(sigmaphi) %*% t(theta - theta0 ) ) ^ ( beta))   )
}
dnomfunmd <- function(theta, theta0, gthetamd, beta, sigmaphi, ...){
    if(beta > 0){
        phitheta0 <- 1
    }else{
         phitheta0 <- 0
     }
    (phitheta0 - phimd(theta, theta0, sigmaphi, beta)) * gthetamd(theta, ...)
}

dnomfunmd1 <- function(theta, theta0, gthetamd, beta, sigmaphi, ...){
    if(beta > 0){
        phitheta0 <- 1
    }else{
         phitheta0 <- 0
     }
    (phitheta0 - phimd(theta, theta0, sigmaphi, beta)) * gthetamd(theta, ...)/dmvnorm(theta, theta0, sigmaphi)
}


gtheta1md <- function(theta, ...){

   a <- ...[[1]]
   b <- ...[[2]]
   
  theta[, 1] < a[2] &  theta[, 1] >= a[1] & theta[, 2] < b[2] & theta[, 2] >= b[1]
   
    
}
gtheta2md <- function(theta, ...){

    gamma <- ...[[1]]
    beta <- ...[[2]]
  diag(exp(- ((theta) %*% solve( gamma) %*% t(theta)) ^ (beta)))
    
}


likfunmd <- function(theta, Sigma,  X){
  dmvnorm(X, theta, Sigma)
  
}
m1funmd<- function(theta, w, theta0, gthetamd,beta, sigmaphi, ...){
   
    dnomfunmd(theta, theta0, gthetamd,beta, sigmaphi, ...) * w 
}


guess.kmeans1 <- function(x)
{
    den<- density(x)
    kmeancut <- den$x[which(den$y == max(den$y))]
    out <- kmeans(as.vector(x),centers=c(kmeancut, quantile(x, 1)))
    list(t1=max(x[out$cluster==1]),t2=min(x[out$cluster==2]))
}
guess.kmeans <- function(x)
{
    out <- kmeans(as.vector(x),centers=c(min(x),mean(x),max(x)))
    list(t1=max(x[out$cluster==1]),t2=max(x[out$cluster==2]))
}




nLimgFun1 <- function(img, sigmaphi, beta, blur, s = 1, gaussian = FALSE){
##Arguments:
    ##img: input image
    ##sigmaphi: the prior sigma parameter, default to be diag(c(1, 1), 2, 2). If the edges are not well connected, choose a larger number, e.g. diag(c(3, 3), 2, 2)
    ##beta: the beta parameter control the number of edge points. The larger |beta| is, the less the edge points are.
    ##blur: The smoothing parameter, which can be the standard deviation for the Guassian filter or the moving average. It correponses to the sigma parameter in isoblur in imager package
    ## s: The sample standard deviation, default to be 1
    ##gaussian: logic; whether to use gaussian filter for smoothing
##Values:
    ##edge: The resulting edge from Bayes detector
    ##mag1: Array valued edge. The same as edge
    ##mag: The raw canny edge result before maximum suppression and thresholding
    ##img: input image
    ##rimg: smoothed image
    ##gr: The gradient used in postpocessing.
    img <- grayscale(img)
    if(blur!= 0){
    im <- grayscale(img) %>% isoblur(blur, gaussian =  gaussian)
}else{
    im <- img
}
    gr <- imgradient(im,"xy")

    mag <- with(gr,sqrt(x^2+y^2))
    slice <- im[, , 1, 1]
    slice <- slice/sd(slice)
    edge <- slice
    vslice <- as.vector(slice)
    nr <- nrow(slice)
    nc <- ncol(slice)
    k <-1
    print(s)
    s <- diag(var(as.vector(slice)) * 2, 2, 2)
    sigmaphi <- sigmaphi
    ss <- solve(s)
    n <- 3
    first <- log((det(2 * pi * s))^(-n/2) * det(2 * pi * s/n)^{1/2})
    taug <- mgauss.hermite(intnum, c(0, 0), sigmaphi)
    tau <- sum(dnomfunmd1(taug$points, c(0, 0), gtheta2md, beta, sigmaphi = sigmaphi, list(diag(1, 2, 2), c(1))) * taug$weights)
  
    for(i in 1:nrow(slice)){
        for(j in 1:ncol(slice)){
            vslicex <- as.vector(c((slice[min(i+ 1, nr), j ] - slice[max(i- 1, 1), j ]), slice[min(i+ 1, nr), min(j+ 1, nc) ] - slice[max(i-1, 1), min(j +1, nc)], slice[min(i+ 1, nr), max(j -1, 1)] - slice[max(i - 1, 1), max(j -1, 1)]))
            vslicey <- c((slice[i, min(j + 1, nc)] - slice[i, max(j - 1, 1)]),  (slice[min(i + 1, nr), min(j + 1, nc)] - slice[min(i+ 1,nr), max(j - 1, 1)]),  (slice[max(i-1, 1), min(j + 1, nc)] - slice[max(i-1, 1), max(j - 1, 1)]))
            n <- length(vslicex)
            x <- cbind(vslicex, vslicey)
   
            mu <- apply(t(x), 1, mean)
         
            mq <- mgauss.hermite(intnum, mu, s/n)
            m1 <- m1funmd(mq$points, w = mq$weights, theta0= c(0, 0), gthetamd = gtheta2md, beta = beta,  sigmaphi = sigmaphi, list(diag(1, 2, 2), c(1)))
           multi <- first +  (n/2 * t(mu) %*%ss%*% mu - n/2 *sum(diag(x %*% ss %*% t(x)))/n)
          
             m0 <- likfunmd(theta = c(0, 0),Sigma = s, X =x)
   
            edge[i, j]<- log(sign(beta) * sum(m1)) - sum(log(m0)) + multi - log(sign(beta) * tau)
        }
       print(i)
    }
    mag1 <- mag
    mag1[, , 1, 1] <- edge
 
    res <- list(edge = edge/3, mag1 = mag1/3, mag = mag, img = img, rimg = im, gr = gr)
    return(res)
}

     

postprocess <- function(img, edge,gr,  cutoff = 1, bayes = 1){
##Arguments:
    ##img: input image
    ##edge: the edge result from nLimgFun1. Must be a cimg object defined in imager package
    ##gr: the gradient output from the nLimgFun1
    ##cutoff: 1  one cutoff, no weak edge; 2 has two cutoff including weak edge and strong edge; The larger values are edges
    ##bayes: whether processing the Bayesian detector result
##Values:
    ##out: the resulting thining edge
    ##t: the cutoff in the thresholding procedure. 
    
    
    mag1 <- edge
    alpha <- 1
    gr <- imgradient(img,"xy")
    norm = sqrt(gr$x^2 + gr$y^2)
    nX <- Xc(img) + gr$x/norm
    nY <- Yc(img) + gr$y/norm
    val.fwd <- interp(edge,data.frame(x=as.vector(nX),y=as.vector(nY)))
    nX <- Xc(img) - gr$x/norm
    nY <- Yc(img) - gr$y/norm
    val.bwd <- interp(edge,data.frame(x=as.vector(nX),y=as.vector(nY)))
    throw <- (edge < val.bwd) | (edge < val.fwd)
    if(bayes == 1){
         den<- density(edge)
         kmeancut <- den$x[which(den$y == max(den$y))]
        edge[throw] <- kmeancut
    }
    else{
        edge[throw] <- 0
    }

    if(cutoff == 1){
        guess <- try(guess.kmeans(as.vector(edge)))
        t1 <-alpha*guess$t1
        t2 <- alpha*guess$t2
        t1 <- t2 <- max(c(t1, t2))
    }else if(cutoff == 2){
        guess <- guess.kmeans(as.vector(edge))
        t1 <-alpha*guess$t1
        t2 <- alpha*guess$t2
    }else{
        t1 <- t2 <- cutoff
    }
    strong <- as.cimg(edge>t2)
    weak <- as.cimg(edge %inr% c(t1,t2))
    out <- rescueFill(strong,weak)
    attr(out,"thresholds") <- c(t1,t2)
    as.pixset(out)
    return(list(out, t1))
}
