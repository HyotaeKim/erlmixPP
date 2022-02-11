
# erlmixPP

<!-- badges: start -->
<!-- badges: end -->

Erlang mixture modeling for Poisson process intensities.

## Installation

You can install the package with the following R-commnand:

``` r
remotes::install_github("HyotaeKim/erlmixPP")
```

## Example

``` r
library(erlmixPP)

# Decreasing example given in Section 3.1 of Kim and Kottas (2022)
# synthetic data generation
seed <- 61
set.seed(seed)
maxT <- 2e1
param <- NULL
param$shape <- 0.5
param$scale <- 8e-5
samp <- sampler.nhpp(maxT,param,"D")

# fitting samples to the erlang mixture model
nsim <- 1e3
nsim <- 1e5
nsim <- 1e7

pri.param <- NULL
pri.param$J <- 50
pri.param$a.theta <- 2; pri.param$b.theta <- 1; pri.param$sig2.theta <- 0.0009
pri.param$a.c0.H <- 1/10; pri.param$sig2.c0.H <- 0.15
pri.param$a.b.H0 <- length(samp)/maxT; pri.param$sig2.b.H0 <- 0.1

system.time(post.samp <- do.mcmc.nhpp.erlmix(nsim,samp,maxT,pri.param))
print(c(post.samp$ar.theta,post.samp$ar.c0.H,post.samp$ar.b.H0)) # acceptance ratios of each parameter

# trace plots for the model parameters
par(mfrow=c(2,2))
ts.plot(post.samp$theta)
ts.plot(post.samp$c0.H)
ts.plot(post.samp$b.H0)

# evaluating prior and posterior intensities
library(bayesmeta)
post.samp$pri.theta <- post.samp$pri.H <- NULL
for(i in 1:length(post.samp$theta)){
  pri.theta <- rlomax(1,pri.param$a.theta,pri.param$b.theta)
  pri.c0.H <- rexp(1,pri.param$a.c0.H)
  pri.b.H0 <- rexp(1,pri.param$a.b.H0)
  pri.H0.H <- pri.theta/pri.b.H0
  pri.H <- rgamma(pri.param$J,pri.c0.H*pri.H0.H,pri.c0.H)
  post.samp$pri.theta <- c(post.samp$pri.theta,pri.theta)
  post.samp$pri.H <- rbind(post.samp$pri.H,pri.H)
}

lambda <- function(grid,pri.param,post.samp){
  res <- NULL
  for(i in grid){
    pri.lambda <- post.lambda <- NULL
    for(k in 1:nrow(post.samp$H)){
      pri.lambda <- c(pri.lambda,sum(post.samp$pri.H[k,]*dgamma(i,1:pri.param$J,scale=post.samp$pri.theta[k])))
      post.lambda <- c(post.lambda,sum(post.samp$H[k,]*dgamma(i,1:pri.param$J,scale=post.samp$theta[k])))
    }
    res$pri <- cbind(res$pri,pri.lambda)
    res$post <- cbind(res$post,post.lambda)
  }
  return(res)
}

grid <- seq(0,maxT,1e-1)
grid <- grid[-1]
insty <- lambda(grid,pri.param,post.samp)
pri.est <- post.est <- NULL
pri.est$lbd <- apply(insty$pri,2,quantile,probs=0.025)
pri.est$ubd <- apply(insty$pri,2,quantile,probs=0.975)
pri.est$mean <- apply(insty$pri,2,mean)
post.est$lbd <- apply(insty$post,2,quantile,probs=0.025)
post.est$ubd <- apply(insty$post,2,quantile,probs=0.975)
post.est$mean <- apply(insty$post,2,mean)

# summary plot
pdf(".\\figs\\dec_estimates.pdf",family="Helvetica")
par(mar=c(5.1,4.6,4.1,2.1))
layout(matrix(c(1,1,5,1,1,5,2,3,4),ncol=3,byrow=FALSE))
plot(grid,param$shape/param$scale*(grid/param$scale)^(param$shape-1),typ="l",lwd=2,
     ylim=c(0,150),ylab="Intensity",xlab="Time",bty="n",cex.lab=1.5,cex.main=1.5)
polygon(c(grid,rev(grid)),c(post.est$lbd,rev(post.est$ubd)),col="gray80",lty=0)
lines(grid,post.est$mean,lty=4,lwd=2,col="gray40")
lines(grid,param$shape/param$scale*(grid/param$scale)^(param$shape-1),lwd=2)
hist(post.samp$theta,probability=TRUE,xlab=expression(theta),main="",cex.lab=1.5,cex.main=1.5)
lines(seq(0,10,0.001),dlomax(seq(0,10,0.001),pri.param$a.theta,pri.param$b.theta),lty=2,lwd=2,col="blue")
hist(post.samp$c0.H,probability=TRUE,xlab=expression(c[0]),main="",cex.lab=1.5,cex.main=1.5)
lines(seq(0,1,0.001),dexp(seq(0,1,0.001),pri.param$a.c0.H),lty=2,lwd=2,col="blue")
hist(post.samp$b,probability=TRUE,xlab=expression(b),main="",cex.lab=1.5,cex.main=1.5)
lines(seq(0,1,0.001),dexp(seq(0,1,0.001),pri.param$a.b.H0),lty=2,lwd=2,col="blue")
plot(samp,rep(0,length(samp)),pch="|",xlim=c(0,20),ylab="",xlab="Observed point pattern",yaxt="n",bty="n",cex.lab=1.5,cex.main=1.5)
dev.off()
```
