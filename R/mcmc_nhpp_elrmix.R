do.mcmc.nhpp.erlmix <- function(nsim,t,maxT,pri.param,init=NULL){

  H0 <- function(J,theta,a.H0,b.H0){
    u <- (1:J)*theta
    l <- c(0,u[-J])
    res <- (u/b.H0)^a.H0-(l/b.H0)^a.H0
    return(res)
  }

  acpt.prob.theta <- function(theta.p,theta.c,gamma,H,c0.H,b.H0){
    H0.H.p <- H0(J,theta.p,a.H0,b.H0)
    H0.H.c <- H0(J,theta.c,a.H0,b.H0)

    res <- sum(c0.H*(H0.H.p-H0.H.c)*log(c0.H)-(lgamma(c0.H*H0.H.p)-lgamma(c0.H*H0.H.c))+c0.H*((H0.H.p-H0.H.c)*log(H)))-
      sum(gamma*(log(theta.p)-log(theta.c))+t*(1/theta.p-1/theta.c))
    numer <- -sum(H*pgamma(maxT,1:J,scale=theta.p))-(a.theta+1)*log(b.theta+theta.p)+log(theta.p)
    denom <- -sum(H*pgamma(maxT,1:J,scale=theta.c))-(a.theta+1)*log(b.theta+theta.c)+log(theta.c)
    return(res+numer-denom)
  }

  acpt.prob.b <- function(b.H0.p,b.H0.c,H,c0.H){
    H0.H.p <- H0(J,theta,a.H0,b.H0.p)
    H0.H.c <- H0(J,theta,a.H0,b.H0.c)

    res <- sum(c0.H*(H0.H.p-H0.H.c)*log(c0.H)-(lgamma(c0.H*H0.H.p)-lgamma(c0.H*H0.H.c))+c0.H*((H0.H.p-H0.H.c)*log(H)))
    numer <- -a.b.H0*b.H0.p+log(b.H0.p)
    denom <- -a.b.H0*b.H0.c+log(b.H0.c)

    return(res+numer-denom)
  }

  acpt.prob.c0 <- function(c0.H.p,c0.H.c,H,b.H0){
    H0.H <- H0(J,theta,a.H0,b.H0)

    res <- sum(H0.H*(c0.H.p*log(c0.H.p)-c0.H.c*log(c0.H.c))-(lgamma(c0.H.p*H0.H)-lgamma(c0.H.c*H0.H))+
                 H0.H*((c0.H.p-c0.H.c)*log(H))-H*(c0.H.p-c0.H.c))
    numer <- -a.c0.H*c0.H.p+log(c0.H.p)
    denom <- -a.c0.H*c0.H.c+log(c0.H.c)

    return(res+numer-denom)
  }

  # prior specification
  J <- pri.param$J
  a.theta <- pri.param$a.theta
  b.theta <- pri.param$b.theta
  sig2.theta <- pri.param$sig2.theta
  a.c0.H <- pri.param$a.c0.H
  sig2.c0.H <- pri.param$sig2.c0.H
  # a.H0 <- pri.param$a.H0
  a.H0 <- 1
  a.b.H0 <- pri.param$a.b.H0
  sig2.b.H0 <- pri.param$sig2.b.H0

  # initialization for MCMC
  n <- length(t)
  if(is.null(init)){
    theta <- b.theta/(a.theta-1)
    c0.H <- 1/a.c0.H
    b.H0 <- 1/a.b.H0
  }else{
    theta <- init$theta
    c0.H <- init$c0.H
    b.H0 <- init$b.H0
  }
  H <- H0.H <- H0(J,theta,a.H0,b.H0)
  gamma <- sample(1:J,n,replace=TRUE,prob=H/sum(H))

  out <- list()
  out$time <- paste("MCMC iteration begins at ",date())
  out$ar.theta <- 0
  out$ar.c0.H <- 0
  out$ar.b.H0 <- 0
  print(out$time)

  for(index in 1:nsim){

    #update of augxiliary variable gamma
    gamma <- rep(0,n)
    for(i in 1:n){
      log.prob <- dgamma(t[i],1:J,scale=theta,log=TRUE)+log(H)
      if(all(exp(log.prob)==0)){
        log.prob <- log.prob-max(log.prob)
      }
      prob <- exp(log.prob)
      gamma[i] <- findInterval(runif(1),cumsum(prob/sum(prob)))+1
    }

    #update of theta
    prop.log.theta <- rnorm(1,log(theta),sqrt(sig2.theta))
    acpt.ratio.theta <- acpt.prob.theta(exp(prop.log.theta),theta,gamma,H,c0.H,b.H0)
    if(log(runif(1))<min(0,acpt.ratio.theta)){
      theta <- exp(prop.log.theta)
      # out$ar.theta <- c(out$ar.theta,1)
      out$ar.theta <- out$ar.theta+1
    }

    #update of c0.H
    prop.log.c0.H <- rnorm(1,log(c0.H),sqrt(sig2.c0.H))
    acpt.ratio.c0.H <- acpt.prob.c0(exp(prop.log.c0.H),c0.H,H,b.H0)
    if(log(runif(1))<min(0,acpt.ratio.c0.H)){
      c0.H <- exp(prop.log.c0.H)
      # out$ar.c0.H <- c(out$ar.c0.H,1)
      out$ar.c0.H <- out$ar.c0.H+1
    }

    #update of b.H0
    prop.log.b.H0 <- rnorm(1,log(b.H0),sqrt(sig2.b.H0))
    acpt.ratio.b.H0 <- acpt.prob.b(exp(prop.log.b.H0),b.H0,H,c0.H)
    if(log(runif(1))<min(0,acpt.ratio.b.H0)){
      b.H0 <- exp(prop.log.b.H0)
      # out$ar.b.H0 <- c(out$ar.b.H0,1)
      out$ar.b.H0 <- out$ar.b.H0+1
    }

    #update of H
    H0.H <- H0(J,theta,a.H0,b.H0)
    H <- rgamma(J,as.numeric(table(c(1:J,gamma))-1)+c0.H*H0.H,pgamma(maxT,1:J,scale=theta)+c0.H)
    H[H==0] <- ((as.numeric(table(c(1:J,gamma))-1)+c0.H*H0.H)/(pgamma(maxT,1:J,scale=theta)+c0.H))[H==0]

    # outputs
    out$H <- rbind(out$H,H)
    out$theta <- c(out$theta,theta)
    out$c0.H <- c(out$c0.H,c0.H)
    out$b.H0 <- c(out$b.H0,b.H0)

    if(index%%(nsim/5)==0){
      print(paste("number of iterations = ",index," ",index/nsim*100,"% MCMC iteration has been completed. ",date()))
      out$time <- rbind(out$time,paste("number of iterations = ",index," ",date()))
    }
  }
  out$ar.theta <- out$ar.theta/nsim
  out$ar.c0.H <- out$ar.c0.H/nsim
  out$ar.b.H0 <- out$ar.b.H0/nsim

  return(out)
}
