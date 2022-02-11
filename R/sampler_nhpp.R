sampler.nhpp <- function(maxT,param,type){

  if(type=="D"){

    # Weibull hazard function: Decreasing intensity function
    a <- param$shape
    b <- param$scale
    if(a>=1) warning("For a decreasing intensity, the shape parameter in the Weibull hazard function must be less than 1.")
    n <- rpois(1,(maxT/b)^a)
    u <- runif(n)
    samp <- u^(1/a)*maxT

  }else if(type=="I"){

    # Weibull hazard function: Increasing intensity function
    a <- param$shape
    b <- param$scale
    if(a<=1) warning("For an increasing intensity, the shape parameter in the Weibull hazard function must be greater than 1.")
    n <- rpois(1,(maxT/b)^a)
    u <- runif(n)
    samp <- u^(1/a)*maxT

  }else if(type=="B"){

    # proportional to mixture of Weibull distributions: Bimodal intensity function
    a <- param$shape
    b <- param$scale
    w <- param$weight
    Lambda <- w[1]*pweibull(maxT,a[1],b[1])+w[2]*pweibull(maxT,a[2],b[2])
    n <- rpois(1,Lambda)
    s <- findInterval(runif(n),c(w[1]*pweibull(maxT,a[1],b[1])/Lambda,1))+1
    u <- runif(n,0,pweibull(maxT,a[1],b[1]))
    u2 <- runif(n,0,pweibull(maxT,a[2],b[2]))
    u[s==2] <- u2[s==2]
    samp <- qweibull(u,a[s],b[s])

  }else{
    warning("The type parameter must be either D, I or B.")
  }

  return(samp)
}
