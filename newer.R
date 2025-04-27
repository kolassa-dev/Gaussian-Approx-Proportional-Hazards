library(parallel)
library(foreach)
library(doParallel)
library(moments)
cl <- makeCluster(detectCores() - 2, outfile="log.txt")
teststats <- function(cur){
      teststatistics<-vector(length(cur$selist),mode="list")
      for(jj in seq(length(cur$selist))){
         teststatistics[[jj]] <- (cur$x[1,]-cur$beta)/cur$selist[[jj]]
      }
      names(teststatistics) <- names(cur$selist)
      return(teststatistics)
}

clusterEvalQ(cl, {
  library("survival")
  library("bootstrap")
  library("PHInfiniteEstimates")
  library("mlogit")
  library("dfidx")
  
  setCoreSeed <- function(k) {
    set.seed(k)
  }
  
  mcinfo <- function(xmat, lambda, beta, rateC, innersize=100) {
    p <- dim(xmat)[2]
    N <- dim(xmat)[1]
    mysum <- 0
    fishinfo <- array(0, c(p, p))
    for (j in 1:innersize) {
      # see sim for comments
      #     v <- runif(n = N)
      #     Tlat <- (-log(v) / (lambda * exp(xmat %*% rep(beta, p))))
      
      #     C <- rexp(n = N, rate = rateC)
      #     C <- C + log(2) / lambda
      
      #     dat <- list(time=pmin(Tlat, C),status= as.numeric(Tlat <= C))
      dat <- sim(N, p, lambda, beta, rateC, x=xmat)
      
      # update sum of E(l'')
      mysum <- mysum + pllk(rep(beta, p), xmat[order(dat$time),], dat$status[order(dat$time)], cc=TRUE)$d2
      # update sum of Fisher information
      fishinfo <- fishinfo - pllk(rep(beta, p), xmat[order(dat$time),], dat$status, cc=TRUE)$d2
    }
    return(list(avg=mysum[1,1]/innersize, sol=sqrt(solve(fishinfo/innersize)[1,1])))
  }
  
  sim <- function(N, p, lambda, beta, rateC, binary=FALSE, x=NULL) {
    if(is.null(x)){
      if (binary) {
        # (one) binary covariate
        # x <- sample(c(0, 1), N, replace = TRUE, prob = c(0.5, 0.5))
        # covariates --> N binary, p parameters
        x <- array(sample(c(0, 1), N * p, replace=TRUE, prob=c(0.5, 0.5)), c(N, p))
      } else {
        # (one) covariate --> N Normal
        # x <- rnorm(n = N)
        # covariates --> N Normal, p parameters
        x <- array(rnorm(n = N * p), c(N, p))
      }
    }#end if is.null(x)
    
    # Weibull latent event times
    v <- runif(n = N)
    Tlat <- (-log(v) / (lambda * exp(x %*% rep(beta, p))))
    
    # censoring times
    C <- rexp(n=N, rate=rateC)
    # censoring mechanism
    C <- C + log(2) / lambda
    
    # follow-up times and event indicators
    time <- pmin(Tlat, C)
    status <- as.numeric(Tlat <= C)
    
    # data set
    data.frame(id=1:N, time=time, status=status, x=x)
  }
  
  sim_trials <- function(num=1000, p=1, lambda=0.01, beta=-0.8, rateC=0.003, N=30) {
    se <- sejk <- sol <- avg <- sebo <- rep(NA, num)
    x <- array(NA, c(p, num))
    
    theta <- function(ind, dat) {
      f <- coxph(Surv(time, status) ~ x.1 + x.2 + x.3 + x.4 + x.5, data=dat, subset=ind, x=TRUE)
      return(f$coef[1])
    }
    
    for (i in 1:num) {
      if (i %% 10 == 0) {
        cat("sim ")
        cat(i)
        cat("\n")
      }
      dat <- sim(N, p, lambda, beta, rateC)
      fit <- coxph(Surv(time, status) ~ x.1 + x.2 + x.3 + x.4 + x.5, data=dat, x=TRUE)
      coeff <- summary(fit)$coef
      xmat <- cbind(dat$x.1, dat$x.2, dat$x.3, dat$x.4, dat$x.5)
      obsllk <- pllk(coeff[,1], xmat[order(dat$time),], dat$status[order(dat$time)], cc = TRUE)$d2
      x[,i] <- coeff[,1]
      se[i] <- coeff[1,3]
      sejk[i] <- jackknife(1:N, theta, dat)$jack.se
      sebo[i] <- sd(bootstrap(1:N,999,theta, dat)$thetastar)
      mcout <- mcinfo(xmat, lambda, beta, rateC)
      avg[i] <- mcout$avg
      sol[i] <- mcout$sol
    }
    
    # number of times Standard Error < Standard Deviation
    vx <- sqrt(var(t(x))[1,1])
    total1 <- sum(se < vx)
    total2 <- sum(sejk < vx)
    total3 <- sum(sol < vx)
    selist <- list(se=se,sejk=sejk,sebo=sebo,sol=sol,vx=vx)
    names(selist) <- c("Naive", "Jackknife", "Information", "True")
    return(list(totals = c(raw = total1, jack = total2),
                num = num, p = p, lambda= lambda, beta= beta, rateC = rateC, N = N,	
                rawsummary = summary(se),
                jacksummary = summary(sejk), 
                true = mean(vx), 
                avg = avg, 
                vx = vx,
                x = x,
                se = se,
                sejk = sejk,
                sol = sol,
                selist = selist,
                total3 = total3))
  }#End simtrials
  
  sim_cc <- function(num=3000, p=1, lambda=0.01, beta=0, rateC=0.003, N=30) {
    pcount <- array(0, c(2, 3))
    dimnames(pcount) <- list(cc=c("Yes","No"), lev=c("0.01","0.05","0.10"))
    
    theta <- function(ind, dat) {
      f <- coxph(Surv(time, status) ~ x.1 + x.2 + x.3 + x.4 + x.5, data=dat, subset=ind, x=TRUE)
      return(f$coef[1])
    }
    
    for (i in 1:num) {
      if (i %% 10 == 0) {
        cat("cc ")
        cat(i)
        cat("\n")
      }
      dat <- sim(N, p, lambda, beta, rateC)
      fit <- coxph(Surv(time, status) ~ x.1 + x.2 + x.3 + x.4 + x.5, data=dat, x=TRUE)
      firstfit <- bestbeta(fit, usecc = FALSE)
      p2 <- summary(firstfit)$coefficients[1, 5]
      secondfit <- bestbeta(fit, usecc = -sign(firstfit$coef[1])/2, start=firstfit$coef)
      p1 <- summary(secondfit)$coefficients[1, 5]
      if (p1 < 0.01) pcount[1, 1] <- pcount[1, 1] + 1
      if (p1 < 0.05) pcount[1, 2] <- pcount[1, 2] + 1
      if (p1 < 0.10) pcount[1, 3] <- pcount[1, 3] + 1
      if (p2 < 0.01) pcount[2, 1] <- pcount[2, 1] + 1
      if (p2 < 0.05) pcount[2, 2] <- pcount[2, 2] + 1
      if (p2 < 0.10) pcount[2, 3] <- pcount[2, 3] + 1
    }
    
    return(pcount = pcount)
  }
})

registerDoParallel(cl)

nmc <- 1000
sampsize <- 100
ratevector <- c(0.001, 0.005, 0.015, 0.250)

current_seed <- 1
result <- foreach(rate = ratevector) %dopar% {
  setCoreSeed(current_seed)
  current_seed <- current_seed + 1
  simout <- sim_trials(num=nmc, p=5, rateC=rate, N=sampsize)
  return(simout)
}
print("result") ; print(result)
names(result)<-as.character(ratevector)

current_seed <- 1
resultcc <- foreach(rate = ratevector) %dopar% {
  setCoreSeed(current_seed)
  current_seed <- current_seed + 1
  simoutcc <- sim_cc(num=nmc, p=5, rateC=rate, N=sampsize)
  return(simoutcc)
}
plotresults <- function(result){
# finding plot dimensions
   xmn <- 0
   xmx <- 0
   ymn <- 0
   ymx <- 0
   print("results"); print(result)
   for (ii in seq(length(result))) {
      cur <- result[[ii]]
      teststatistics<-teststats(cur)
      if(ii==1){
         plotdata <- array(NA,c(length(result),length(teststatistics),2))
         dimnames(plotdata)<-list(names(result),names(teststatistics), c("skewness","kurtosis"))
      }
      for(jj in seq(length(teststatistics))){
         aa <- qqnorm(teststatistics[[jj]], plot.it=FALSE)
         xmn <- min(xmn, aa$x)
         xmx <- max(xmx, aa$y)
         ymn <- min(ymn, aa$y-aa$x)
         ymx <- max(ymx, aa$y-aa$x)
     }
   }

# plots
   for (ii in seq(length(result))) {
      cur <- result[[ii]]
      print(sum(cur$sol < sqrt(var(t(cur$x))[1,1])))
      print(cur$totals)
      teststatistics<-teststats(cur)
      plot(c(xmn - 1, xmx + 1), c(ymn - 1, ymx + 1), type = "n",
         xlab = "Theoretical Quantile", ylab = "Deviation from Expected",
         main = paste("Censoring Rate = ", cur$rate))
      for(jj in seq(length(teststatistics))){
         points(aa$x, aa$y-aa$x, pch=jj, col=jj, cex=.6)
      }
      abline(0, 0)
      nplot <- length(dimnames(plotdata)[[2]])
      legend("bottomright", 
#         c("Naive", "Jackknife", "Information", "Exact SD"), 
             dimnames(plotdata)[[2]],
             pch=seq(nplot), col=seq(nplot))
#     browser()
   }
   png(filename="stop.png")
   return(plotdata)
}
stopCluster(cl)

library(moments)
momentmat <- plotresults(result)
library(xtable)
print(xtable(momentmat[,,1],digits=4,type="latex"),file="skewness.tex")
print(xtable(momentmat[,,1]-3,digits=4,type="latex"),file="excesskurtosis.tex")
