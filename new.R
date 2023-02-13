library(parallel)
library(foreach)
library(doParallel)
library(moments)
cl <- makeCluster(detectCores() - 2, outfile="log.txt")

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
    se <- sejk <- sol <- avg <- rep(NA, num)
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
      
      mcout <- mcinfo(xmat, lambda, beta, rateC)
      avg[i] <- mcout$avg
      sol[i] <- mcout$sol
    }
    
    # number of times Standard Error < Standard Deviation
    vx <- sqrt(var(t(x)))
    print("vx"); print(vx)
    total1 <- sum(se < vx[1,1])
    total2 <- sum(sejk < vx[1,1])
    count <- sum(sol < vx[1,1])
    
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
                count = count))
  }
  
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
names(result)<-as.character(ratevector)

current_seed <- 1
resultcc <- foreach(rate = ratevector) %dopar% {
  setCoreSeed(current_seed)
  current_seed <- current_seed + 1
  simoutcc <- sim_cc(num=nmc, p=5, rateC=rate, N=sampsize)
  return(simoutcc)
}
plotresults<-function(result){
# finding plot dimensions
   xmn <- 0
   xmx <- 0
   ymn <- 0
   ymx <- 0
   out<-array(NA,c(length(result),4,2))
   dimnames(out)<-list(names(result),
      c("Asymptotic","Jackknife","Expected Informaton","True"),
      c("skewness","kurtosis"))
   for (ii in seq(length(result))) {
      cur <- result[[ii]]
      teststatistics=list((cur$x-cur$beta)/cur$se,
                          (cur$x-cur$beta)/cur$sejk,
                          (cur$x-cur$beta)/cur$sol,
                          (cur$x-cur$beta)/cur$vx[1,1])
      aa <- qqnorm(teststatistics[[1]], plot.it=FALSE)
      bb <- qqnorm(teststatistics[[2]], plot.it=FALSE)
      cc <- qqnorm(teststatistics[[3]], plot.it=FALSE)
      dd <- qqnorm(teststatistics[[4]], plot.it=FALSE)
      xmn <- min(xmn, range(c(aa$x, bb$x, cc$x, dd$x))[1])
      xmx <- max(xmx, range(c(aa$x, bb$x, cc$x, dd$x))[2])
      ymn <- min(ymn, range(c(aa$y-aa$x, bb$y-bb$x, cc$y-cc$x, dd$y-dd$x))[1])
      ymx <- max(ymx, range(c(aa$y-aa$x, bb$y-bb$x, cc$y-cc$x, dd$y-dd$x))[2])
   }

# plots
   for (ii in seq(length(result))) {
      cur <- result[[ii]]
      print(sum(cur$sol < sqrt(var(t(cur$x))[1,1])))
      print(cur$totals)
      teststatistics=list((cur$x-cur$beta)/cur$se,
                          (cur$x-cur$beta)/cur$sejk,
                          (cur$x-cur$beta)/cur$sol,
                          (cur$x-cur$beta)/cur$vx[1,1])
      aa <- qqnorm(teststatistics[[1]], plot.it=FALSE)
      bb <- qqnorm(teststatistics[[2]], plot.it=FALSE)
      cc <- qqnorm(teststatistics[[3]], plot.it=FALSE)
      dd <- qqnorm(teststatistics[[4]], plot.it=FALSE)
 
      plot(c(xmn - 1, xmx + 1), c(ymn - 1, ymx + 1), type = "n",
         xlab = "Theoretical Quantile", ylab = "Deviation from Expected",
         main = paste("Censoring Rate = ", cur$rate))
      points(aa$x, aa$y-aa$x, pch=1, col=1, cex=.6)
      points(bb$x, bb$y-bb$x, pch=2, col=2, cex=.6)
      points(cc$x, cc$y-cc$x, pch=3, col=3, cex=.6)
      points(dd$x, dd$y-dd$x, pch=4, col=4, cex=.6)
      abline(0, 0)
      legend("bottomright", c("Naive", "Jackknife", "Information", "Exact SD"), 
             pch=1:4, col=1:4)
#     browser()
   }
}
stopCluster(cl)

teststats <- function(cur){
      teststatistics<-list((cur$x[1,]-cur$beta)/cur$se,
                          (cur$x[1,]-cur$beta)/cur$sejk,
                          (cur$x[1,]-cur$beta)/cur$sol,
                          (cur$x[1,]-cur$beta)/cur$vx[1,1])
      names(teststatistics) <- c("Naive", "Jackknife", "Information", "Exact SD")
      return(teststatistics)
}
library(moments)
plotresults<-function(result){
# finding plot dimensions
   xmn <- 0
   xmx <- 0
   ymn <- 0
   ymx <- 0
   qqdata<-vector(mode="list",length=length(result))
   for (ii in seq(length(result))) {
      cur <- result[[ii]]
      teststatistics<-teststats(cur)
      if(ii==1){
         out<-array(NA,c(length(result),length(teststatistics),2))
         dimnames(out)<-list(names(result),names(teststatistics),c("skewness","kurtosis"))
      }
      qqdata[[ii]]<-vector(mode="list",length=4)
      for(jj in seq(length(teststatistics))){
#        print(teststatistics[[jj]])
         out[ii,jj,1]<-skewness(teststatistics[[jj]])
         out[ii,jj,2]<-kurtosis(teststatistics[[jj]])
         qqdata[[ii]][[jj]] <- qqnorm(teststatistics[[jj]], plot.it=FALSE)
         xmn <-min(xmn, qqdata[[ii]][[jj]]$x)
         xmx <-max(xmx, qqdata[[ii]][[jj]]$x)
         ymn <-min(ymn, qqdata[[ii]][[jj]]$y-qqdata[[ii]][[jj]]$x)
         ymx <-max(ymx, qqdata[[ii]][[jj]]$y-qqdata[[ii]][[jj]]$x)
         cat("xmn,xmx,ymn,ymx",xmn,xmx,ymn,ymx,"\n")
      }
   }
# plots
#  pdf()
   cat("seq(length(result))",seq(length(result)),"\n")
   for (ii in seq(length(result))) {
      cur <- result[[ii]]
      print(sum(cur$sol < sqrt(var(t(cur$x))[1,1])))
      print(cur$totals)
      png(filename=paste("new",ii,".png",sep=""),height=300)
      par(mar=c(5.1,4.1,1,2.1))
      cat("New plot",ii,"\n")
      plot(c(xmn , xmx ), c(ymn , ymx ), type = "n",
#        main = paste("Censoring Rate = ", cur$rate),
         xlab = "Theoretical Quantile", ylab = "Deviation from Expected"
      )
      for(jj in seq(length(qqdata[[ii]]))){
         points(qqdata[[ii]][[jj]]$x, qqdata[[ii]][[jj]]$y-qqdata[[ii]][[jj]]$x, 
            pch=jj, col=jj, cex=.6)
      }
      abline(0, 0)
      legend("bottomright", names(teststatistics), pch=1:4, col=1:4)
#     browser()
   }
   png(filename="stop.png")
   return(out)
}
out<-plotresults(result)
library(xtable)
print(xtable(out[,,1],digits=4,type="latex"),file="skewness.tex")
print(xtable(out[,,1]-3,digits=4,type="latex"),file="excesskurtosis.tex")

