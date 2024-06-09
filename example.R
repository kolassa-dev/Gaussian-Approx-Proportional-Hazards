library(survival)
(out<-coxph(Surv(time,status)~x,data=aml,x=TRUE))
library(PHInfiniteEstimates)#For bestbeta
bestbeta(out,usecc=TRUE)
myfit<-function(x,mydata) coxph(Surv(time,status)~x,
                                data=mydata[x,])$coefficients 
library(bootstrap)#For jackknife
jackknife(seq(dim(aml)[1]),myfit,aml)