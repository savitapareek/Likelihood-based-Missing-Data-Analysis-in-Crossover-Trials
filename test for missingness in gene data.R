rm(list=ls())
dat=read.csv(file.choose())
library(naniar)
library(misty)

mcar_test(dat)

dat$Response
attach(dat)

## data frame
respos=dat$Response
m=10;p=3;s=3;t=3
no=c(6,6,5)
## fit model with original or log responses
treatment1=c(rep(rep(c("trt2","trt1","trt3"),each=m),no[1]),
             rep(rep(c("trt3","trt2","trt1"),each=m),no[2]),
             rep(rep(c("trt1","trt3","trt2"),each=m),no[3]))
subject=rep(paste("sub",1:sum(no)),each=(p*m))#+subject*gene
Gene=rep(paste("gene",1:m),p*sum(no))
Period=rep(rep(paste("per",1:p),each=m),sum(no))
dat11=data.frame(respos,Period,Gene,treatment1,subject)
mcar_test(dat11)
na.test(dat11)
#na.test(airquality)
# rr=rep(0,length(Response))
# data("airquality")
# 
# rr[which(is.na(Response)==T)]=1
# cbind(Response,rr)
# t.test(Response,rr,paired = F)
# glm(rr~Response)
# mc_resp=data.frame(Response,rr)
# mcar_test(mc_resp)
# mcar_test(airquality)
# mcar_test(riskfactors)
# 
# chisq.test(data)
# 
# r1=factor(rr)
# t.test(Response~r1)
# str(r1)
