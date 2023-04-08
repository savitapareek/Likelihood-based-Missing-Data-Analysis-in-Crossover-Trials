## cleaning the enviornment
rm(list=ls())
gc()
strt=Sys.time()

## required libraies
library(psych) 
library(doParallel)
library(MASS)
library(nlme) ## for lme function
library(lmeInfo)## for varcomp_vcov function
library(insight) ##get_variance_random function
library(PDFEstimator)## estimate pdf
library(moments)## for skewness
library(ggplot2)## for plots geom_density
library(ggpubr)## for ggarrange
library(mice)## for mice
###*******function for evaluating missing data estimates*****
updateestimates=function(beth,sigmah,sigmash,booty,bootx,bootz,ymisgivyobs){
  iD=Eb=list()#initial estimate of D
  iD=mclapply(n,"*",sigmash)
  iS=mcMap("+",mcMap("%*%",mcMap("%*%",bootz,iD),lapply(bootz,t)),mclapply(pmn,'*',sigmah))
  b0i=sigma0i=fulyminxbet=xbet=zboi=fulyminxbtminzb=list()
  sigma0i=mclapply(mcMap("+",mclapply(mclapply(bootz, crossprod), "*",sigmah**-1),mclapply(iD, solve)),solve)
  xbet=mclapply(bootx, "%*%",beth)
  sigma0iz=mcMap("%*%",sigma0i,lapply(bootz, t))
  ### missing data is replaced with generated data num is missing obs and no reprsents no of objects in each sequence
  sumresd=bb=ymisobs=list()
  for (i  in 1:s) {
    ymisobs[[i]]= fulyminxbet[[i]]=matrix(nrow=num[[i]],ncol = (p*m*(no[i])))
    b0i[[i]]=matrix(nrow=num[[i]],ncol = ((1)*(no[i])))
    bb[[i]]=array()
    k=1
    for (j in 1:num[i]) {
      ymisobs[[i]][k,]=  replace(booty[[i]], which(misindctr[[i]]==1),  ymisgivyobs[[i]][j,])
      fulyminxbet[[i]][k,]=ymisobs[[i]][k,]-xbet[[i]]
      b0i[[i]][k,]=sigma0iz[[i]]%*%fulyminxbet[[i]][k,]/(sigmah)
      bb[[i]][k]=crossprod(b0i[[i]][k,])
      k=k+1
    }
  }
  zb0i= mcMap("%*%",bootz,lapply(b0i,t))
  yminzb=mcMap("-",ymisobs,lapply(zb0i, t))
  sumyminzb=mcMap("/",mclapply(yminzb, colSums),num)
  xyminzb=mcMap("%*%",lapply(bootx, t),sumyminzb)
  bethat=solve(Reduce('+',mclapply(bootx, crossprod)))%*%Reduce('+',xyminzb)
  xbett=list()
  xbett=mclapply(bootx, "%*%",bethat)
  for (i  in 1:s) {
    fulyminxbtminzb[[i]]=matrix(nrow=num[[i]],ncol = (p*m*(no[i])))
    sumresd[[i]]=array()
    k=1
    for (j in 1:num[i]) {
      fulyminxbtminzb[[i]][k,]=(ymisobs[[i]][k,])-xbett[[i]]-(bootz[[i]]%*%b0i[[i]][k,])
      sumresd[[i]][k]=crossprod( fulyminxbtminzb[[i]][k,])
      k=k+1
    }
  }
  trzs=Reduce('+',mclapply(mcMap("%*%",mclapply(bootz, crossprod),sigma0i),tr))
  sumsigmas1=Reduce("+", mclapply(sumresd, mean))
  fstsigma=(sumsigmas1+trzs)/(p*m*sum(no))
  trs=Reduce('+',mclapply(sigma0i,tr))
  sigssum2=Reduce("+", mclapply(bb, mean))
  fstsigmas=(trs+sigssum2)/sum(no)
  return(c(bethat,fstsigma,fstsigmas))
}

### complete case analysis each sequence has 4 subjects removing missing subject completely
s=3#no of sequences
p=t=3#no of periods and no of treatments
m=10#no of  genes
n=pmn=n1=list()## no of objects in different sequences
no=c(4,4,4)
for (i in 1:s) {
  n[[i]]=diag(no[i])
  pmn[[i]]=diag(no[i]*p*m)
  n1[[i]]=rep(1,no[i])
}
cat("enter no of generated missing observations in sequence 1,2,..,",s)
treatment=list()
## we assume 1: placebo; 2: 10 mg; 3: 25 mg
treatment[[1]]=c(2,1,3)
treatment[[2]]=c(3,2,1)
treatment[[3]]=c(1,3,2)
Trt=list()
for (i in 1:s) {
  k=0; k3=m
  Trt[[i]]=matrix(nrow=p*m,ncol = (t-1))
  arr=array(dim = t-1)
  for (j in 1:p) {
    if(treatment[[i]][j]==1)
      arr[1:(t-1)]=0
    for (tt in 2:t) {
      if(treatment[[i]][j]==tt)
      {arr[tt-1]=1
      arr[-(tt-1)]=0
      }
    }
    Trt[[i]][((k+1):k3),]=kronecker(rep(1,m),rbind(arr))
    k=j*m
    k3=(j+1)*m
  }
}

x1k=x1r=list()

lu= 5 ## which estimates p value you want: say for \beta_1
for (i in 1:s) {
  x1k[[i]]=cbind(rep(1,p*m),rbind(matrix(data=0,nrow = m,ncol = p-1),kronecker(diag(p-1),rep(1,m))),Trt[[i]],kronecker(rep(1,p),rbind(rep(0,m-1),diag(m-1))))#,carry[[i]]
  x1r[[i]]=x1k[[i]][,-lu]
}
x=mcMap(kronecker,n1,x1k)

#making z matrix
z1k=cbind(kronecker(rep(1,p),rep(1,m)))
z=mclapply(n,kronecker,z1k)

#read anova 3 way allseq.csv file: complete case analysis
data11=read.csv(file.choose(),header = T)
attach(data11)
kk=0
sub=array()
for (i in 1:(nrow(data11)/30)) {
  sub[((kk+1):(kk+30))]=rep(i,30)
  kk=kk+30
}
seq=c(rep(1,120),rep(2,120),rep(3,120))
data11$subject=as.factor(sub)
seq=as.factor(seq)
data11$seq=seq

## fit model with original or log responses
treatment1=c(rep(rep(c("trt2","trt1","trt3"),each=m),no[1]),
             rep(rep(c("trt3","trt2","trt1"),each=m),no[2]),
             rep(rep(c("trt1","trt3","trt2"),each=m),no[3]))
subject=rep(paste("sub",1:sum(no)),each=(p*m))#+subject*gene


model1=lme((Response_allseq)~ factor(Period)+factor(Treatment)+factor(Gene),data=data11,random = ~ 1|subject,method = 'ML')# intercept vary randomly with the subjects
#summary(model1)
#for log responses
#model1=lme(log(Response_allseq)~ Period+ Treatment+Gene,data=data11,random = ~ 1|subject,method = 'ML')# intercept vary randomly with the subjects
isigma=round(summary(model1)$sigma^2,5)
isigmas=as.numeric(VarCorr(model1)[1,1])
ii=round(as.numeric((summary(model1)$coefficients)$fixed),3)#round(data.table(summary(model1)$coefficients)[,1]$Estimate,3)
ibet=ii[c(1:5,7:14,6)]
iD=mclapply(n,"*",isigmas)
sig0i=mclapply(mcMap("+",mclapply(mclapply(z, crossprod), "*",isigma**-1),mclapply(iD, solve)),solve)
ymisobs=list(Response_allseq[1:120],
             Response_allseq[121:240],
             Response_allseq[241:360])
b0i= lapply(1:s, function(i) (
  sig0i[[i]]%*%t(z[[i]])%*%(ymisobs[[i]]-x[[i]]%*%ibet)/isigma
))
L1=Reduce("+", lapply(1:s, function(i) (
  -.5*(p*m*no[i]*log(isigma)+(1/isigma)*(tr(crossprod(z[[i]]%*%sig0i[[i]]))
                                         +crossprod(ymisobs[[i]]-x[[i]]%*%ibet
                                                     -z[[i]]%*%b0i[[i]])) 
       +no[i]*log(isigmas)+tr(solve(iD[[i]])%*%sig0i[[i]])
       + t(b0i[[i]])%*%solve(iD[[i]])%*%b0i[[i]] ))))
## compute aic of model
aic=-2*(L1) + 2*(p+t+m)
bic=-2 * L1 + log(sum(no)*p*m) *(p+t+m)
c(aic,bic)
## table 1: complete cases; estimate, se,p-value rmse is taken from summary(model)
summary(model1)
sqrt(diag(varcomp_vcov(model1)))## for standard error of variance components
estii=cbind(c(ibet,isigma,isigmas))
sqrt(diag(varcomp_vcov(model1)))[c(2,1)]
stderr=as.numeric(c(sqrt(diag(vcov(model1)))[c(1:5,7:14,6)],sqrt(diag(varcomp_vcov(model1)))[c(2,1)]))
pvalu=c(summary(model1)$tTable[,5][c(1:5,7:14,6)],0,0)
comp_case1=cbind(estii,stderr,pvalu)
rmse_comp1=summary(model1)$sigma

## ###Imputation analysis AIC%%%%%********############
dat=read.csv(file.choose())## read all responses in gene data 3,3
respos=dat[,1]
no=c(6,6,5)

for (i in 1:s) {
  n[[i]]=diag(no[i])
  pmn[[i]]=diag(no[i]*p*m)
  n1[[i]]=rep(1,no[i])
}
## fit model with original or log responses
treatment1=c(rep(rep(c("trt2","trt1","trt3"),each=m),no[1]),
             rep(rep(c("trt3","trt2","trt1"),each=m),no[2]),
             rep(rep(c("trt1","trt3","trt2"),each=m),no[3]))
subject=rep(paste("sub",1:sum(no)),each=(p*m))#+subject*gene
Gene=rep(paste("gene",1:10),p*sum(no))
Period=rep(rep(paste("per",1:3),each=m),sum(no))
dat11=data.frame(respos,Period,Gene,treatment1,subject)
rom=which(is.na(respos)==T)
d2=na.omit(dat11)
imputed_Data <- mice(dat11, m=5, maxit = 50, method = 'pmm', seed = 500)
summary(imputed_Data)
completeData <- rowMeans(sapply(1:5, function(i) complete(imputed_Data,i)$respos))
#length(completeData[,1])
#nnewy2=list(completeData[1:(p*m*no[1])],completeData[((p*m*no[1])+1):(p*m*(no[1]+no[2]))],completeData[(p*m*(no[1]+no[2])+1):(p*m*sum(no))])

nrow(dat11)
nresp=replace(respos, which(is.na(respos)==T), completeData[rom])
dat11=data.frame(nresp,Period,Gene,treatment1,subject)

model1=lme((nresp)~ factor(Period)+factor(treatment1)+factor(Gene),data=dat11,random = ~ 1|subject,method = 'ML')# intercept vary randomly with the subjects
extractAIC(model1)
## AIC values
isigma=round(summary(model1)$sigma^2,5)
isigmas=as.numeric(VarCorr(model1)[1,1])
ii=round(as.numeric((summary(model1)$coefficients)$fixed),3)#round(data.table(summary(model1)$coefficients)[,1]$Estimate,3)
ibet=ii[c(1:5,7:14,6)]

## estimates std error and p-value from mice are
estii=cbind(c(ibet,isigma,isigmas))
sqrt(diag(varcomp_vcov(model1)))[c(2,1)]
stderr=as.numeric(c(sqrt(diag(vcov(model1)))[c(1:5,7:14,6)],sqrt(diag(varcomp_vcov(model1)))[c(2,1)]))
pvalu=c(summary(model1)$tTable[,5][c(1:5,7:14,6)],0,0)
mic=cbind(estii,stderr,pvalu)
rmse_mic=summary(model1)$sigma
compmic=cbind(comp_case1,mic)

## for AIC BIC values
iD=mclapply(n,"*",isigmas)
ymisobs=list(na.omit(nresp[1:180]),
             na.omit(nresp[181:360]),
             na.omit(nresp[361:510]))
ress=list(nresp[1:180],nresp[181:360],nresp[361:510])
xx=mcMap(kronecker,n1,x1k)
nx=lapply(1:s, function(i) xx[[i]][-which(is.na(ress[[i]])==T),])
#nnewy=lapply(1:s, function(i) newy[[i]][-make11[[i]]])
z1kk=cbind(kronecker(rep(1,p),rep(1,m)))
zz=mclapply(n,kronecker,z1kk)
nz=lapply(1:s, function(i) zz[[i]][-which(is.na(ress[[i]])==T),])
sig0i=mclapply(mcMap("+",mclapply(mclapply(nz, crossprod), "*",isigma**-1),mclapply(iD, solve)),solve)

b0i= lapply(1:s, function(i) (
  sig0i[[i]]%*%t(nz[[i]])%*%(ymisobs[[i]]-nx[[i]]%*%ibet)/isigma
))
L1=Reduce("+", lapply(1:s, function(i) (
  -.5*(p*m*no[i]*log(isigma)+(1/isigma)*(tr(crossprod(nz[[i]]%*%sig0i[[i]]))
                                         +crossprod(ymisobs[[i]]-nx[[i]]%*%ibet
                                                    -nz[[i]]%*%b0i[[i]])) 
       +no[i]*log(isigmas)+tr(solve(iD[[i]])%*%sig0i[[i]])
       + t(b0i[[i]])%*%solve(iD[[i]])%*%b0i[[i]] ))))
## compute aic of model
aic=-2*(L1) + 2*(p+t+m)
bic=-2 * L1 + log(sum(no)*p*m) *(p+t+m)
c(aic,bic)


##MAR case when we consider the whole data number of subjects are 6, 6, 5
data1=read.csv(file.choose(),header = T)# read seq1_9 for 9%, seq1_21 for 21%, seq1_24 for 24% missing
data2=read.csv(file.choose(),header = T)#read seq2_9 for 9%, seq2_21 for 21%, seq2_24 for 24% missing
data3=read.csv(file.choose(),header = T)# read seq3_9 for 9%, seq3_21 for 21%, seq3_24 for 24% missing

n=pmn=n1=y=list()## no of subjects in different sequences
no=c(6,6,5)
for (i in 1:s) {
  n[[i]]=diag(no[i])
  pmn[[i]]=diag(no[i]*p*m)
  n1[[i]]=rep(1,no[i])
}
cat("enter no of generated missing observations in sequence 1,2,..,",s)
num=c(2000,2000,2000)

x=mcMap(kronecker,n1,x1k)
xr=mcMap(kronecker,n1,x1r)
z=mclapply(n,kronecker,z1k)

### write log if you want log responses
y[[1]]=(data1[,2])# log(data1[,2])
y[[2]]=(data2[,2])
y[[3]]=(data3[,2])
newy=y
prob=array(dim=s)
for (i  in 1:s) {
  prob[i]=  length(which(is.na(y[[i]])==T))/length(y[[i]])
}
print(prob)## missing probailities
mean(prob)#average missing response

## now we take initial estimates as complete case final estimates
isigma=as.numeric(comp_case1[,1][15])
isigmas=as.numeric(comp_case1[,1][16])
ibet=as.numeric(comp_case1[,1][1:14])
#####****incorporating missing****
## estimate under H1
beth=ibet
sigmah=isigma
sigmash=isigmas
missubindctr=misindctr=misx=obsx=misz=obsz=countmis=countobs=obsy=misy=w=missub= per=D=S=list()
for (i in 1:s) {
  w[[i]]=which(is.na(y[[i]])==T)
  misindctr[[i]]=rep(0,p*m*no[i])
  misindctr[[i]]= replace(misindctr[[i]],w[[i]],1)
  countmis[[i]]=diag(length(w[[i]]))
  countobs[[i]]=diag((p*m*no[i]-length(w[[i]])))
  misx[[i]]=subset(x[[i]],misindctr[[i]]==1)
  misz[[i]]=subset(z[[i]],misindctr[[i]]==1)
  obsx[[i]]=subset(x[[i]],misindctr[[i]]==0)
  obsz[[i]]=subset(z[[i]],misindctr[[i]]==0)
  obsy[[i]]=subset(y[[i]],misindctr[[i]]==0)
  misy[[i]]=subset(y[[i]],misindctr[[i]]==1)
}
missub=mclapply(mclapply(w,"/",p*m), ceiling)#which subjects are missing
iD=iS=list()
iD=mclapply(n,"*",sigmash)
iS=mcMap("+",mcMap("%*%",mcMap("%*%",z,iD),lapply(z,t)),mclapply(pmn,'*',sigmah))
mumis=mclapply(misx, "%*%",beth)
muobs=mclapply(obsx, "%*%",beth)
S11=mcMap("+",mcMap("%*%",mcMap("%*%",misz,iD),lapply(misz, t)),mclapply(countmis, "*",sigmah))
S12=mcMap("%*%",mcMap("%*%",misz,iD),lapply(obsz, t))
S21=mcMap("%*%",mcMap("%*%",obsz,iD),lapply(misz, t))
S22=mcMap("+",mcMap("%*%",mcMap("%*%",obsz,iD),lapply(obsz, t)),mclapply(countobs, "*",sigmah))
varymisgivobsy=meanymisgivobsy=list()
meanymisgivobsy=mcMap("+",mumis,mcMap("%*%",mcMap("%*%", S12,mclapply(S22, solve)),mcMap("-",mclapply(obsy,matrix),muobs)))
varymisgivobsy=mcMap("-",S11,mcMap("%*%",mcMap("%*%", S12,mclapply(S22, solve)),S21))
ymisgivyobs=list()
#number of missing observations generated
for (i in 1:s) {
  ymisgivyobs[[i]]=mvrnorm(n=num[i],meanymisgivobsy[[i]],varymisgivobsy[[i]])
}
estti=array()
estti=updateestimates(beth,sigmah,sigmash,newy,x,z,ymisgivyobs)
zz=0
####################final algorithmmmmmmmm
repeat {
  zz=zz+1
  beth=as.numeric(estti[1:(p+t+m-2)])
  sigmah=as.numeric(estti[p+t+m-1])
  sigmash=as.numeric(estti[p+t+m])
  estti=array()
  estti=  updateestimates(beth,sigmah,sigmash,newy,x,z,ymisgivyobs)
  if(max(abs(beth-estti[1:(p+t+m-2)]))<5*10**-4&abs(sigmah-estti[p+t+m-1] ) <5*10**-4&abs(sigmash-estti[p+t+m] ) <5*10**-4 )
    break
}
zz# no of iterations to converge
## using these estimates generate one value of ymis,i and bi and use complete data estimation(for generating one obs gibbs sampling is used)
#### to find parameter estimates and there variance
#### repeat the process b=100 times
ntimes=100
bet=varbet=array()
sigmaesq=sigmassq=varsigmaesq=varsigmassq=array(dim =ntimes )
estimatess=matrix(ncol = 2*(p+t+m),nrow = (ntimes))
cores = detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)
estimatess=foreach(imp = 1:ntimes, .combine = rbind,.packages = c('parallel','MASS','psych') )%dopar% {
  beth=estti[1:(p+t+m-2)];sigmah=estti[p+t+m-1];sigmash=estti[p+t+m]
  iD=mclapply(n,"*",sigmash)
  iS=mcMap("+",mcMap("%*%",mcMap("%*%",z,iD),lapply(z,t)),mclapply(pmn,'*',sigmah))
  mumis=mclapply(misx, "%*%",beth)
  S11=mcMap("+",mcMap("%*%",mcMap("%*%",misz,iD),lapply(misz, t)),mclapply(countmis, "*",sigmah))
  muobs=mclapply(obsx, "%*%",beth)
  varymisgivobsy=mclapply(countmis, "*",sigmah)
  ymisgivyobs=  ymisobs=bgivy=nnewy=meanbgivy=varbgivy=newyy=yymisy=yminxbet=list()
  varbgivy=mclapply(mcMap("+",mclapply(mclapply(z, crossprod), "*",sigmah**-1),mclapply(iD, solve)),solve)
  for (i in 1:s) {
    yymisy[[i]]=mvrnorm(n=1,mumis[[i]],S11[[i]])
    newyy[[i]]= replace(newy[[i]], which(misindctr[[i]]==1),  yymisy[[i]])
  }
  
  #number of missing oservations generated 
  # gibbs sampling from bi,ymis,i joint distribution
  nny=yminxbt=list()
  for (burnin in 2:num[1]) {
    for (i in 1:s) {
      meanbgivy[[i]]=(1/sigmah)*varbgivy[[i]]%*%t(z[[i]])%*%(newyy[[i]]-x[[i]]%*%beth)
      bgivy[[i]]=mvrnorm(n=1,meanbgivy[[i]],varbgivy[[i]])
      meanymisgivobsy[[i]]=mumis[[i]]+misz[[i]]%*%(as.matrix(bgivy[[i]]))
      ymisgivyobs[[i]]=mvrnorm(n=1,meanymisgivobsy[[i]],varymisgivobsy[[i]])
      nnewy[[i]]=  replace(newy[[i]], which(misindctr[[i]]==1),  ymisgivyobs[[i]])
      yminxbet[[i]]= nnewy[[i]]-x[[i]]%*%beth
    }
    nny[[burnin]]=nnewy
    yminxbt[[burnin]]=yminxbet
  }
  ### get our final one sample as average of burnin1+1 to num[1] samples
  #for better estimates and sample
  burnin1=.20*num[1]
  nny1=Reduce("+",lapply(nny, unlist)[(burnin1+1):num[1]])/(num[1]-burnin1)
  ymin1=Reduce("+",lapply(yminxbt, unlist)[(burnin1+1):num[1]])/(num[1]-burnin1)
  yminxbet=list()
  yminxbet[[1]]=ymin1[1:180];yminxbet[[2]]=ymin1[181:360];yminxbet[[3]]=ymin1[-c(1:360)]
  nnewy=list()
  nnewy[[1]]=nny1[1:180];nnewy[[2]]=nny1[181:360];nnewy[[3]]=nny1[-c(1:360)]  #= lapply(1:s, function(i) nny1[((i-1)*length(newy[[i]])+1):(i*length(newy[[i]]))])
  
  xtrxinv= solve( Reduce("+",mclapply(x,crossprod)))
  yminzb= Map("-",nnewy,Map("%*%",z,bgivy))
  xtryminzb=Reduce("+",Map("%*%",lapply(x, t),yminzb))
  bet=xtrxinv%*%xtryminzb
  yminxbetminzb=Map("-",yminxbet,Map("%*%",z,bgivy))
  sigmaesq= Reduce("+",mclapply(yminxbetminzb,crossprod))/(p*m*sum(no))
  sigmassq=Reduce("+",mclapply(bgivy,crossprod))/sum(no)
  varbet=as.numeric(sigmaesq)*diag(xtrxinv)
  varsigmaesq=2*sigmaesq**2/(p*m*sum(no))
  varsigmassq=2*sigmassq**2/(sum(no))
  c(bet,sigmaesq,sigmassq,varbet,varsigmaesq,varsigmassq)
}
stopCluster(cl)
est=cbind(estimatess[,(1:(p+t+m))])
varest=cbind(estimatess[,-(1:(p+t+m))])
fvarbeth=colMeans(varest)+(1+(1/ntimes))%*%apply(est, 2, var)
fstderr=sqrt(as.vector(abs(fvarbeth)))
festierr1=cbind(estti,fstderr)

## now save it in a table 
write.csv(cbind(comp_case1,mic,festierr1),"D:\ real1.csv")

###  likelihood under H_1 for finding the aic bic of proposed analysis
esttif=estti
beth=esttif[1:(p+t+m-2)];sigmae=esttif[p+t+m-1];sigmash=esttif[p+t+m]
iD=mclapply(n,"*",sigmash)
sig0i=mclapply(mcMap("+",mclapply(mclapply(z, crossprod), "*",sigmae**-1),mclapply(iD, solve)),solve)
mumis=mclapply(misx, "%*%",beth)
muobs=mclapply(obsx, "%*%",beth)
S11=mcMap("+",mcMap("%*%",mcMap("%*%",misz,iD),lapply(misz, t)),mclapply(countmis, "*",sigmae))
S12=mcMap("%*%",mcMap("%*%",misz,iD),lapply(obsz, t))
S21=mcMap("%*%",mcMap("%*%",obsz,iD),lapply(misz, t))
S22=mcMap("+",mcMap("%*%",mcMap("%*%",obsz,iD),lapply(obsz, t)),mclapply(countobs, "*",sigmae))
varymisgivobsy=meanymisgivobsy=ymisobs=list()
meanymisgivobsy=mcMap("+",mumis,mcMap("%*%",mcMap("%*%", S12,mclapply(S22, solve)),mcMap("-",mclapply(obsy,matrix),muobs)))
varymisgivobsy=mcMap("-",S11,mcMap("%*%",mcMap("%*%", S12,mclapply(S22, solve)),S21))

set.seed(11) ## to reproduce the results
ymisgivyobs=lapply(1:s, function(i) mvrnorm(n=num[i],meanymisgivobsy[[i]],varymisgivobsy[[i]]))
ymisobs=lapply(1:s, function(i) t(sapply(1:num[i],function(j) (
  replace(newy[[i]], which(misindctr[[i]]==1),  ymisgivyobs[[i]][j,])
)) ))
b0i= lapply(1:s, function(i) t(sapply(1:num[i],function(j) (
  sig0i[[i]]%*%t(z[[i]])%*%(ymisobs[[i]][j,]-x[[i]]%*%beth)/sigmae
)) ))
## for residual analysis
re= lapply(1:s, function(i) rowMeans(sapply(1:num[i],
                                            function(j) ymisobs[[i]][j,]-x[[i]]%*%beth-z[[i]]%*%b0i[[i]][j,])))
residu= c(re[[1]],re[[2]],re[[3]])
rmse=sqrt(var(residu))
rmse #root mean square value

## for L1 
L1=Reduce("+", lapply(1:s, function(i) (
  -.5*(p*m*no[i]*log(sigmae)+(1/sigmae)*(tr(crossprod(z[[i]]%*%sig0i[[i]]))
                                         +mean(sapply(1:num[i], function(j) 
                                           crossprod(ymisobs[[i]][j,]-x[[i]]%*%beth
                                                     -z[[i]]%*%b0i[[i]][j,])) ))
       +no[i]*log(sigmash)+tr(solve(iD[[i]])%*%sig0i[[i]])
       +mean(sapply(1:num[i], function(k) t(b0i[[i]][k,])%*%solve(iD[[i]])%*%b0i[[i]][k,] )))
))) 
## compute aic of model
aic=-2*(L1) + 2*(p+t+m)  
bic=-2 * L1 + log(sum(no)*p*m) *(p+t+m)
c(aic,bic)
Sys.time()-strt # time to get estimates and std errors
