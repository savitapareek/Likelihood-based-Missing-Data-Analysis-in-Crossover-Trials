## cleaning the enviornment
rm(list=ls())
gc()

## install required packages
library(psych) ## for trace of a matrix
library(doParallel) ## for parallelization
library(MASS) ## for mvtnorm

### function definition for parameter estimation with missing
updateestimates=function(beth,sigmah,sigmash,booty,bootx,bootz,ymisgivyobs){
  iD=Eb=list()
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

##########---- cretaing matrices and data generation
strt=Sys.time()
s=2#no of sequences
p=t=2#no of periods and no of treatments
m=4#no of  genes

## true values for simulation
tsigma=1.2^2	;tsigmas=.7^2
tbet=array()
cat("enter true value of fixed effects total fixed effects are",p+t+m-2)
tbet=c(4.5,	.2,	1.06,	.46,	1.09,	.50	)
truest=c(tbet,tsigma,tsigmas)
no=array()
cat("enter missing probabilities in each sequence, and total seq=",s)
prob=c(.15,.15)
#prob=scan(nmax=s) 
n=pmn=n1=list()## no of subjects in different sequences
no=c(50,50)
for (i in 1:s) {
  cat("no of subjects in seq",i,"\n")
  # no[i]=as.integer(readline(prompt="Enter an integer: "))
  n[[i]]=diag(no[i])
  pmn[[i]]=diag(no[i]*p*m)
  n1[[i]]=rep(1,no[i])
}
cat("enter no of generated missing observations in sequence 1,2,..,",s)
#num= scan(nmax=s)
num=c(2000,2000)
treatment=list()
treatment[[1]]=c(1,2)#,3,4)
treatment[[2]]=c(2,1)#,3,4)

### making X matrix
Trt=list()
for (i in 1:s) {
  k=0; k3=m
  Trt[[i]]=matrix(nrow=p*m,ncol = (t-1))
  arr=array(dim = t-1)
  for (j in 1:p) {
    if(treatment[[i]][j]==t)
      arr[1:(t-1)]=0
    for (tt in 1:(t-1)) {
      if(treatment[[i]][j]==tt)
      {arr[tt]=1
      arr[-tt]=0
      }
    }
    Trt[[i]][((k+1):k3),]=kronecker(rep(1,m),rbind(arr))
    k=j*m
    k3=(j+1)*m
  }
}
x1k=list()
for (i in 1:s) {
  x1k[[i]]=cbind(kronecker(rep(1,p),rep(1,m)),rbind(kronecker(diag(p-1),rep(1,m)),matrix(data=0,nrow = m,ncol = p-1)),Trt[[i]],kronecker(rep(1,p),rbind(diag(m-1),rep(0,m-1))))
}
x=z=list()
x=mcMap(kronecker,n1,x1k)

#making z matrix
z1k=cbind(kronecker(rep(1,p),rep(1,m)))
z=mclapply(n,kronecker,z1k)

### In oreder to perform 500 simulations we read data 
nrep=500
myy=list()
myy[[1]]=read.csv(file.choose(),header = T) ## read 'sample in seq1.csv'
myy[[2]]=read.csv(file.choose(),header = T) ## read 'sample in seq2.csv'

##****** incorporating missing*******
## estimates are saved in festi matrix
festi=matrix(nrow = nrep,ncol =2*(p+t+m) )

## for loop to run 500 simulations one by one
for(sim in 1:nrep){
  ## initial estimates for each simulation
  beth=tbet-.2 
  sigmah=tsigma-.1
  sigmash=tsigmas+.2
  missubindctr=misindctr=misx=obsx=misz=obsz=countmis=countobs=obsy=misy=w=missub= per=newy=D=S=list()
  for (i in 1:s) {
    newy[[i]]=array(dim  = p*m*no[i])#
    newy[[i]]=as.numeric(myy[[i]][sim,-1])
  }
  # newy=lapply(newy,t) 
  repeat{
    for (i in 1:s) {
      missubindctr[[i]]=rbinom((p*no[[i]]),1,prob[i]) 
      w[[i]]=which(missubindctr[[i]]==1)
    }
    if(length(which(lapply(w, length)==0))==0)
      break
  }
  for (i in 1:s) {
    misindctr[[i]]=rep(0,p*m*no[i])
    for (o in 1:length(w[[i]])) {
      misindctr[[i]][(m*(w[[i]][o]-1)+1):(m*(w[[i]][o]-1)+m)]=1
    }
    countmis[[i]]=diag(table(misindctr[[i]])[2])
    countobs[[i]]=diag(table(misindctr[[i]])[1])
    misx[[i]]=subset(x[[i]],misindctr[[i]]==1)
    misz[[i]]=subset(z[[i]],misindctr[[i]]==1)
    obsx[[i]]=subset(x[[i]],misindctr[[i]]==0)
    obsz[[i]]=subset(z[[i]],misindctr[[i]]==0)
    obsy[[i]]=subset(newy[[i]],misindctr[[i]]==0)
    misy[[i]]=subset(newy[[i]],misindctr[[i]]==1)
  }
  missub=mclapply(mclapply(w,"/",p), ceiling)
  per=mclapply(w, "%%",p)
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
  #number of missing oservations generated
  for (i in 1:s) {
    ymisgivyobs[[i]]=mvrnorm(n=num[i],meanymisgivobsy[[i]],varymisgivobsy[[i]])
  }
  estti=array()
  estti=updateestimates(beth,sigmah,sigmash,newy,x,z,ymisgivyobs)
  ####################final algorithmmmmmmmm
  repeat {
    beth=as.numeric(estti[1:(p+t+m-2)])
    sigmah=as.numeric(estti[p+t+m-1])
    sigmash=as.numeric(estti[p+t+m])
    estti=array()
    estti=  updateestimates(beth,sigmah,sigmash,newy,x,z,ymisgivyobs)
    if(max(abs(beth-estti[1:(p+t+m-2)]))<5*10**-4&abs(sigmah-estti[p+t+m-1] ) <5*10**-4&abs(sigmash-estti[p+t+m] ) <5*10**-4 )
      break
  }
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
    iD=iS=list()
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
    yminxbet=lapply(1:s, function(i) matrix(ymin1,ncol = s)[,i])
    nnewy= lapply(1:s, function(i) matrix(nny1,ncol = s)[,i])
    
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
  festierr=c(estti,fstderr)
  print(festierr)
  festi[sim,] =festierr
}  

###******calculation for coverage probability****
covprob=array(dim = length(truest))
covprob[1:(p+t+m)]=0
bias=meansq=matrix(nrow = nrep,ncol = (p+t+m))
for (i in 1:nrep) {
  mytable=cbind(truest,festi[i,(1:(p+t+m))],festi[i,-(1:(p+t+m))])
  bias[i,]= (festi[i,(1:(p+t+m))]-truest)
  meansq[i,]=((festi[i,-(1:(p+t+m))])**2)+ (festi[i,(1:(p+t+m))]-truest)**2
  lowerci=mytable[,2]-1.96*mytable[,3]
  upperci=mytable[,2]+1.96*mytable[,3]
  for (j in 1:(p+t+m)) {
    if(lowerci[j]<truest[j]&truest[j]<upperci[j])
      covprob[j]=covprob[j]+1
  }
}
estimat1=cbind(colMeans(festi[1:nrep,(1:(p+t+m))]),colMeans(festi[1:nrep,-(1:(p+t+m))]) ,abs(colMeans(bias)),colMeans(meansq),covprob/nrep)
colnames(estimat1)=c("estimate","stderror", "bias","meansq","cprob")
ae=round(cbind(truest,estimat1),3)
ae
Sys.time()-strt
#write.csv(festi,"D:\ for 15% missing all sim estimates and se.csv")
#writ.csv(ae,"D:\ 15% and comp est.csv")

