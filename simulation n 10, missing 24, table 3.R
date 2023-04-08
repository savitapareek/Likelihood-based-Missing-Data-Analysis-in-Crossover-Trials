## cleaning the enviornment
rm(list=ls())
gc()

## install required packages
library(psych) ## for trace of a matrix
library(doParallel) ## for parallelization
library(MASS) ## for mvtnorm
library(mice)
### function definition for parameter estimation with missing
updateestimates=function(beth,sigmah,sigmash,x,y,z){
  iD=mclapply(n,"*",sigmash)
  iS=mcMap("+",mcMap("%*%",mcMap("%*%",z,iD),lapply(z,t)),mclapply(pmn,'*',sigmah))
  b0i=fulyminxbet=fulyminxbtminzb= sumresd=bb=ymisobs=list()
  sigma0i=mclapply(mcMap("+",mclapply(mclapply(z, crossprod), "*",sigmah**-1),mclapply(iD, solve)),solve)
  xbet=mclapply(x, "%*%",beth)
  sigma0iz=mcMap("%*%",sigma0i,lapply(z, t))
  ### missing data is replaced with generated data num is missing obs and no reprsents no of objects in each sequence
  
  for (i  in 1:s) {
    ymisobs[[i]]= fulyminxbet[[i]]=matrix(nrow=num[[i]],ncol = (p*m*(no[i])))
    b0i[[i]]=matrix(nrow=num[[i]],ncol = ((1)*(no[i])))
    bb[[i]]=array()
    k=1
    for (j in 1:num[i]) {
      ymisobs[[i]][k,]=  replace(y[[i]], which(misindctr[[i]]==1),  ymisgivyobs[[i]][j,])
      fulyminxbet[[i]][k,]=ymisobs[[i]][k,]-xbet[[i]]
      b0i[[i]][k,]=sigma0iz[[i]]%*%fulyminxbet[[i]][k,]/(sigmah)
      bb[[i]][k]=crossprod(b0i[[i]][k,])
      k=k+1
    }
  }
  zb0i= mcMap("%*%",z,lapply(b0i,t))
  yminzb=mcMap("-",ymisobs,lapply(zb0i, t))
  sumyminzb=mcMap("/",mclapply(yminzb, colSums),num)
  xyminzb=mcMap("%*%",lapply(x, t),sumyminzb)
  bethat=solve(Reduce('+',mclapply(x, crossprod)))%*%Reduce('+',xyminzb)
  xbett=list()
  xbett=mclapply(x, "%*%",bethat)
  for (i  in 1:s) {
    fulyminxbtminzb[[i]]=matrix(nrow=num[[i]],ncol = (p*m*(no[i])))
    sumresd[[i]]=array()
    k=1
    for (j in 1:num[i]) {
      fulyminxbtminzb[[i]][k,]=(ymisobs[[i]][k,])-xbett[[i]]-(z[[i]]%*%b0i[[i]][k,])
      sumresd[[i]][k]=crossprod( fulyminxbtminzb[[i]][k,])
      k=k+1
    }
  }
  trzs=Reduce('+',mclapply(mcMap("%*%",mclapply(z, crossprod),sigma0i),tr))
  sumsigmas1=Reduce("+", mclapply(sumresd, mean))
  fstsigma=(sumsigmas1+trzs)/(p*m*sum(no))
  trs=Reduce('+',mclapply(sigma0i,tr))
  sigssum2=Reduce("+", mclapply(bb, mean))
  fstsigmas=(trs+sigssum2)/sum(no)
  return(c(bethat,fstsigma,fstsigmas))
}

### function definition for parameter estimation without missing
updateestimate1=function(ibet,isigma,isigmas,x,y,z){
  iD=mclapply(nn,"*",isigmas)
  yminxbt=mcMap("-", y,mclapply(x, "%*%" , ibet) )
  iS=mcMap("+",mcMap("%*%",mcMap("%*%",z,iD),lapply(z,t)),mclapply(pmn1,'*',isigma))
  xtrxinv= solve( Reduce("+",mclapply(x,crossprod)))
  sigma0i=b0i=yminzb=list()
  sigma0i=mclapply(mcMap("+",mclapply(mclapply(z, crossprod), "*",isigma**-1),mclapply(iD, solve)),solve)
  b0i=mclapply(mcMap("%*%",mcMap("%*%",sigma0i,lapply(z,t)),yminxbt),"/",isigma)
  xyminzb=mcMap("%*%",lapply(x,t),mcMap("-",y,mcMap("%*%",z,b0i)))
  bethat=xtrxinv%*%Reduce("+",xyminzb)
  yminxbt=mcMap("-", y,mclapply(x, "%*%" , bethat))
  t1=  Reduce("+",mclapply(yminxbt,crossprod))-  Reduce("+",mclapply(mcMap("%*%",lapply(yminxbt,t),mcMap("%*%",z,b0i)),"*",2))
  t2=Reduce("+",mcMap("%*%",mcMap("%*%",lapply(b0i,t),mclapply(z,crossprod)),b0i))+Reduce("+",mclapply(mcMap("%*%",mclapply(z,crossprod),sigma0i),tr))
  fstsigma=(t1+t2) /(p*m*sum(no1))## check the denominator
  fstsigmas=(Reduce("+",mclapply(b0i,crossprod))+Reduce("+",mclapply(sigma0i,tr)))/sum(no1)
  return(c(bethat,fstsigma,fstsigmas))
}

##########---- cretaing matrices and data generation
strt=Sys.time()
s=3#no of sequences
p=t=3#no of periods and no of treatments
m=4#no of  genes

## true values for simulation
tsigma=1.2^2	;tsigmas=.7^2
tbet=array()
cat("enter true value of fixed effects total fixed effects are",p+t+m-2)
tbet=c(2.5,	.4,	1.06,	.26,	.32,	.50,.70,.60	)
truest=c(tbet,tsigma,tsigmas)
no=array()
cat("enter missing probabilities in each sequence, and total seq=",s)
#prob=c(.15,.15,.15)
#prob=scan(nmax=s) 
n=pmn=n1=list()## no of subjects in different sequences
no=c(10,10,10)
for (i in 1:s) {
  cat("no of subjects in seq",i,"\n")
  # no[i]=as.integer(readline(prompt="Enter an integer: "))
  n[[i]]=diag(no[i])
  pmn[[i]]=diag(no[i]*p*m)
  n1[[i]]=rep(1,no[i])
}
cat("enter no of generated missing observations in sequence 1,2,..,",s)
#num= scan(nmax=s)
num=c(2000,2000,2000)
treatment=list()
treatment[[1]]=c(1,2,3)
treatment[[2]]=c(2,1,3)
treatment[[3]]=c(3,2,1)
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

### In order to perform 300 simulations we generate the complete data 
nrep=300
y=b=e=list()
D=mclapply(n,"*",tsigmas)
covmate=mclapply(pmn,'*',tsigma)
set.seed(11)## to reproduce the results
for (i in 1:s) {
  b[[i]]=mvrnorm(nrep,rep(0,no[i]),D[[i]])
  e[[i]]=mvrnorm(nrep,rep(0,p*m*no[i]),covmate[[i]])
  y[[i]]= matrix(ncol =  p*m*no[i],nrow = nrep)
  for (jj in 1:nrep) {
    y[[i]][jj,]= x[[i]]%*%tbet+ z[[i]]%*%b[[i]][jj,]+e[[i]][jj,]
  }
}

##****** incorporating missing*******
## estimates are saved in festi matrix
festi=matrix(nrow = nrep,ncol =12*(p+t+m))
misper=matrix(nrow=nrep,ncol=s)
## for loop to run 100 simulations one by one
for(sim in 1:nrep){
  ## initial estimates for each simulation
  beth=ibeth=tbet-.12 
  sigmah=isigmah=tsigma-.1
  sigmash=isigmash=tsigmas+.12
  misx=obsx=misz=obsz=countmis=countobs=obsy=misy=newy=list()
  for (i in 1:s) {
    newy[[i]]=array(dim  = p*m*no[i])#
    newy[[i]]=as.numeric(y[[i]][sim,])
  }
  
  ##missing mechanism is MAR
  fn=function(arr){
    sapply(1:length(arr), function(i) rbinom(1,1,exp(arr[i])/(1+exp(arr[i]))))
  }
  misindctr=lapply(1:s, function(i) rep(0,p*m*no[i]))
  
  phi0=0.1;phi1=-.71;phi2=.1
 # repeat{
  repeat{
    logipr=lapply(1:s, function(i) sapply(1:no[i], function(j) 
      c(0, sapply(sapply(2:p, function(k)
        ifelse(k==2,phi0+phi1*newy[[i]][((j-1)*p*m)+((k-2)*m)+1],
               phi0+phi1*newy[[i]][((j-1)*p*m)+((k-2)*m)+1]
               +phi2*newy[[i]][((j-1)*p*m)+((k-3)*m)+1])),fn))))
    fl=  lapply(logipr,colSums)
    misper[sim,]= 1-sapply(1:s, function(i ) length(which(fl[[i]]==0))/no[i])
    if(min(misper[sim,])>0&max(misper[sim,])<1) break
  }
  cat("missing in each sequence is",misper[sim,])
#  if(mean(misper[sim,])<.45) break
 # }
  subind=fl
  fnn=function(arr,j){
    if(arr[2]==1&arr[3]==0) c((((j-1)*p*m)+m+1):(((j-1)*p*m)+2*m))
    else if(arr[3]==1&arr[2]==0) c((((j-1)*p*m)+2*m+1):(((j-1)*p*m)+3*m))
    else if(arr[3]==1&arr[2]==1) c((((j-1)*p*m)+m+1):(((j-1)*p*m)+3*m))
  }
  
  indd= lapply(1:s, function(i) sapply(1:no[i],function(j) fnn(logipr[[i]][,j],j)))
  ind1=lapply(indd,unlist)   
  for (i in 1:s) {
    misindctr[[i]][ind1[[i]]]=1
  }
  
  for (i in 1:s) {
    countmis[[i]]=diag(table(misindctr[[i]])[2])
    countobs[[i]]=diag(table(misindctr[[i]])[1])
    misx[[i]]=subset(x[[i]],misindctr[[i]]==1)
    misz[[i]]=subset(z[[i]],misindctr[[i]]==1)
    obsx[[i]]=subset(x[[i]],misindctr[[i]]==0)
    obsz[[i]]=subset(z[[i]],misindctr[[i]]==0)
    obsy[[i]]=subset(newy[[i]],misindctr[[i]]==0)
    misy[[i]]=subset(newy[[i]],misindctr[[i]]==1)
  }
  
  iD=mclapply(n,"*",sigmash)
  iS=mcMap("+",mcMap("%*%",mcMap("%*%",z,iD),lapply(z,t)),mclapply(pmn,'*',sigmah))
  mumis=mclapply(misx, "%*%",beth)
  muobs=mclapply(obsx, "%*%",beth)
  S11=mcMap("+",mcMap("%*%",mcMap("%*%",misz,iD),lapply(misz, t)),mclapply(countmis, "*",sigmah))
  S12=mcMap("%*%",mcMap("%*%",misz,iD),lapply(obsz, t))
  S21=mcMap("%*%",mcMap("%*%",obsz,iD),lapply(misz, t))
  S22=mcMap("+",mcMap("%*%",mcMap("%*%",obsz,iD),lapply(obsz, t)),mclapply(countobs, "*",sigmah))
  meanymisgivobsy=mcMap("+",mumis,mcMap("%*%",mcMap("%*%", S12,mclapply(S22, solve)),mcMap("-",mclapply(obsy,matrix),muobs)))
  varymisgivobsy=mcMap("-",S11,mcMap("%*%",mcMap("%*%", S12,mclapply(S22, solve)),S21))
  ymisgivyobs=list()
  #number of missing oservations generated
  for (i in 1:s) {
    ymisgivyobs[[i]]=mvrnorm(n=num[i],meanymisgivobsy[[i]],varymisgivobsy[[i]])
  }
  estti=updateestimates(beth,sigmah,sigmash,x,newy,z)
  ####################final algorithmmmmmmmm
  kj=0
  repeat {
    kj=kj+1
    bethh=as.numeric(estti[1:(p+t+m-2)])
    sigmahh=as.numeric(estti[p+t+m-1])
    sigmashh=as.numeric(estti[p+t+m])
    estti=  updateestimates(bethh,sigmahh,sigmashh,x,newy,z)
    if(abs(max(estti-c(bethh,sigmahh,sigmashh)))<5*10**(-3))
      break
  }
  kj
  cbind(estti,c(tbet,tsigma,tsigmas))
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
  festipa=c(estti,fstderr,(estti-truest)/truest,fstderr^2+(estti-truest)^2)
  print(cbind(estti,fstderr))
  
  ##********* now estimate and stderr using complete data ***************
  no1=unlist(lapply(1:s, function(i) length(which(subind[[i]]==0))))
  make1= lapply(1:s, function(i) sapply(1:length(which(subind[[i]]==1|subind[[i]]==2)), function(j) 
    ((((which(subind[[i]]==1|subind[[i]]==2)[j])-1)*p*m)+1):  
      ((which(subind[[i]]==1|subind[[i]]==2)[j])*p*m)))
  make11=lapply(make1,as.vector)   
  nx=lapply(1:s, function(i) x[[i]][-make11[[i]],])
  nnewy=lapply(1:s, function(i) newy[[i]][-make11[[i]]])
  nz=lapply(1:s, function(i) z[[i]][-make11[[i]],(1:no1[i])])
  nn=pmn1=nn1=list()
  for (i in 1:s) {
    nn[[i]]=diag(no1[i])
    pmn1[[i]]=diag(no1[i]*p*m)
    nn1[[i]]=rep(1,no1[i])
  }
  festii=updateestimate1(ibeth,isigmah,isigmash,nx,nnewy,nz)
  zz=0
  #  further steps
  repeat {
    zz=zz+1
    beth=as.numeric(festii[1:(p+t+m-2)])
    sigmah=as.numeric(festii[p+t+m-1])
    sigmash=as.numeric(festii[p+t+m])
    festii=updateestimate1(beth,sigmah,sigmash,nx,nnewy,nz)
    if(max(abs(festii-c(beth,sigmah,sigmash)))<5*10**-4)
      break
  }
  c(festii)
  zz
  
  ## for standard error we use bootstrap method by taking 
  #100 bootstrap samples of same size as original sample
  ntimes=100; D=lapply(1:s, function(i) sigmash*diag(no1[i]))
  sig0i=lapply(1:s, function(i) solve((1/sigmah)*crossprod(nz[[i]])+solve(D[[i]])))
  
  ## to generate 2000 obs and we use the average as the final sample
  ngibs=2000
  bgivy1=lapply(1:s, function(i) mvrnorm(ngibs,sig0i[[i]]%*%t(nz[[i]])%*%(nnewy[[i]]-nx[[i]]%*%beth)/sigmah, sig0i[[i]]) )
  bgivy= lapply(1:s, function(i) colMeans(bgivy1[[i]]))
  
  ## now estimate the parameters as if one has complete data
  beti=t(sapply(1:ntimes, function(k) solve(Reduce("+", lapply(1:s, function(i) crossprod(nx[[i]]))))%*%Reduce("+",lapply(1:s, function(i) (
    t(nx[[i]])%*%(nnewy[[i]]-nz[[i]]%*%bgivy[[i]]))))))
  sigmai=sapply(1:ntimes, function(k) Reduce("+",lapply(1:s, function(i) (1/(p*m*no1[i]))*crossprod(nnewy[[i]]-nx[[i]]%*%beth-nz[[i]]%*%bgivy[[i]]))))
  sigmasi=sapply(1:ntimes, function(k) Reduce("+",lapply(1:s, function(i) (1/no[i])*crossprod(bgivy[[i]]))))
  varbeti= t(sapply(1:ntimes, function(k) diag(sigmai[k]* solve(Reduce("+", lapply(1:s, function(i) crossprod(nx[[i]])))))))
  varsigmai=sapply(1:ntimes, function(i) (2*sigmai[i]**2)/(p*m*sum(no1)))
  varsigmasi=sapply(1:ntimes, function(i) (2*sigmasi[i]**2)/(sum(no1)))
  fvars=c(colMeans(varbeti),mean(varsigmai),mean(varsigmasi))+(1+(1/ntimes))*(apply(cbind(beti,sigmai,sigmasi),2,var))
  festic=c(festii,sqrt(fvars),(festii-truest)/truest,fvars+(festii-truest)^2)
  
  ### imputation analysis
  no1=no#unlist(lapply(1:s, function(i) length(which(subind[[i]]==0))))
  farr=lapply(1:s, function(i) numeric(length=p*m*no[i]))
  for (i in 1:s) {
    for (j in 1:no[i]) {
      if(length(which(logipr[[i]][,j]==1))==1) {
        farr[[i]][((which(logipr[[i]][,j]==1)-1)*m+1+((j-1)*p*m)):(((j-1)*p*m)+(which(logipr[[i]][,j]==1))*m)]= 1 } 
      
      if(length(which(logipr[[i]][,j]==1))==2) {
        farr[[i]][(((j-1)*p*m)+m+1):((j)*p*m)]=1}
    }
    
  }
  make11=farr
  
  ##response for mice package
  nnewy1=lapply(1:s, function(i) replace(newy[[i]],which(make11[[i]]==1),NA))
  respos=c(nnewy1[[1]],nnewy1[[2]],nnewy1[[3]])
  #no=c(10,10,10)
  ## fit model with original or log responses
  treatment1=c(rep(rep(c("trt1","trt2","trt3"),each=m),no[1]),
               rep(rep(c("trt2","trt1","trt3"),each=m),no[2]),
               rep(rep(c("trt3","trt2","trt1"),each=m),no[3]))
  subject=rep(paste("sub",1:sum(no)),each=(p*m))#+subject*gene
  Gene=rep(paste("gene",1:m),p*sum(no))
  Period=rep(rep(paste("per",1:p),each=m),sum(no))
  dat11=data.frame(respos,Period,Gene,treatment1,subject)
  ##mice imputation MI
  imputed_Data <- mice(dat11, m=5, maxit = 50, method = 'pmm', seed = 500)
  summary(imputed_Data)
  completeData <- rowMeans(sapply(1:5, function(i) complete(imputed_Data,i)$respos))
  #length(completeData[,1])
  nnewy2=list(completeData[1:(p*m*no[1])],completeData[((p*m*no[1])+1):(p*m*(no[1]+no[2]))],completeData[(p*m*(no[1]+no[2])+1):(p*m*sum(no))])
  nn=pmn1=nn1=list()
  for (i in 1:s) {
    nn[[i]]=diag(no[i])
    pmn1[[i]]=diag(length(nnewy2[[i]]))
    nn1[[i]]=rep(1,no[i])
  }
  festii=updateestimate1(ibeth,isigmah,isigmash,x,nnewy2,z)
  zz=0
  #  further steps
  repeat {
    zz=zz+1
    beth=as.numeric(festii[1:(p+t+m-2)])
    sigmah=as.numeric(festii[p+t+m-1])
    sigmash=as.numeric(festii[p+t+m])
    festii=updateestimate1(beth,sigmah,sigmash,x,nnewy2,z)
    if(max(abs(festii-c(beth,sigmah,sigmash)))<5*10**-4)
      break
  }
  c(festii)
  zz
  
  ## for standard error we use bootstrap method by taking 
  #100 bootstrap samples of same size as original sample
  ntimes=100; D=lapply(1:s, function(i) sigmash*diag(no[i]))
  sig0i=lapply(1:s, function(i) solve((1/sigmah)*crossprod(z[[i]])+solve(D[[i]])))
  
  ## to generate 1000 obs and we use the average as the final sample
  ngibs=2000
  bgivy1=lapply(1:s, function(i) mvrnorm(ngibs,sig0i[[i]]%*%t(z[[i]])%*%(nnewy2[[i]]-x[[i]]%*%beth)/sigmah, sig0i[[i]]) )
  bgivy= lapply(1:s, function(i) colMeans(bgivy1[[i]]))
  
  ## now estimate the parameters as if one has complete data
  beti=t(sapply(1:ntimes, function(k) solve(Reduce("+", lapply(1:s, function(i) crossprod(x[[i]]))))%*%Reduce("+",lapply(1:s, function(i) (
    t(x[[i]])%*%(nnewy2[[i]]-z[[i]]%*%bgivy[[i]]))))))
  sigmai=sapply(1:ntimes, function(k) Reduce("+",lapply(1:s, function(i) (1/(p*m*no[i]))*crossprod(nnewy2[[i]]-x[[i]]%*%beth-z[[i]]%*%bgivy[[i]]))))
  sigmasi=sapply(1:ntimes, function(k) Reduce("+",lapply(1:s, function(i) (1/no[i])*crossprod(bgivy[[i]]))))
  varbeti= t(sapply(1:ntimes, function(k) diag(sigmai[k]* solve(Reduce("+", lapply(1:s, function(i) crossprod(x[[i]])))))))
  varsigmai=sapply(1:ntimes, function(i) (2*sigmai[i]**2)/(p*m*sum(no)))
  varsigmasi=sapply(1:ntimes, function(i) (2*sigmasi[i]**2)/(sum(no)))
  fvars=c(colMeans(varbeti),mean(varsigmai),mean(varsigmasi))+(1+(1/ntimes))*(apply(cbind(beti,sigmai,sigmasi),2,var))
  festiim=c(festii,sqrt(fvars),(festii-truest)/truest,fvars+(festii-truest)^2)
  festi[sim,] =c(festipa,festic,festiim)
}  
mean(colMeans(misper[1:(sim-1),]))## on an average missing subjects in each sequence
aa=matrix(colMeans(festi[1:(sim-1),]),nrow=length(truest),ncol=12)
write.csv(aa,"D:\ ree_n1021mis.csv")
write.csv(festi,"D:\ all simulation,n 10,missing 24.csv")


