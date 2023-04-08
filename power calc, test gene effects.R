## cleaning the enviornment
rm(list=ls())
gc()

## install required packages
library(psych) ## for trace of a matrix
library(doParallel) ## for parallelization
library(MASS) ## for mvtnorm

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

##########---- cretaing matrices and data generation
strt=Sys.time()
s=3#no of sequences
p=t=3#no of periods and no of treatments
m=4#no of  genes

## true values for simulation
tsigma=1.2^2	;tsigmas=.7^2
tbet=array()
cat("enter true value of fixed effects total fixed effects are",p+t+m-2)
tbet=c(2.5,	.4,	1.06,	.26,	.32,	.50,.7,.6	)
truest=c(tbet,tsigma,tsigmas)
no=array()
cat("enter missing probabilities in each sequence, and total seq=",s)
n=pmn=n1=list()## no of subjects in different sequences
no=c(3,3,4)
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
treatment[[1]]=c(1,2,3)#,3,4)
treatment[[2]]=c(2,1,3)#,3,4)
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

### In order to perform 500 simulations we generate the complete data 
nrep=1000
y=D=b=e=list()
D=mclapply(n,"*",tsigmas)
covmate=mclapply(pmn,'*',tsigma)
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
## for loop to run 500 simulations one by one
cl1=makeCluster(detectCores()-1)
registerDoParallel(cl1)
festi=foreach(sim = 1:nrep,.combine = rbind,.packages = c("psych","parallel","MASS"))%dopar%{
  ## initial estimates for each simulation
  co=0
  beth=ibeth=tbet-.12 
  sigmah=isigmah=tsigma-.1
  sigmash=isigmash=tsigmas+.12
  misindctr=misx=obsx=misz=obsz=countmis=countobs=obsy=misy=w=newy=D=S=list()
  for (i in 1:s) {
    newy[[i]]=array(dim  = p*m*no[i])#
    newy[[i]]=as.numeric(y[[i]][sim,])
  }
  
  ##missing mechanism is MAR
  fn=function(arr){
    sapply(1:length(arr), function(i) rbinom(1,1,exp(arr[i])/(1+exp(arr[i]))))
  }
  misindctr=lapply(1:s, function(i) rep(0,p*m*no[i]))
  
  phi0=0.1;phi1=-.41;phi2=.1
  repeat{
    logipr=lapply(1:s, function(i) sapply(1:no[i], function(j) 
      c(0, sapply(sapply(2:p, function(k)
        ifelse(k==2,phi0+phi1*newy[[i]][((j-1)*p*m)+((k-2)*m)+1],
               phi0+phi1*newy[[i]][((j-1)*p*m)+((k-2)*m)+1]
               +phi2*newy[[i]][((j-1)*p*m)+((k-3)*m)+1])),fn))))
    fl=  lapply(logipr,colSums)
    misper= 1-sapply(1:s, function(i ) length(which(fl[[i]]==0))/no[i])
    if(min(misper)>0&max(misper)<1) break
  }
  cat("missing in each sequence is",misper)
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
  #number of missing observations generated
  for (i in 1:s) {
    ymisgivyobs[[i]]=mvrnorm(n=num[i],meanymisgivobsy[[i]],varymisgivobsy[[i]])
  }
  esttif=updateestimates(beth,sigmah,sigmash,x,newy,z)
  ####################final algorithmmmmmmmm
  kj=0
  repeat {
    kj=kj+1
    bethh=as.numeric(esttif[1:(p+t+m-2)])
    sigmahh=as.numeric(esttif[p+t+m-1])
    sigmashh=as.numeric(esttif[p+t+m])
    esttif=  updateestimates(bethh,sigmahh,sigmashh,x,newy,z)#ymisgivyobs)
    if(max(abs(bethh-esttif[1:(p+t+m-2)]))<5*10**-4&abs(sigmahh-esttif[p+t+m-1] ) <5*10**-4&abs(sigmashh-esttif[p+t+m] ) <5*10**-4 )
      break
  }
  kj
  cbind(esttif,c(tbet,tsigma,tsigmas))

    ## likelihood under H1
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
  ymisgivyobs=lapply(1:s, function(i) mvrnorm(n=num[i],meanymisgivobsy[[i]],varymisgivobsy[[i]]))
  ymisobs=lapply(1:s, function(i) t(sapply(1:num[i],function(j) (
    replace(newy[[i]], which(misindctr[[i]]==1),  ymisgivyobs[[i]][j,])
  )) ))
  b0i= lapply(1:s, function(i) t(sapply(1:num[i],function(j) (
    sig0i[[i]]%*%t(z[[i]])%*%(ymisobs[[i]][j,]-x[[i]]%*%beth)/sigmae
  )) ))
  L1=Reduce("+", lapply(1:s, function(i) (
    -.5*(p*m*no[i]*log(sigmae)+(1/sigmae)*(tr(crossprod(z[[i]]%*%sig0i[[i]]))
                                           +mean(sapply(1:num[i], function(j) 
                                             crossprod(ymisobs[[i]][j,]-x[[i]]%*%beth
                                                       -z[[i]]%*%b0i[[i]][j,])) ))
         +no[i]*log(sigmash)+tr(solve(iD[[i]])%*%sig0i[[i]])
         +mean(sapply(1:num[i], function(k) t(b0i[[i]][k,])%*%solve(iD[[i]])%*%b0i[[i]][k,] )))
  ))) 
 
  ####estimation under H0
  # make a note: for reduced model: to test response varaite 1,2, 3
  hh1=c(7,8)
  beth=tbet[-hh1]-.12
  sigmah=tsigma-.1
  sigmash=tsigmas+.2
  hh=2; df=2
misx1=obsx1=misz=obsz=countmis=countobs=obsy=misy=list()
  x1=lapply(1:s, function(i) x[[i]][,-hh1])
 
  for (i in 1:s) {
    countmis[[i]]=diag(table(misindctr[[i]])[2])
    countobs[[i]]=diag(table(misindctr[[i]])[1])
    misx1[[i]]=subset(x1[[i]],misindctr[[i]]==1)
    misz[[i]]=subset(z[[i]],misindctr[[i]]==1)
    obsx1[[i]]=subset(x1[[i]],misindctr[[i]]==0)
    obsz[[i]]=subset(z[[i]],misindctr[[i]]==0)
    obsy[[i]]=subset(newy[[i]],misindctr[[i]]==0)
    misy[[i]]=subset(newy[[i]],misindctr[[i]]==1)
  }
  iD=iS=list()
  iD=mclapply(n,"*",sigmash)
  iS=mcMap("+",mcMap("%*%",mcMap("%*%",z,iD),lapply(z,t)),mclapply(pmn,'*',sigmah))
  mumis=mclapply(misx1, "%*%",beth)
  muobs=mclapply(obsx1, "%*%",beth)
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
  estti=updateestimates(beth,sigmah,sigmash,x1,newy,z)
  ####################final algorithmmmmmmmm
  repeat {
    beth=as.numeric(estti[1:(p+t+m-2-hh)])# make a note
    sigmah=as.numeric(estti[p+t+m-hh-1])
    sigmash=as.numeric(estti[p+t+m-hh])
    estti=array()
    estti=  updateestimates(beth,sigmah,sigmash,x1,newy,z)
    if(max(abs(c(beth,sigmah,sigmash)-estti))<5*10**-4)
      break
  }
  
  ##likelihood under H0
  beth=estti[1:(p+t+m-2-hh)];sigmae=estti[p+t+m-hh-1];sigmash=estti[p+t+m-hh]
  iD=mclapply(n,"*",sigmash)
  sig0i=mclapply(mcMap("+",mclapply(mclapply(z, crossprod), "*",sigmae**-1),mclapply(iD, solve)),solve)
  
  mumis=mclapply(misx1, "%*%",beth)
  muobs=mclapply(obsx1, "%*%",beth)
  S11=mcMap("+",mcMap("%*%",mcMap("%*%",misz,iD),lapply(misz, t)),mclapply(countmis, "*",sigmae))
  S12=mcMap("%*%",mcMap("%*%",misz,iD),lapply(obsz, t))
  S21=mcMap("%*%",mcMap("%*%",obsz,iD),lapply(misz, t))
  S22=mcMap("+",mcMap("%*%",mcMap("%*%",obsz,iD),lapply(obsz, t)),mclapply(countobs, "*",sigmae))
  varymisgivobsy=meanymisgivobsy=ymisobs=list()
  meanymisgivobsy=mcMap("+",mumis,mcMap("%*%",mcMap("%*%", S12,mclapply(S22, solve)),mcMap("-",mclapply(obsy,matrix),muobs)))
  varymisgivobsy=mcMap("-",S11,mcMap("%*%",mcMap("%*%", S12,mclapply(S22, solve)),S21))
  ymisgivyobs=lapply(1:s, function(i) mvrnorm(n=num[i],meanymisgivobsy[[i]],varymisgivobsy[[i]]))
  ymisobs=lapply(1:s, function(i) t(sapply(1:num[i],function(j) (
    replace(newy[[i]], which(misindctr[[i]]==1),  ymisgivyobs[[i]][j,])
  )) ))
  b0i= lapply(1:s, function(i) t(sapply(1:num[i],function(j) (
    sig0i[[i]]%*%t(z[[i]])%*%(ymisobs[[i]][j,]-x1[[i]]%*%beth)/sigmae
  )) ))
  L0=Reduce("+", lapply(1:s, function(i) (
    -.5*(p*m*no[i]*log(sigmae)+(1/sigmae)*(tr(crossprod(z[[i]]%*%sig0i[[i]]))
                                           +mean(sapply(1:num[i], function(j) 
                                             crossprod(ymisobs[[i]][j,]-x1[[i]]%*%beth
                                                       -z[[i]]%*%b0i[[i]][j,])) ))
         +no[i]*log(sigmash)+tr(solve(iD[[i]])%*%sig0i[[i]])
         +mean(sapply(1:num[i], function(k) t(b0i[[i]][k,])%*%solve(iD[[i]])%*%b0i[[i]][k,] )))
  )))
  cr=qchisq(1-.05,df) ##is the degree of freedom 
  if(-2*(L0-L1)>cr) co=co+1
 c(co,misper)
} 
stopCluster(cl1)
colMeans(festi[,-1])## on an average missing subjects in each sequence
powr=as.numeric((sum(festi[,1]))/(nrep)) # writing power
round(powr,2)
Sys.time()-strt

