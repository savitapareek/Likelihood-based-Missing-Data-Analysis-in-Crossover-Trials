## cleaning the enviornment
rm(list=ls())
gc()

# installing required librraries
library(psych)
library(doParallel)
library(MASS)

strt=Sys.time()

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
s=2#no of sequences
p=t=2#no of periods and no of treatments
m=4#no of  genes
## true values for simulation
tsigma=1.2^2	;tsigmas=.7^2
cat("enter true value of fixed effects total fixed effects are",p+t+m-2)
tbet=c(4.5,	0.2,	1.06,	.46,	1.09,	.50	)
truest=c(tbet,tsigma,tsigmas)
cat("enter missing probabilities in each sequence, and total seq=",s)
prob=c(.15,.15)
n=pmn=n1=list()
no=c(5,5)## no of objects in different sequences
for (i in 1:s) {
  n[[i]]=diag(no[i])
  pmn[[i]]=diag(no[i]*p*m)
  n1[[i]]=rep(1,no[i])
}
cat("enter no of generated missing observations in sequence 1,2,..,",s)
num=c(2000,2000)
treatment=list()
treatment[[1]]=c(1,2)
treatment[[2]]=c(2,1)

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

## X matrix under reduced model: to test response varaite 1, 3
xr=lapply(1:s, function(i) x[[i]][,-c(4,6)])

## generating data from true values
nrep=1000
y=D=S=yr=list()
D=mclapply(n,"*",tsigmas)
covmate=mclapply(pmn,'*',tsigma)
b=e=list()
for (i in 1:s) {
  b[[i]]=mvrnorm(nrep,rep(0,no[i]),D[[i]])
  e[[i]]=mvrnorm(nrep,rep(0,p*m*no[i]),covmate[[i]])
  y[[i]]= matrix(ncol =  p*m*no[i],nrow = nrep)
  for (jj in 1:nrep) {
    y[[i]][jj,]= x[[i]]%*%tbet+ z[[i]]%*%b[[i]][jj,]+e[[i]][jj,]
  }
}

######*********** incorporating missing**********
prob=c(.15,.15)
cl1=makeCluster(detectCores()-1)
registerDoParallel(cl1)
festi=foreach(sim= 1:nrep,.combine = rbind,.packages = c('psych','parallel','MASS'))%dopar% {
  co=0
  newy=lapply(1:s,function(i) y[[i]][sim,])
  #estimation under H0
   # make a note: for reduced model: to test response varaite 1, 3
  beth=tbet[-c(4,6)]-.2
  sigmah=tsigma-.1
  sigmash=tsigmas+.2
  missubindctr=misindctr=misx1=obsx1=misz=obsz=countmis=countobs=obsy=misy=w=missub= per=D=S=list()
  x1=xr
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
    misx1[[i]]=subset(x1[[i]],misindctr[[i]]==1)
    misz[[i]]=subset(z[[i]],misindctr[[i]]==1)
    obsx1[[i]]=subset(x1[[i]],misindctr[[i]]==0)
    obsz[[i]]=subset(z[[i]],misindctr[[i]]==0)
    obsy[[i]]=subset(newy[[i]],misindctr[[i]]==0)
    misy[[i]]=subset(newy[[i]],misindctr[[i]]==1)
  }
  missub=mclapply(mclapply(w,"/",p), ceiling)
  per=mclapply(w, "%%",p)
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
  estti=updateestimates(beth,sigmah,sigmash,newy,x1,z,ymisgivyobs)
  ####################final algorithmmmmmmmm
  repeat {
    beth=as.numeric(estti[1:(p+t+m-4)])# make a note
    sigmah=as.numeric(estti[p+t+m-3])
    sigmash=as.numeric(estti[p+t+m-2])
    estti=array()
    estti=  updateestimates(beth,sigmah,sigmash,newy,x1,z,ymisgivyobs)
    if(max(abs(c(beth,sigmah,sigmash)-estti))<5*10**-4)
      break
  }
  
  ## estimate under H1 when data is generated from model
  beth=tbet-.2 # make a note
  sigmah=tsigma-.1
  sigmash=tsigmas+.2;misx=obsx=list()
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
  esttif=updateestimates(beth,sigmah,sigmash,newy,x,z,ymisgivyobs)
  ####################final algorithmmmmmmmm
  repeat {
    beth=as.numeric(esttif[1:(p+t+m-2)])# make a note
    sigmah=as.numeric(esttif[p+t+m-1])
    sigmash=as.numeric(esttif[p+t+m])
    esttif=  updateestimates(beth,sigmah,sigmash,newy,x,z,ymisgivyobs)
    if(max(abs(c(beth,sigmah,sigmash)-esttif))<5*10**-4)
      break
  }
  
  ## likelihood under H0 and under H1
  #under H1
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
  
  ##under H0
  beth=estti[1:(p+t+m-4)];sigmae=estti[p+t+m-3];sigmash=estti[p+t+m-2]
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
  cr=qchisq(1-.05,2) ##2 is the degree of freedom 
  if(-2*(L0-L1)>cr) co=co+1
  c(esttif,co)
}  
stopCluster(cl1)
# writing power
powr=as.numeric((colSums(festi)[(p+t+m)+1])/nrep) 
powr
Sys.time()-strt
