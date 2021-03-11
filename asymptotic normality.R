rm(list = ls())
library(ggplot2)
library(EnvStats)

## read the simulation results 
sim1=read.csv(file.choose(),header = T)

## store the results; one for 100 simulations and other 500 simulations
sim100=sim1[1:100,c(2:9)]
sim500=sim1[1:500,c(2:9)]

## first we will work on 100 simulations
esti=sim100

## standardize the estimates
trep=nrow(esti)
obs=matrix(nrow = 8,ncol=nrow(esti))
for (i in 1:8) {
  obs[i,]=(esti[,i]-mean(esti[,i]))/sqrt(var(esti[,i]))
}

## Find the empirical density corresponding to each parameter
x <- qemp(p = seq(0, 1, len = trep), obs[1,]) 
y <- demp(x, obs[1,]) 
x2 <- qemp(p = seq(0, 1, len = trep), obs[2,]) 
y2 <- demp(x2, obs[2,]) 
x3 <- qemp(p = seq(0, 1, len = trep), obs[3,]) 
y3 <- demp(x3, obs[3,]) 
x4 <- qemp(p = seq(0, 1, len = trep), obs[4,]) 
y4 <- demp(x4, obs[4,]) 
x5 <- qemp(p = seq(0, 1, len = trep), obs[5,]) 
y5 <- demp(x5, obs[5,]) 
x6 <- qemp(p = seq(0, 1, len = trep), obs[6,]) 
y6 <- demp(x6, obs[6,]) 
x7 <- qemp(p = seq(0, 1, len = trep), obs[7,]) 
y7 <- demp(x7, obs[7,]) 
x8 <- qemp(p = seq(0, 1, len = trep), obs[8,]) 
y8 <- demp(x8, obs[8,]) 

## data frame and ggplot graph
df33 <-data.frame(parameter=rep(c("Intercept","Period1","Trt1","Res1","Res2","Res3","sigma_e^2","sigma_s^2"), each=trep),
                 parameter_value=c(x,x2,x3,x4,x5,x6,x7,x8),
                 relative_frequency=c(y,y2,y3,y4,y5,y6,y7,y8))
p1=ggplot(df33, aes(parameter_value,relative_frequency,group=parameter)) +
  geom_line(aes(linetype=parameter),lwd=.09)+ theme_classic()
p1

