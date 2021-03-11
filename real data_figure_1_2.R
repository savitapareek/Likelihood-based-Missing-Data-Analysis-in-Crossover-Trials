## cleaning the enviornment
rm(list=ls())
gc()

## calling the required libraries
library(ggpubr)
library(ggplot2)

### complete case analysis
##each sequence has 4 subjects removing missing subject completely
data11=read.csv(file.choose(),header = T)#read 'anova 3 way allseq.csv' file
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

### for interaction plot of each sequence
## seq 1
dat=data11[c(1:120),]
subject=as.factor(paste("sub",rep(1:4,each=30),sep=""))
Period=dat$period;Treatment=dat$treatment;Gene=dat$gene; mean_response=dat$Response_allseq
df2=data.frame(Period,Treatment,Gene,mean_response,subject)
p11=ggplot(df2, aes(Period,mean_response,group=Gene,linetype=Gene)) +
  stat_summary(fun.y = mean, geom = "line")+
  stat_summary(fun.y = mean, geom = "point",lwd=.2)+theme_classic()

p12=ggplot(df2, aes(Treatment,mean_response,group=Gene,linetype=Gene)) +
  stat_summary(fun.y = mean, geom = "point",lwd=.2)+
  stat_summary(fun.y = mean, geom = "line")+ theme_classic()

p13=ggplot(df2, aes(subject,mean_response,group=Gene,linetype=Gene)) +
  stat_summary(fun.y = mean, geom = "line")+
  stat_summary(fun.y = mean, geom = "point",lwd=.2)+theme_classic()

## seq 2
dat=data11[c(121:240),]
Period=dat$period;Treatment=dat$treatment;Gene=dat$gene; mean_response=dat$Response_allseq
df2=data.frame(Period,Treatment,Gene,mean_response,subject)
p21=ggplot(df2, aes(Period,mean_response,group=Gene,linetype=Gene)) +
  stat_summary(fun.y = mean, geom = "point",lwd=.2)+
  stat_summary(fun.y = mean, geom = "line")+ theme_classic()

p22=ggplot(df2, aes(Treatment,mean_response,group=Gene,linetype=Gene)) +
  stat_summary(fun.y = mean, geom = "point",lwd=.2)+
  stat_summary(fun.y = mean, geom = "line")+ theme_classic()

p23=ggplot(df2, aes(subject,mean_response,group=Gene,linetype=Gene)) +
  stat_summary(fun.y = mean, geom = "point",lwd=.2)+
  stat_summary(fun.y = mean, geom = "line")+ theme_classic()

## seq 3
dat=data11[-c(1:240),]
Period=dat$period;Treatment=dat$treatment;Gene=dat$gene; mean_response=dat$Response_allseq
df2=data.frame(Period,Treatment,Gene,mean_response,subject)
p31=ggplot(df2, aes(Period,mean_response,group=Gene,linetype=Gene)) +
  stat_summary(fun.y = mean, geom = "point",lwd=.2)+
  stat_summary(fun.y = mean, geom = "line")+ theme_classic()

p32=ggplot(df2, aes(Treatment,mean_response,group=Gene,linetype=Gene)) +
  stat_summary(fun.y = mean, geom = "point",lwd=.2)+
  stat_summary(fun.y = mean, geom = "line")+ theme_classic()

p33=ggplot(df2, aes(subject,mean_response,group=Gene,linetype=Gene)) +
  stat_summary(fun.y = mean, geom = "point",lwd=.2)+
  stat_summary(fun.y = mean, geom = "line")+
  theme_classic()                          

fig=ggarrange(p13,p23,p33,common.legend = T,ncol=3,legend = "bottom",labels = c("A","B","C"))
fig
## combine all three sequence now
fp=ggarrange(p11,p12,p21,p22,p31,p32,labels = c("A","D","B","E","C","F"),nrow = 3,ncol=2,common.legend = T,legend = "bottom")
fp

ggexport(fp, filename = "interaction_allseq.png")
ggexport(fig, filename = "interacsub_allseq.png")
