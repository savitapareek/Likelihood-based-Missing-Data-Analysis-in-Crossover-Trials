## clean the environment
rm(list=ls())
## install the relevant libraries
library(dplyr)
library(ggplot2)

#library(broom)## for augment function
## read the data
dat=read.csv(file.choose())## read anova 3 way all seq from csv files from vol D
head(dat)
subb=paste("Subject",1:12)
Subject=factor(rep(paste("Subject",1:12),each=30))#,labels=subb)
#sub1=rep(1:12,each=30)
Sequence=factor(rep(c("Sequence 1","Sequence 2","Sequence 3"),each=120))
Treatment=factor(c(rep(rep(c("10mg","Placebo","25mg"),each=10),4),
                   rep(rep(c("25mg","10mg","Placebo"),each=10),4),
                   rep(rep(c("Placebo","25mg","10mg"),each=10),4)))
Gene=factor(rep(paste("Gene",1:10),36))#,labels = c(paste("Gene",1:9),"Gene 10"))
Period=factor(rep(rep(paste("Period",1:3),each=10),12))
Response=dat$Response_allseq
## data frame having all the variables period, treatment gene and responses
data=data.frame(Response,Subject,Sequence,Treatment,Gene,Period)
attach(data)
head(data)
View(data)
str(data)

df=data
## subject by gene, period by gene and treatment by gene interaction as in normal paper
dff=df[1:120,];dff1=df[121:240,];dff2=df[241:360,]
head(dff)
#View(dff)
g1=ggplot(data = dff, mapping = aes(x =Period, y = Response,group=Gene,linetype=Gene)) +
  stat_summary(fun = mean, geom = "line")+
  theme_classic2()+
  theme(axis.title=element_text(size=7))+
  font("y.text", size = 7)+ font("x.text", size = 7)+
  theme(legend.title = element_text(size=7), #change legend title font size
        legend.text = element_text(size=7))

g2=ggplot(data = dff1, mapping = aes(x =Period, y = Response,group=Gene,linetype=Gene)) +
  stat_summary(fun = mean, geom = "line")+
  theme_classic2()+
  theme(axis.title=element_text(size=7))+
  font("y.text", size = 7)+ font("x.text", size = 7)+
  theme(legend.title = element_text(size=7), #change legend title font size
        legend.text = element_text(size=7))
g3=ggplot(data = dff2, mapping = aes(x =Period, y = Response,group=Gene,linetype=Gene)) +
  stat_summary(fun = mean, geom = "line")+
  theme_classic2()+
  theme(axis.title=element_text(size=7))+
  font("y.text", size = 7)+ font("x.text", size = 7)+
  theme(legend.title = element_text(size=7), #change legend title font size
        legend.text = element_text(size=7))

g4=ggplot(data = dff, mapping = aes(x =Treatment, y = Response,group=Gene,linetype=Gene)) +
  stat_summary(fun = mean, geom = "line")+
  theme_classic2()+
  theme(axis.title=element_text(size=7))+
  font("y.text", size = 7)+ font("x.text", size = 7)+
  theme(legend.title = element_text(size=7), #change legend title font size
        legend.text = element_text(size=7))
g5=ggplot(data = dff1, mapping = aes(x =Treatment, y = Response,group=Gene,linetype=Gene)) +
  stat_summary(fun = mean, geom = "line")+
  theme_classic2()+
  theme(axis.title=element_text(size=7))+
  font("y.text", size = 7)+ font("x.text", size = 7)+
  theme(legend.title = element_text(size=7), #change legend title font size
        legend.text = element_text(size=7))
g6=ggplot(data = dff2, mapping = aes(x =Treatment, y = Response,group=Gene,linetype=Gene)) +
  stat_summary(fun = mean, geom = "line")+
  theme_classic2()+
  theme(axis.title=element_text(size=7))+
  font("y.text", size = 7)+ font("x.text", size = 7)+
  theme(legend.title = element_text(size=7), #change legend title font size
        legend.text = element_text(size=7))

ggarrange(g1+rremove("xlab"),g4+rremove("ylab")+rremove("xlab"),
          g2+rremove("xlab"),g5+rremove("ylab")+rremove("xlab"),g3,g6+rremove("ylab"),
          nrow=3,ncol=2,labels = c("A","D","B","E","C","F"),
          legend="bottom",common.legend = T,font.label=list(color="black",size=9))

Subjects=factor(rep(paste("Sub",1:12),each=30))
ndf=data.frame(df,Subjects)
dff=ndf[1:120,];dff1=ndf[121:240,];dff2=ndf[241:360,]
g7=ggplot(data = dff, mapping = aes(x =Subjects, y = Response,group=Gene,linetype=Gene)) +
  stat_summary(fun = mean, geom = "line")+
  theme_classic2()+
  theme(axis.title=element_text(size=8))+
  font("y.text", size = 7)+font("x.text", size = 7)+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=9))
g8=ggplot(data = dff1, mapping = aes(x =Subjects, y = Response,group=Gene,linetype=Gene)) +
  stat_summary(fun = mean, geom = "line")+
  theme_classic2()+
  theme(axis.title=element_text(size=8))+font("x.text", size = 7)+
  font("y.text", size = 7)+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=9))
g9=ggplot(data = dff2, mapping = aes(x =Subjects, y = Response,group=Gene,linetype=Gene)) +
  stat_summary(fun = mean, geom = "line")+
  theme_classic2()+
  theme(axis.title=element_text(size=8))+
  font("y.text", size = 7)+font("x.text", size = 7)+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=9))
ggarrange(g7,g8+rremove("ylab"),g9+rremove("ylab"),
          nrow=1,ncol=3,labels = c("A","B","C"),
          legend="bottom",common.legend = T,font.label=list(color="black",size=9))

