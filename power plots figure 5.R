rm(list=ls())
library(ggplot2)
library(tidyverse)
library(ggpubr)
### for treatment testing H1 =.234
x=c(10,20,50)
#yc=c(0.720, 
#     0.962, 
#     1.000 
 #    )
All_reponses=c(0.61,	0.87,	1 )
Pairwise_res_1_2=c(0.63,	0.87,	0.99)
Pairwise_res_1_3=c(0.49,	0.78,	1)
Pairwise_res_2_3=c(0.63	,0.89,	1)
Pairwise_trt_1_2=c(	0.18,	0.32	,0.61	)
#ymis45=c(	.672,.927,.994)
# plot(x,yc,type="l",xlim=c(25,150),ylim=c(0,1),xlab="Sample Size",ylab="Power",lty=1,col="cyan")
# lines(x,ymis15,lty=2,col="red")
# lines(x,ymis25,lty=3,col="blue")
# lines(x,ymis35,lty=4,col="purple")
# lines(x,ymis45,lty=5,col="green")
# legend("topright",legend=c("complete data","15% missing","25% missing","35% missing","45% missing"), col=c("cyan","red","blue","purple","green"),
#        lty =c(1,2,3,4,5),seg.len = 1, ncol=1,cex=.71)


df2 <- data.frame(Hypothesis=rep(c( "All repsonses (i)",
                                    "Pairwise responses 1, 2 (ii)",
                                    "Pairwise responses 1, 3 (ii)",
                                    "Pairwise responses 2, 3 (ii)",
                                    "Pairwise treatments 1, 2 (iii)"), each=length(x)),
                  Total_subjects=rep(x,5),
                  Power=c(All_reponses,Pairwise_res_1_2,Pairwise_res_1_3,Pairwise_res_2_3,Pairwise_trt_1_2))
head(df2)
p1=ggplot(df2, aes(Total_subjects,Power,group=Hypothesis)) +
  geom_line(aes(linetype=Hypothesis))+
 geom_point()
 # scale_color_brewer(palette="Paired")
#p + theme(legend.position="bottom")
p1=p1+theme_classic()
p1

getwd()
ggexport(p1, filename = "figure_hyp.png")
