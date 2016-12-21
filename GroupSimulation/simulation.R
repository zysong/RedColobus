library(dplyr)
library(lattice)
source("groupSizeSimulation.R")
df.fbd<-gsSimulation() #female-biased dispersal
df.mbd<-gsSimulation(d.m=0.1, d.f=0.001) #male-biased dispersal
df.fbd.wv<-gsSimulation(p=0.025, v.m=.1, v.f=.02) #female-biased dispersal + weak virulence
df.mbd.wv<-gsSimulation(p=0.025, v.m=.1, v.f=.02, d.m = 0.1, d.f = 0.001) #male-biased dispersal + weak virulence
df.fbd["dispersal"]<-rep("fbd", 50)
df.fbd["virulence"]<-rep("high", 50)
df.mbd["dispersal"]<-rep("mbd", 50)
df.mbd["virulence"]<-rep("high", 50)
df.fbd.wv["dispersal"]<-rep("fbd", 50)
df.fbd.wv["virulence"]<-rep("low", 50)
df.mbd.wv["dispersal"]<-rep("mbd", 50)
df.mbd.wv["virulence"]<-rep("low", 50)
df.summary<-rbind(df.fbd, df.mbd, df.fbd.wv, df.mbd.wv)
df.summary$dispersal<-as.factor(df.summary$dispersal)
df.summary$virulence<-as.factor(df.summary$virulence)

boxplot(mean.g~dispersal*virulence, data=df.summary, ylab="Mean group size")
boxplot(var.g.w~dispersal*virulence, data=df.summary, ylab="Within-group variance")
boxplot(var.g.b~dispersal*virulence, data=df.summary, ylab="Between-group variance")
boxplot(cov.g.w~dispersal*virulence, data=df.summary, ylab="Within-group CoV")
boxplot(cov.g.b~dispersal*virulence, data=df.summary, ylab="Between-group CoV")
