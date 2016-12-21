library(dplyr)
library(forecast)

InGroupData<-read.csv("../Data/InGroupData.csv", check.names = FALSE)
AgeData<-read.csv("../Data/AgeData.csv", check.names = FALSE)
VirusData<-read.csv("../Data/VirusData.csv", check.names = FALSE)
VirusData<-subset(VirusData, select=-Age)
VirusData<-na.omit(VirusData)

#The number of immigrants, emigrants, and deaths every year by Sex
sum.immig<-tapply(InGroupData$Immig, list(InGroupData$Immig, InGroupData$Sex), length)
sum.immig<-sum.immig[-nrow(sum.immig),] #remove incomplete 2015 data 
sum.emig<-tapply(InGroupData$Emig, list(InGroupData$Emig, InGroupData$Sex), length)
sum.death<-tapply(InGroupData$Death, list(InGroupData$Death, InGroupData$Sex), length)
sum.death<-sum.emig[-nrow(sum.death),] #remove incomplete 2015 data 

sums.InGroup<-colSums(subset(InGroupData, select = as.character(2006:2014)))
sex.ratio<-colSums(subset(InGroupData, Sex=="F", select = as.character(2006:2014)))/
  colSums(subset(InGroupData, Sex=="M", select = as.character(2006:2014)))
mean.age<-colSums(subset(AgeData, select = as.character(2006:2014)))/sums.InGroup
AdultMat.f<-subset(AgeData, Sex=="F", select = as.character(2006:2014))>3
AdultMat.m<-subset(AgeData, Sex=="M", select = as.character(2006:2014))>4
InfantMat<-AgeData[2:10]==1
JuvenileMat<-AgeData[2:10]==2 | AgeData[2:10]==3 
sums.adult.f<-colSums(subset(InGroupData, Sex=="F", select = as.character(2006:2014))*AdultMat.f)
sums.adult.m<-colSums(subset(InGroupData, Sex=="M", select = as.character(2006:2014))*AdultMat.m)
sums.infant<-colSums(InGroupData[5:13]*InfantMat)
sums.juvenile<-colSums(InGroupData[5:13]*JuvenileMat)
ts.demo<-ts(cbind(Sum= sums.InGroup, AdF=sums.adult.f, AdM=sums.adult.m, Juv=sums.juvenile, Ift=sums.infant), start = 2006)
plot(ts.demo)

#Compute ratios
sex.ratio.adult<-sums.adult.f/sums.adult.m
sex.ratio.adult.m<-1/sex.ratio.adult
infant.ratio<-sums.infant/sums.adult.f #birth rate
juvenile.ratio<-sums.juvenile/sums.adult.f
ts.ratio<-ts(cbind(AdF=sex.ratio.adult, AdM=sex.ratio.adult.m, Juv=juvenile.ratio, Ift=infant.ratio), start = 2006)
plot(ts.ratio)

#population characteristics is likely to have lagged effect on life events like migration and birth
ts.demo.lag1<-lag(ts.demo)
ts.ratio.lag1<-lag(ts.ratio)
sums.adult.m.lag1<-sums.adult.m[-length(sums.adult.m)]
sums.adult.f.lag1<-sums.adult.f[-length(sums.adult.f)]
sums.juvenile.lag1<-sums.juvenile[-length(sums.juvenile)]
sex.ratio.adult.lag1<-sex.ratio.adult[-length(sex.ratio.adult)]
dif.m.lag1<-diff(sums.adult.m)
dif.f.lag1<-diff(sums.adult.f)
#ndiffs(ts.demo)
ts.dif1<-diff(ts.demo, 1)
Acf(ts.demo)

#birth rate is negatively correlated with the population size
pop.lag1<-sums.adult.f.lag1+sums.adult.m.lag1+sums.juvenile.lag1
sums.adult.lag1<-sums.adult.f.lag1+sums.adult.m.lag1
lm.birth<-lm(infant.ratio[-1]~pop.lag1)
summary(lm.birth)
summary(lm(infant.ratio[-1]~sums.adult.lag1+I(sums.adult.lag1^2)))
summary(lm(infant.ratio[-1]~pop.lag1+I(pop.lag1^2)))
pdf("../Manuscript/figures/Birth.pdf")
plot(pop.lag1, infant.ratio[-1], xlab = "Population size", ylab = "Birth rate")
abline(lm.birth)
dev.off()

#adult female change as a response to the numbers of adult males and adult females in the previous year
summary(lm(dif.f.lag1~sums.adult.m.lag1+sums.adult.f.lag1)) #both are significant
summary(lm(dif.f.lag1~sums.adult.m[-1]+sums.adult.f[-1])) #neither is significant
summary(lm(sums.adult.f[-1]~sums.adult.m.lag1+sums.adult.f.lag1)) #only correlated with the males

#immigration of females ~ males + males/females + interaction
summary(lm(sum.immig[,'F']~sums.adult.m+sums.adult.f))
summary(lm(sum.immig[-1,'F']~pop.lag1))
summary(lm(sum.immig[-1,'F']~sums.adult.m.lag1+sums.adult.f.lag1)) #neither is significant
lm.immig<-lm(sum.immig[-1,'F']~sums.adult.m.lag1+sums.adult.f.lag1+I(sums.adult.m.lag1/sums.adult.f.lag1)+I(sums.adult.m.lag1^2/sums.adult.f.lag1))
summary(lm.immig)#interaction is significantly and positively correlated with immigration
step(lm.immig)
#This is the best model:
summary(lm(sum.immig[-1,'F']~sums.adult.m.lag1+I(sums.adult.m.lag1^2/sums.adult.f.lag1)))
#All are significant
lm.immig1<-lm(sum.immig[-1,'F']~I(sums.adult.m.lag1^2/sums.adult.f.lag1))
summary(lm.immig1)#Interestingly, this is not significant
summary(lm(sum.immig[-1,'F']~ sums.adult.m.lag1+I(sums.adult.m.lag1^2)))
pdf("../Manuscript/figures/Immig.pdf")
plot(sums.adult.m.lag1/sex.ratio.adult.lag1, sum.immig[-1,'F'], xlab = "M(t-1)^2/F(t-1)", ylab = "Female immigration")
abline(lm.immig1)
dev.off()

#emigration of females ~ females + females/males + interaction
summary(lm(sum.emig[,'F']~sums.adult.m.lag1+sums.adult.f.lag1))
summary(lm(sum.emig[,'F']/sums.adult.f.lag1~sex.ratio.adult.lag1))
summary(lm(sum.emig[,'F']~I(1/sums.adult.m.lag1)*sums.adult.f.lag1)) 
#correlated positively with females and negatively with males, uncorrelated with sex ratio
lm.emig<-lm(sum.emig[,'F']~sums.adult.f.lag1+I(sums.adult.f.lag1*sex.ratio.adult.lag1))
summary(lm.emig)
#The interaction is positively and significantly correlated with emigration, 
#but female number and sex ratio are only negatively correlated and marginally significant.
vif(lm.emig) #variance inflation factor
plot(hatvalues(lm.emig)) #leverage
lm.emig1<-lm(sum.emig[,'F']~sums.adult.f.lag1+ sums.adult.m.lag1+I(sums.adult.f.lag1*sex.ratio.adult.lag1))
summary(lm.emig1) #significant correlation
step(lm.emig1)
par(mfrow=c(2,2))
plot(lm.emig1)
par(mfrow=c(1,1))
pdf("../Manuscript/figures/Emig.pdf")
plot(sums.adult.f.lag1*sex.ratio.adult.lag1, sum.emig[,'F'], xlab = "F(t-1)^2/M(t-1)", ylab = "Female emigration")
abline(lm.emig1)
dev.off()

#plot results
pdf("../Manuscript/figures/GroupSize.pdf")
plot(startObs:endYear,sums.InGroup, type="l", lwd=3, ylim=c(0,150), 
     xlab = "Year", ylab = "Number of Individuals") #group size increases
lines(startObs:endYear,sums.adult.f, col="red", lwd=3)
lines(startObs:endYear,sums.adult.m, col="blue", lwd=3)
lines(startObs:endYear,sums.infant, col="green", lwd=3)
lines(startObs:endYear,sums.juvenile, col="orange", lwd=3)
legend(2006, 150, c("All", "FA", "MA", "Inf", "Juv"), 
       col = c("black", "red", "blue", "green", "orange"), lty = 1, lwd =3)
dev.off()

pdf("../Manuscript/figures/SexRatio.pdf")
plot(startObs:endYear, sex.ratio, type="l", lwd=3, ylim=c(1, 3.5),
     xlab = "Year", ylab = "Sex ratio (F/M)") #sex ratio decreases and then increases, but does not recover
lines(startObs:endYear, sex.ratio.adult, type="l", lwd=3, lty=2,
      xlab = "Year", ylab = "Adult sex ratio (F/M)")
legend(2013, 3.5, c("All", "Adult"), lty=c(1,2), lwd=3)
dev.off()

plot(startObs:endYear, mean.age)

pdf("../Manuscript/figures/BirthRate.pdf")
plot(startObs:(endYear-1), infant.ratio[1:(endYear-startObs)], type="l", lwd=3, xlab = "Year", ylab = "Birth rate")
dev.off()

plot(startObs:endYear, juvenile.ratio)

pdf("../Manuscript/figures/Immigration.pdf")
barplot(t(sum.mig[,2:3]), xlab = "Year", ylab = "Immigration")
legend(1, 6, c("Female", "Male"), fill = c("black", "gray"))
dev.off()

attach(VirusData)

infection<-subset(VirusData, select = SIV:GBVC)
#correlation between dichotomous variables
tetrachoric(infection)

set.seed(1)
glm.SIV<-glm(SIV ~ Sex + MinAge + Immigrant, family = binomial(link = "logit"), data = VirusData)
summary(glm.SIV)
cv.SIV<-cv.glm(VirusData, glm.SIV)
train<-sample(nrow(VirusData), (nrow(VirusData)*2)%/%3)
lda(SIV ~ Sex + MinAge + Immigrant, data = VirusData, subset = train)

glm.SHFV1<-glm(SHFV1 ~ .-SHFV1, family = binomial(link = "logit"), data = VirusData)
summary(glm.SHFV1)
step.SHFV1<-step(glm.SHFV1)
summary(glm(SHFV1~ Sex + MinAge + Immigrant, family = binomial(link = "logit"), data = VirusData))
# Sex is significant
glm.SHFV2<-glm(SHFV2 ~ .-SHFV2, family = binomial(link = "logit"), data = VirusData)
summary(glm.SHFV2)
step.SHFV2<-step(glm.SHFV2)
summary(glm(step.SHFV2$formula, family = binomial(link = "logit"), data = VirusData))
summary(glm(SHFV2 ~ Sex + MinAge + Immigrant, family = binomial(link = "logit"), data = VirusData)) 
# Sex is significant
summary(glm(GBVC ~ Sex + MinAge + Immigrant, family = binomial(link = "logit"), data = VirusData)) #Sex is marginally significant


TotalVirus<-SIV + SHFV1 + SHFV2 + GBVC
summary(lm(TotalVirus ~ Sex + MinAge + Immigrant)) #Sex is significant

sum.sample<-tapply(!is.na(SIV), list(Sex, MinAge), sum)
sum.SIV<-tapply(SIV, list(Sex, MinAge), sum)
percent.SIV<-sum.SIV/sum.sample
sum.SHFV1<-tapply(SHFV1, list(Sex, MinAge), sum)
percent.SHFV1<-sum.SHFV1/sum.sample
sum.SHFV2<-tapply(SHFV2, list(Sex, MinAge), sum)
percent.SHFV2<-sum.SHFV2/sum.sample
sum.GBVC<-tapply(GBVC, list(Sex, MinAge), sum)
percent.GBVC<-sum.GBVC/sum.sample

pdf("../Manuscript/figures/VirusPlot%01d.pdf", onefile=FALSE)
barplot(sum.sample, beside = TRUE, xlab = "Age", ylab = "Sample size", cex.lab = 1.5, cex.axis = 1.2, cex.names = 1.2)
legend(0, 6, c("Female", "Male"), fill = c("black", "gray"))
barplot(percent.SIV, beside = TRUE, xlab = "Age", ylab = "Proportion of infection", main = "SIV infection", cex.lab = 1.5, cex.axis = 1.2, cex.names = 1.2)
legend(0, 1, c("Female", "Male"), fill = c("black", "gray"))
barplot(percent.SHFV1, beside = TRUE, xlab = "Age", ylab = "Proportion of infection", main = "SHFV1 infection", cex.lab = 1.5, cex.axis = 1.2, cex.names = 1.2)
legend(0, 1, c("Female", "Male"), fill = c("black", "gray"))
barplot(percent.SHFV2, beside = TRUE, xlab = "Age", ylab = "Proportion of infection", main = "SHFV2 infection", cex.lab = 1.5, cex.axis = 1.2, cex.names = 1.2)
legend(0, 1, c("Female", "Male"), fill = c("black", "gray"))
barplot(percent.GBVC, beside = TRUE, xlab = "Age", ylab = "Proportion of infection", main = "GBVC infection", cex.lab = 1.5, cex.axis = 1.2, cex.names = 1.2)
legend(0, 1, c("Female", "Male"), fill = c("black", "gray"))
dev.off()

detach(VirusData)
