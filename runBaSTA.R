library(BaSTA)
library(snowfall)
#some extra code for more plots
#source("plotFancyBaSTA.R")
set.seed(101)

#the whole chuck of code is to reorganize the data for BaSTA
oriData<-read.csv("../Colin/forBaSTA201503.csv", head=T)
nInd<-nrow(oriData)
startYear<-min(oriData$Start)
startObs<-2006
endYear<-max(oriData$End)
nYear<-endYear-startYear+1
##some function
obsFill<-function(time1, time2, Im) {
  if(time1>=max(startObs, Im))
    c(rep(0,time1-startYear),rep(1, time2-time1+1),rep(0, endYear-time2))
  else
    c(rep(0,time1-startYear), 1, rep(0,max(startObs, Im)-time1-1), rep(1, time2-max(startObs, Im)+1),rep(0, endYear-time2))
}
##the full observation record
obsMat<- t(mapply(obsFill, oriData$Start, oriData$End, oriData$Migration))
##add birth and death years
birthDeath<-cbind(oriData$DOB,oriData$Death)
##permutate the sex of the infants/juveniles with unknown sex
#PermSexN<-sample(c(0,1), sum(oriData$SexN), replace=T)
#oriData$SexM[oriData$SexN==1]<-PermSexN
#oriData$SexF[oriData$SexN==1]<-1-PermSexN
sexMat<-cbind(oriData$SexM, oriData$SexF)
##combine all components into a matrix
inputMat<-cbind(1:nInd, birthDeath, obsMat, sexMat)
dimnames(inputMat)<-list(1:nInd, c("ID", "Birth", "Death", startYear:endYear, "SexM", "SexF"))
inputData<-as.data.frame(inputMat)
inputData<-subset(inputData, SexM + SexF == 1)
n.input<-nrow(inputData)
##observation must be 0 for birth year as requred by BaSTA
for(i in 1:n.input){
  if(inputData[i,]$Birth>0) inputData[i,][as.character(inputData[i,]$Birth)]<-0 
  if(inputData[i,]$Death>0) inputData[i,][as.character(inputData[i,]$Death)]<-0
}

#data check (if silent = F)
newData<-DataCheck(inputData, studyStart=2000, studyEnd=2015, autofix = rep(1, 7), silent=F)
#run basta to estimate the parameters for the chosen model (Here I chose Gombertz and bthtub shape)
out0<-basta(object = inputData, studyStart=2000, studyEnd=2015, model = "GO", shape = "simple", 
            recaptTrans=c(2000, 2006), niter=20000, burnin=2000, thinning=40, nsim = 4, parallel = TRUE, ncups = 4)
#view the summary of the result
summary(out0, digits=3)

#plot the bathtub mortality curve
plot(out0, plot.trace=F, noCI=T) #note the high error bars because we do not have age data of old individuals
