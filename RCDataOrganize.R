#the whole chuck of code is to reorganize the data of RC monkeys from the 2015 March record
library(dplyr)
data.demo<-read.csv("../Colin/RCDemo201503.csv", header = T, na.strings = c(""))
data.infection<-read.csv("../Tony/RCinfection.csv", header = T, na.strings = c(""))
#convert date info to year
data.infection$Collection.Year<-as.integer(format(as.Date(data.infection$Collection.date, "%m/%d/%y"), "%Y"))
#select function is masked by MASS!
data.infection<-dplyr::select(data.infection, Animal.ID, Collection.Year, SIV=SIVkrc, SHFV1=SHFV.krc1, SHFV2=SHFV.krc2, GBVC=GBV.Ckrc)
#join the two data frames
Data0<-left_join(data.demo, data.infection, by=c("VirusID" = "Animal.ID"))

nInd<-nrow(Data0)
startObs<-2006
endYear<-max(Data0$End)
Emig<-Data0$RightCensor
Emig[Emig==2015]<-'NA'
Data0$Emig<-Emig

#Migration by sex and year
ImmigValue<-!is.na(Data0$Migration)
EmigValue<-!is.na(Data0$Emig)

##Functions to transform the life table
GroupStatusFill<-function(In, End) {
  c(rep(0,In-startObs),rep(1, End-In+1),rep(0, endYear-End))
}

MinAgeFill<-function(Start, End){
  c(rep(0, max(Start-startObs,0)), (max(startObs-Start,0)+1):(End-Start+1), rep(0, endYear-End))
}

##the in-group record
InGroupStart<- pmax(Data0$Migration, Data0$Start, startObs, na.rm = TRUE)
InGroupMat<- t(mapply(GroupStatusFill, InGroupStart, Data0$End))

##the age at blood sample collection
age.collection.min<-Data0$Collection.Year-Data0$Start
age.collection<-Data0$Collection.Year-Data0$DOB

##the age record
AgeMat<- t(mapply(MinAgeFill, Data0$Start, Data0$End))
##combine all components into a matrix
InGroupData<-cbind(Data0$Sex, Data0$Migration, Data0$Emig, Data0$Death, InGroupMat)
AgeData<-cbind(Data0$Sex, AgeMat*InGroupMat)
VirusData<-cbind(Data0$Sex, ImmigValue, EmigValue, age.collection, age.collection.min, Data0$SIV, Data0$SHFV1, Data0$SHFV2, Data0$GBVC)

InGroupData<-as.data.frame(InGroupData)
AgeData<-as.data.frame(AgeData)
VirusData<-as.data.frame(VirusData)

#rename columns
names(InGroupData)<-c("Sex", "Immig", "Emig", "Death", startObs:endYear)
names(AgeData)<-c("Sex", startObs:endYear)
names(VirusData)<-c("Sex", "Immigrant", "Emigrant", "Age", "MinAge", "SIV", "SHFV1", "SHFV2", "GBVC")

InGroupData<-InGroupData[,-ncol(InGroupData)]
AgeData<-AgeData[,-ncol(AgeData)]

InGroupData$Sex<-factor(InGroupData$Sex, levels = 1:2, labels = c("F", "M"))
AgeData$Sex<-factor(AgeData$Sex, levels = 1:2, labels = c("F", "M"))
VirusData$Sex<-factor(VirusData$Sex, levels = 1:2, labels = c("F", "M"))

write.csv(InGroupData, "InGroupData.csv", row.names=FALSE)
write.csv(AgeData, "AgeData.csv", row.names = FALSE)
write.csv(VirusData, "VirusData.csv", row.names=FALSE)

