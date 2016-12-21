#the whole chuck of code is to reorganize the data of RC monkeys from the 2015 March record
library(dplyr)
library(stringr)

#read in the two files
data.demo<-read.csv("../Data/RCDemo201503.csv", header = T, na.strings = c(""))
data.infection<-read.csv("../Data/RCinfection.csv", header = T, na.strings = c(""))
data.basic<-select(data.demo, Name, VirusID, Sex, DOB, Death)

#visualize the pattern of missing data
library(mice)
md.pattern(data.basic)
library(VIM)
aggr(data.basic, prop=TRUE, numbers=TRUE)
matrixplot(data.basic)

#split the names into three columns (generations)
pedigree<-str_split_fixed(data.demo$Name, "-", 3)
data.demo$G1<-pedigree[,1]
data.demo$G2<-pedigree[,2]
data.demo$G3<-pedigree[,3]
#count the number of offspring of each female
data.demo$n.offspring<-data.demo$n.g2<-data.demo$n.g3<-0
n.g2<-data.demo %>% filter(G2!="") %>% group_by(G1) %>% summarise(n.g2=n_distinct(G2))
n.g3<-data.demo %>% filter(G3!="") %>% group_by(G2) %>% summarise(n.g3=n_distinct(G3))
list.g1<-data.frame(G1=data.demo$G1[data.demo$G2==""]) 
list.g1<-left_join(list.g1, n.g2, by="G1")
list.g1$n.g2[is.na(list.g1$n.g2)]<-0
list.g2<-data.frame(G2=data.demo$G2[data.demo$G3==""&data.demo$G2!=""]) 
list.g2<-left_join(list.g2, n.g3, by="G2")
list.g2$n.g3[is.na(list.g2$n.g3)]<-0
data.demo$n.g2[data.demo$G2==""]<-list.g1$n.g2
data.demo$n.g3[data.demo$G3==""&data.demo$G2!=""]<-list.g2$n.g3
data.demo$n.offspring<-rowSums(select(data.demo, n.g2, n.g3), na.rm=TRUE)
hist(data.demo$n.offspring[data.demo$Sex=="F"], main = "Offspring of females", 
     xlab = "Number of offspring")

#convert date info to year
data.infection$Collection.Year<-as.integer(format(as.Date(data.infection$Collection.date, "%m/%d/%y"), "%Y"))
data.infection<-select(data.infection, Animal.ID, Collection.Year, SIV=SIVkrc, SHFV1=SHFV.krc1, SHFV2=SHFV.krc2, GBVC=GBV.Ckrc)
#join the two data frames
Data0<-left_join(data.demo, data.infection, by=c("VirusID" = "Animal.ID"))

nInd<-nrow(Data0)
startObs<-2006
endYear<-max(Data0$End)
Emig<-Data0$RightCensor
Emig[Emig==2015]<-'NA'
Data0$Emig<-as.factor(Emig)

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
VirusData<-cbind(Data0$Sex, ImmigValue, EmigValue, age.collection, age.collection.min, Data0$n.offspring, Data0$SIV, Data0$SHFV1, Data0$SHFV2, Data0$GBVC)

InGroupData<-as.data.frame(InGroupData)
AgeData<-as.data.frame(AgeData)
VirusData<-as.data.frame(VirusData)

#rename columns
names(InGroupData)<-c("Sex", "Immig", "Emig", "Death", startObs:endYear)
names(AgeData)<-c("Sex", startObs:endYear)
names(VirusData)<-c("Sex", "Immigrant", "Emigrant", "Age", "MinAge", "Offspring", "SIV", "SHFV1", "SHFV2", "GBVC")

#delete the 2015 data because it is incomplete
InGroupData<-InGroupData[,-ncol(InGroupData)]
AgeData<-AgeData[,-ncol(AgeData)]

InGroupData$Sex<-factor(InGroupData$Sex, levels = 1:2, labels = c("F", "M"))
AgeData$Sex<-factor(AgeData$Sex, levels = 1:2, labels = c("F", "M"))
VirusData$Sex<-factor(VirusData$Sex, levels = 1:2, labels = c("F", "M"))

write.csv(InGroupData, "InGroupData.csv", row.names=FALSE)
write.csv(AgeData, "AgeData.csv", row.names = FALSE)
write.csv(VirusData, "VirusData.csv", row.names=FALSE)
