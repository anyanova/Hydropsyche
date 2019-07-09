###################################################
####################################################
###June 2019 HYDROPSYCHE STATS WORKSPACENUMBERS

#Hydropsyche data that was used for BRMS models
setwd("U:/!Hydropsyche2018/brms")
HY7<-read.csv("HY7.csv")
 
##How many samples? # 
dim(HY7) #2194

#How many Hydropsyche?
sum(HY7[,c(2:7)]) #16222

#What percent of samples contained Hydropsyche?
HYmath<-as.data.frame(rowSums(HY7[,c(2:7)]))
rownames(HYmath)<-HY7$BarcodeID
colnames(HYmath)<-"TotalHY"
dim(HYmath) #still 2194 total
head(HYmath)

#How many zeroes?
length(HYmath[HYmath$TotalHY==0,]) #1738 zeroes
length(HYmath[HYmath$TotalHY!=0,]) #456 non zeroes


##Species specific abundance
sum(HY7$HYOC) #9773, 60.25% 
sum(HY7$HYOS) #5147, 31.73%
sum(HY7$HYCO) #148, 0.91%
sum(HY7$HYCA) #1018, 6.28%
sum(HY7$HYQU) #126, 0.78%
sum(HY7$HCOC) #10, 0.06%




##Table of catch rate by SEGCODE, one sp at a time
HYCA<-as.data.frame(with(HY7, tapply(HYCA,list(SEGCODE), FUN = function(x) mean(na.omit(x)))))
colnames(HYCA)<-"HYCA"

HCOC<-as.data.frame(with(HY7, tapply(HCOC,list(SEGCODE), FUN = function(x) mean(na.omit(x)))))
colnames(HCOC)<-"HCOC"

HYCO<-as.data.frame(with(HY7, tapply(HYCO,list(SEGCODE), FUN = function(x) mean(na.omit(x)))))
colnames(HYCO)<-"HYCO"

HYOC<-as.data.frame(with(HY7, tapply(HYOC,list(SEGCODE), FUN = function(x) mean(na.omit(x)))))
colnames(HYOC)<-"HYOC"

HYOS<-as.data.frame(with(HY7, tapply(HYOS,list(SEGCODE), FUN = function(x) mean(na.omit(x)))))
colnames(HYOS)<-"HYOS"

HYQU<-as.data.frame(with(HY7, tapply(HYQU,list(SEGCODE), FUN = function(x) mean(na.omit(x)))))
colnames(HYQU)<-"HYQU"

HYcatch<-cbind(HYCA,HCOC,HYCO,HYOC,HYOS,HYQU)
HYcatch$SegCode<-rownames(HYcatch)
HYcatch

##stacked barplot for Fig1

#melt dataframe
library(reshape)
mHYcatch<-melt(HYcatch, id=c("SegCode"))
colnames(mHYcatch)<-c("SegCode","Spp","CatchRate")
head(mHYcatch)


mHYcatch$SegCode<- factor(mHYcatch$SegCode,levels = c("G1","G2","T1","G3","G4","G5","T2","C1","T3","C2","T4","C3","C4"))


##STACKED BAR PLOT OF CATCH RATES BY SPECIES##
HYbar<-ggplot(mHYcatch,aes(x = SegCode, y=CatchRate, fill=as.factor(mHYcatch$Spp))) + xlab("Segment ID") + ylab("Catch rate") + 
	geom_bar(stat="identity", position="stack") + 
	scale_fill_manual(as.factor(mHYcatch$Spp), values=c("red","orange","yellow","navyblue","springgreen1","magenta"))	+
	theme_bw(base_size=20)
HYbar


				     
