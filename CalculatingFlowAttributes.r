#Calculating flow attributes for Lees Ferry
library(lubridate)

#import csv
	LEESgage<-read.csv('C:/Users/ametcalfe/Desktop/GageData/Gages_QandHeight/CR4_LEES.csv')
	
#Limit to years 2012-2016 (Remove hash to open function, which was breaking everything)
	LEESgage$Date<-as.Date(LEESgage$DateTime, format="%m/%d/%Y")
	LEESgage$Year<-year(LEESgage$Date)
	LEESgage<-LEESgage[LEESgage$Year<2017,]
		
	
head(LEESgage)
LEESgage_mdd<-tapply(LEESgage$Discharge,LEESgage$Date, mean)
LEESgage_sdd<-tapply(LEESgage$Discharge,LEESgage$Date, sd)
LEESgage$MeanQ<-mean(LEESgage_mdd, na.rm=TRUE)
LEESgage$SDQ<-mean(LEESgage_sdd, na.rm=TRUE)
LEESgage$HI<-mean(LEESgage_sdd/LEESgage_mdd, na.rm=TRUE)
#Mean daily stage change (average of daily max-min over 5 years)
LEESgage_hmax<-tapply(LEESgage$Height,LEESgage$Date,max)
LEESgage_hmin<-tapply(LEESgage$Height,LEESgage$Date,min)
LEESgage$StageChange<-mean(LEESgage_hmax - LEESgage_hmin, na.rm=TRUE)
LEES<-LEESgage[1,c(1,7:10),] #One line data frame for LEES with Site, MeanQ,SDQ, and HI
LEES
