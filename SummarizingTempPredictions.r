##Calculating mean and range temperatures for the segment summary table (Table 1) in the Hydropsyche manuscript
library(matrixStats)

##About the TempPreds data
#This spreadsheet is monthly temp prediction output from Kim Dibble + Charles Yakulic temperature model for the upper basin (Dibble et al. 2018), based on data from USGS temp gages and other temp loggers from 1985-2015
	#I added a column called "SEGCODE" that matches the river segments with our study area.
		#colorado1 = C1
		#colorado1g = C2
		#colorado2 = Split in Grand Canyon  for segments "C3" and "C4" is at temp.station 103_colorado2 just upstream of Phantom Ranch. This division is not important for temperature, but was decided upon for the discharge data because there is a stream gage at Phantom.
		#colorado3 = below Hoover Dam, removed for this subset
		#green1 = G1
		#green2 = SPLIT INTO 3 SECTIONS
			#G2 = FGD to Yampa temp.station 1_green2 to  66_green2
			#G3 = Yampa to Jensen temp.station  67_green2 to 109_green2
			#G4 = Jensen to Green River, UT temp.station 110_green2 to 292_green2
			#G5 = Green River, UT to confluence with CO 293_green2 to 412_green2
		#san.juan = T4 (Navajo Dam to Lake Powell)
		#white = removed for this subset
		#yampa = T1 Deerlodge to confluence
		
#Temp predictions for the Gunnison and Dolores are from Gage data (see gage worksheets)


###Import CSV with Kim temp prediction + Segcodes
setwd("U:/!Hydropsyche2018")
TempPreds<-read.csv("KimTempPredictions_HydropsycheSEGCODES.csv")
head(TempPreds)

#Calculate mean temps + SD by segcode
MeanTemps<-as.data.frame(with(TempPreds, tapply(pred.temp.C,list(SEGCODE), FUN = function(x) mean(x))))
colnames(MeanTemps)<-"MeanTemp"
MeanTemps

SDTemps<-as.data.frame(with(TempPreds, tapply(pred.temp.C,list(SEGCODE), FUN = function(x) sd(x))))
colnames(SDTemps)<-"SDTemp"
SDTemps




#Calculate Range Temp
#Get monthly max and min per segcode 
	MaxTemps<-as.data.frame(with(TempPreds, tapply(pred.temp.C,list(SEGCODE, month), FUN = function(x) max(x))))
	MaxTemps$Max<-rowMaxs(as.matrix(MaxTemps)) 
	MaxTemps

	MinTemps<-as.data.frame(with(TempPreds, tapply(pred.temp.C,list(SEGCODE, month), FUN = function(x) min(x))))
	MinTemps$Min<-rowMins(as.matrix(MinTemps))
	MinTemps

#Get range for each site
RangeTempBySite<-as.data.frame(MaxTemps$Max-MinTemps$Min)
colnames(RangeTempBySite)<-"RangeTemp"
rownames(RangeTempBySite)<-rownames(MinTemps)
RangeTempBySite


##another column with SD of range temp
#dataframe with mean temperature range by month by seg
RangeTempMonthly<-as.data.frame(MaxTemps-MinTemps)
RangeTempMonthly<-RangeTempMonthly[,c(1:12)]
RangeTempMonthly


#standard deviation of max-min per month per site
RangeTempMonthly$RangeTempSD<-rowSds(as.matrix(RangeTempMonthly))
RangeTempMonthly


#Combine into a neat summary
TempSum<-cbind(MeanTemps,SDTemps, RangeTempBySite, RangeTempMonthly$RangeTempSD)
colnames(TempSum)<-c("MeanTemp","SDTemp","RangeTemp","SDRangeTemp")
TempSum

#Save CSV (also available on GitHub)
#write.csv(TempSum, "TempPredsSummary.csv")

