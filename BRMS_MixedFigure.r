##Scatter plots + BRMS model output for figure 3 in Hydropsyche Dist Paper
library(ggplot2)
library(brms)
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

setwd("U:/!Hydropsyche2018/FIGURES")

HY2<-read.csv("U:/!Hydropsyche2018/brms/HY7.csv")
HY3<-read.csv("U:/!Hydropsyche2018/brms/HY7.csv")

##Set zeroes below 0
HY3$SUM<-rowSums(HY3[,c(2:7)])
HY3$black0<-ifelse(HY3$SUM<1, -3,-100)
HY3$HYOC<-ifelse(HY3$HYOC<1,-100, HY3$HYOC)
HY3$HYOS<-ifelse(HY3$HYOS<1,-100, HY3$HYOS)
HY3$HYCO<-ifelse(HY3$HYCO<1,-100, HY3$HYCO)
HY3$HYCA<-ifelse(HY3$HYCA<1,-100, HY3$HYCA)
HY3$HYQU<-ifelse(HY3$HYQU<1,-100, HY3$HYQU)
HY3$HCOC<-ifelse(HY3$HCOC<1,-100, HY3$HCOC)

##Stage and Q
HY3$MeanQ_m3s<-HY2$MeanQ * 0.028316846592 #change to m3s

head(HY3)

############################################################

#Scatterplots, small for paper
ggStage<- ggplot(aes(Stage_cm, HYOC), data=HY3) +
ylab("Hydropsyche catch rate (#/hr)") + ylim(-12,600) +
geom_point(aes(Stage_cm,black0),  col="grey", alpha=0.6, cex=.7)+
geom_point(aes(Stage_cm,HYOS),  col="springgreen1", alpha=0.6, cex=.7)+
geom_point(aes(Stage_cm,HYOC),  col="navyblue", alpha=0.6, cex=.7)+
geom_point(aes(Stage_cm,HYCO),  col="yellow", alpha=0.6, cex=.7)+
geom_point(aes(Stage_cm,HYQU),  col="magenta", alpha=0.6, cex=.7)+
geom_point(aes(Stage_cm,HYCA),  col="red", alpha=0.6, cex=.7)+
geom_point(aes(Stage_cm,HCOC),  col="orange", alpha=0.6, cex=.7)+
theme_bw() + scale_x_continuous("Stage change (cm)", seq(0,120,20))
# + geom_vline(aes(xintercept=13), linetype="dashed", color="red")
ggStage 


ggmt<- ggplot(aes(lat, HYOC), data=HY3) + xlab("Mean Temperature (°C)") + 
ylab("Hydropsyche catch rate (#/hr)") + ylim(-12,600) +
geom_point(aes(lat,black0),  col="grey", alpha=0.6, cex=.7)+
geom_point(aes(lat,HYOS),  col="springgreen1", alpha=0.6, cex=.7)+
geom_point(aes(lat,HYOC),  col="navyblue", alpha=0.6, cex=.7)+
geom_point(aes(lat,HYCO),  col="yellow", alpha=0.6, cex=.7)+
geom_point(aes(lat,HYQU),  col="magenta", alpha=0.6, cex=.7)+
geom_point(aes(lat,HYCA),  col="red", alpha=0.6, cex=.7)+
geom_point(aes(lat,HCOC),  col="orange", alpha=0.6, cex=.7)+
theme_bw()
ggmt

ggrt<- ggplot(aes(RangeTemp, HYOC), data=HY3) + xlab("Temperature Range (°C)") + 
ylab("Hydropsyche catch rate (#/hr)") + ylim(-12,600) +
geom_point(aes(RangeTemp,black0),  col="grey", alpha=0.6, cex=.7)+
geom_point(aes(RangeTemp,HYOS),  col="springgreen1", alpha=0.6, cex=.7)+
geom_point(aes(RangeTemp,HYOC),  col="navyblue", alpha=0.6, cex=.7)+
geom_point(aes(RangeTemp,HYCO),  col="yellow", alpha=0.6, cex=.7)+
geom_point(aes(RangeTemp,HYQU),  col="magenta", alpha=0.6, cex=.7)+
geom_point(aes(RangeTemp,HYCA),  col="red", alpha=0.6, cex=.7)+
geom_point(aes(RangeTemp,HCOC),  col="orange", alpha=0.6, cex=.7)+
theme_bw()
ggrt

multiplot(ggmt,ggrt,ggStage,cols=1)



###PLOTTING BRMS OUTPUT MANUALLY!!!!\
#http://www.flutterbys.com.au/stats/tut/tut7.2b.html

##Load all necessary model outputs into one R window
load("U:/!Hydropsyche2018/brms/Rworkspaces/BestMods2019/brms_HYOC27.RData")
load("U:/!Hydropsyche2018/brms/Rworkspaces/BestMods2019/brms_HYOS27.RData")


HYOC<-HYOC27
HYOS<-HYOS27

#Export desired variable data out of marginal_effects into a usable format for GGplot, and plot multiple species in one plot
summary(marginal_effects(HYOC))
	HYOCmt_man<-marginal_effects(HYOC)$lat
	head(HYOCmt_man)
	HYOClat_man<-marginal_effects(HYOC)$lat
	head(HYOClat_man)
summary(marginal_effects(HYOS))
	HYOSrt_man<-marginal_effects(HYOS)$RangeTemp
	head(HYOSrt_man)
summary(marginal_effects(HYOS))
	HYOSsc_man<-marginal_effects(HYOS)$Stage_cm
	head(HYOSsc_man)
	
summary(marginal_effects(HYOC))
	HYOCDay_man<-marginal_effects(HYOC)$Jul
	head(HYOCDay_man)
summary(marginal_effects(HYOS))
	HYOSDay_man<-marginal_effects(HYOS)$Jul
	head(HYOSDay_man)


#Manually Relabel the X axes with calculations so it is not Z scored (Z = (x-mean(x))/sd(x)
		lat<-mean(HY2$lat)
		sdTEMP<-sd(HY2$lat)
		#ie: What is the Z score of actual value 0
			(0-lat)/sdTEMP 
			
		RTEMP<-mean(HY2$RangeTemp)
		sdRTEMP<-sd(HY2$RangeTemp)
		#ie: What is the Z score of actual value 10
			(2-RTEMP)/sdRTEMP #1.104972

			
		meanSC<-mean(HY2$Stage_cm)
		sdSC<-sd(HY2$Stage_cm)
		#ie: What is the Z score of actual value 2 cm of stage change?
			(80-meanSC)/sdSC 
			
		LAT<-mean(HY2$lat)
		sdLAT<-sd(HY2$lat)
		#ie: What is the Z score of actual value 10
			(37-LAT)/sdLAT #1.104972
			
		DAY<-mean(HY2$Jul)
		sdDAY<-sd(HY2$Jul)
		#ie: What is the Z score of actual value 10
			(1-DAY)/sdDAY #1.104972
			

#Plot BRMS output 


#Day
DAY<-ggplot() + 
	geom_line(data=HYOCDay_man, aes(y = estimate__, x = Jul)) + 
	geom_ribbon(data=HYOCDay_man, aes(y = estimate__, x = Jul, ymin = lower__,ymax = upper__), fill = "navyblue", alpha = 0.6) + 
	scale_y_continuous("Modeled abundance") + scale_x_continuous("Month", breaks=c(-2.686736, -1.803547, -0.8904204,  0.02270638,0.9508025,  1.863929), labels=c("Jan", "Mar", "May", "July", "Sept", "Nov")) + 
	theme_bw()
DAY


#Day
DAY<-ggplot() + 
	geom_line(data=HYOSDay_man, aes(y = estimate__, x = Jul)) + 
	geom_ribbon(data=HYOSDay_man, aes(y = estimate__, x = Jul, ymin = lower__,ymax = upper__), fill = "springgreen1", alpha = 0.6) + 
	scale_y_continuous("Modeled abundance") + scale_x_continuous("Month", breaks=c(-2.686736, -1.803547, -0.8904204,  0.02270638,0.9508025,  1.863929), labels=c("Jan", "Mar", "May", "July", "Sept", "Nov")) + 
	theme_bw()
DAY



#Temp
Temp<-ggplot() + 
	geom_line(data=HYOClat_man, aes(y = estimate__, x = lat)) + 
	geom_ribbon(data=HYOClat_man, aes(y = estimate__, x = lat, ymin = lower__,ymax = upper__), fill = "navyblue", alpha = 0.6) + 
	scale_y_continuous("Modeled abundance") + scale_x_continuous("Latitude", breaks=c(-3.265594, -1.544379, 0.1768356, 1.89805), labels=c("8", "10", "12","14")) + 
	theme_bw()
Temp


#lat
Lat<-ggplot() + 
	geom_line(data=HYOClat_man, aes(y = estimate__, x = lat)) + 
	geom_ribbon(data=HYOClat_man, aes(y = estimate__, x = lat, ymin = lower__,ymax = upper__), fill = "navyblue", alpha = 0.6) + 
	scale_y_continuous("Modeled abundance") + scale_x_continuous("Latitude", breaks=c(-0.6663728, 0.05355411, 0.773481, 1.493408, 2.213335, 2.933262,  3.653189,  4.373115), labels=c("36", "37", "38", "39", "40", "41", "42", "43")) + 
	theme_bw()
Lat



##Range Temp
RTemp<-ggplot() + 
	geom_line(data=HYOSrt_man, aes(y = estimate__, x = RangeTemp)) + 
	geom_ribbon(data=HYOSrt_man, aes(y = estimate__, x = RangeTemp, ymin = lower__, ymax = upper__), fill = "springgreen1", alpha = 0.6) + 
	scale_y_continuous("Modeled abundance") + scale_x_continuous("Range temperature (°C)", breaks=c(-0.5397707, 0.1790343,  0.8978394, 1.616644, 2.335449), labels=c("5", "10", "15", "20", "25")) +
	theme_bw()
RTemp

##Stage
STAGE<-ggplot() + 
	geom_line(data=HYOSsc_man, aes(y = estimate__, x = Stage_cm)) + 
	geom_ribbon(data=HYOSsc_man, aes(y = estimate__, x = Stage_cm, ymin = lower__,ymax = upper__), fill = "springgreen1", alpha = 0.6) + 
	scale_y_continuous("Modeled abundance") + scale_x_continuous("Stage Change (cm)", breaks=c(-1.898737, -1.196596, -0.4944544,  0.207687, 0.9098283, 1.61197), labels=c("0", "20", "40", "60", "80", "100")) +
	theme_bw()
STAGE



	
#Final figure! Put it all together
multiplot(ggmt,ggrt,ggStage, Temp, RTemp, STAGE, cols=2)




#Species richness vs stage change
#Create P/A data species specific and species richness column SPP
HYPA<-HY7
HYPA$HYOC<-ifelse(HYPA$HYOC<1,0, 1)
HYPA$HYOS<-ifelse(HYPA$HYOS<1,0, 1)
HYPA$HYCO<-ifelse(HYPA$HYCO<1,0, 1)
HYPA$HYCA<-ifelse(HYPA$HYCA<1,0, 1)
HYPA$HYQU<-ifelse(HYPA$HYQU<1,0, 1)
HYPA$HCOC<-ifelse(HYPA$HCOC<1,0, 1)
HYPA$SPP<-rowSums(HYPA[,c(2:7)])

#Add total counts sum column
HYcount<-HY3[,c(1,21)]
HYPA<-merge(HYcount, HYPA, by="BarcodeID")
head(HYPA)
table(HYPA$SUM)

ggpa<-ggplot(aes(Stage_cm, SPP), data=HYPA) + xlab("Stage change (cm)") + 
ylab("Hydropsyche species richness") + geom_point() +
geom_smooth(method='lm') + theme_bw() + ylim(0,3)
ggpa 


lm(SUM~Stage_cm, data=HYPA)
abline(60.44, -23.11)






#Plot mean species richness vs. mean stage change by river segment
head(HYPA)
HYspp<-with(HYPA, tapply(SPP,list(SEGCODE), FUN = function(x) mean (unique(na.omit(x)))))
HYsr<-as.data.frame(HYspp)
HYstage<-with(HYPA, tapply(Stage_cm,list(SEGCODE), FUN = function(x) mean (unique(na.omit(x)))))
HYsc<-as.data.frame(round(HYstage,2))
HYplot<-cbind(HYsc,HYsr)
colnames(HYplot)<-c("Stage_cm", "HYsr")
HYplot$SegCode<-rownames(HYplot)


library(ggrepel)

##PLOT of MEAN SPECIES RICHNESS / SEGCODE
ggsr<-ggplot(aes(Stage_cm, HYsr, label=SegCode), data=HYplot) + xlab("Stage change (cm)") + 
ylab("Hydropsyche mean species richness") +
geom_label_repel(aes(label = SegCode), box.padding   = 0.35, point.padding = 0, segment.color = NA)+
geom_smooth(method='lm') + theme_bw()
ggsr



##Plot TOTAL species richness per segment with mean Stage change from table 1 in manuscript
HYtab<-data.frame(rownames(HYsr))
HYtab$TotalSR<-c(3,2,1,1,3,3,2,2,2,2,1,2,3)
HYtab$MeanStage<-c(6.7,5.6,70.61,69.36,2.8,29.5,12.1,10.1,4.8,6.1,3,6.1,9.7)
colnames(HYtab)<-c("SegCode","TotalSR","MeanStage")
HYtab


ggsr<-ggplot(aes(MeanStage, TotalSR, label=SegCode), data=HYtab) + xlab("Stage change (cm)") + 
scale_y_continuous("Hydropsyche species richness", seq(1,3,1)) + geom_point() +
geom_label_repel(aes(label = SegCode), box.padding   = 0.35, point.padding = 0.5, segment.color = "grey15")+
geom_smooth(method='lm', lty=2, se=FALSE, col="grey" ) + theme_bw()
ggsr 









####CALCULATE NEW STAGE CHANGE FOR TABLE 1 IN GRAND CANYON WITH JEFE'S DATA
setwd("U:/!Hydropsyche2018/brms/Rworkspaces/NewStage")
s6 <- read.csv("StageChangePerRM.csv")
s6$SEGCODE<-ifelse(s6$RM<88.6, "C3", "C4")
head(s6)

with(s6, tapply(MeanDeltaStage,list(SEGCODE), FUN = function(x) mean(x)))















