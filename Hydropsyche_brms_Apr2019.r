###APRIL 2019 RERUN MODELS WITH CORRECTED LOESS LINE FOR STAGE CHANGE IN GRAND CANYON
##brms zero-inflated (NEEDS RTOOLS C++ capability)

setwd("U:/!Hydropsyche2018/brms/Rworkspaces/NewStage")
library(brms)
library(rstan)
#speed up processing
	rstan_options (auto_write=TRUE)
	options (mc.cores=parallel::detectCores ())

#load dataset
HY7<-read.csv("HY7.csv")

HY<-HY7
##Split HY into Grand Canyon and non Grand Canyon
	HYgc<-HY[HY$Reach=="CRGrandCanyon",]
	##remove old stage change measurements
	HYgc<-HYgc[,c(1:18)]

	HYub<-HY[HY$Reach!="CRGrandCanyon",]
	dim(HYgc);dim(HYub)

#Import Jeff's Stage Change calculations from Grand Canyon Q model
	s6 <- read.csv("GC_StageChangePerRM.csv")

##Use predict to reduce stochasticity in predictions because the old predictions were all over the place
s5 <- s6[order(s6$RM),]
rng <- seq(round(min(s5$RM, na.rm = TRUE), 2), round(max(s5$RM, na.rm = TRUE), 2), 0.01)
pred <- predict(loess(s5$MeanDeltaStage ~ s5$RM, span = 0.3), rng)
s7 <- data.frame(RM = rng, MeanDeltaStage = pred)
HYgc$Stage_cm <- s6[which.min(abs(s6$RM - HYgc$RiverMile[i])), "MeanStage"]

#Use a for loop to associate modeled stage change with light trap samples using river mile
	HYgc$Stage_cm<- NA
	for(i in 1:dim(HYgc)[1]){
	HYgc$Stage_cm[i] <- s7[which.min(abs(s7$RM - HYgc$RiverMile[i])), "MeanDeltaStage"]
	}

##Recombine GC & UB samples into one dataframe for models
	HY<-rbind(HYgc,HYub)
	dim(HY)
head(HY)
	
#write.csv(HY, "HY8.csv")	##NON STANDARDIZED!!
	
#Standardize model vars with Z scores AFTER UPDATING STAGE!!
	#Z score ENV vars
		HY$Jul<-scale(as.numeric(HY$Jul))
		HY$MeanQ<-scale(HY$MeanQ)
		HY$Stage_cm<-scale(HY$Stage_cm)
		HY$lat<-scale(as.numeric(HY$lat))
		HY$MeanTemp<-scale(HY$MeanTemp)
		HY$RangeTemp<-scale(HY$RangeTemp)
	head(HY)

	
#RUN HYOC STAGE MODELS
	#Stage Change
		HYOC19<-brm(bf(HYOC~Stage_cm, zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC20<-brm(bf(HYOC~1, zi~Stage_cm + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC21<-brm(bf(HYOC~Stage_cm, zi~Stage_cm + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC19); waic(HYOC19)
		summary(HYOC20); waic(HYOC20)
		summary(HYOC21); waic(HYOC21)
		#StageChangeQuad 
		HYOC19q<-brm(bf(HYOC~Stage_cm + I(Stage_cm^2), zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC20q<-brm(bf(HYOC~1, zi~Stage_cm + I(Stage_cm^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC21q<-brm(bf(HYOC~Stage_cm + I(Stage_cm^2), zi~Stage_cm + I(Stage_cm^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC19q); waic(HYOC19q)
		summary(HYOC20q); waic(HYOC20q)
		summary(HYOC21q); waic(HYOC21q)
		#StageChangeMixed
		HYOC19m<-brm(bf(HYOC~Stage_cm + I(Stage_cm^2), zi~Stage_cm + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC20m<-brm(bf(HYOC~Stage_cm, zi~Stage_cm + I(Stage_cm^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC19m); waic(HYOC19m)
		summary(HYOC20m); waic(HYOC20m)
		
###Save workspace for models 16 to 21
save.image("U:/!Hydropsyche2018/brms/Rworkspaces/NewStage/brms_HYOC_newSTAGE_19to21.RData")


##HYOC MULTIVAR	
		HYOC26<-brm(bf(HYOC~MeanTemp + I(MeanTemp^2) + Stage_cm + I(Stage_cm^2), zi~MeanTemp + I(MeanTemp^2) + Stage_cm + I(Stage_cm^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC26); waic(HYOC26)
		HYOC30<-brm(bf(HYOC~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2) + Stage_cm + I(Stage_cm^2), zi~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2) + Stage_cm + I(Stage_cm^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC30); waic(HYOC30)
				
	###Save workspace
	save.image("U:/!Hydropsyche2018/brms/Rworkspaces/NewStage/HYOC/brms_HYOC_newSTAGE_25to30.RData")		

		HYOC31<-brm(bf(HYOC~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2) + lat + Stage_cm + I(Stage_cm^2) + lat + RangeTemp + I(RangeTemp^2), zi~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2) + lat + RangeTemp + I(RangeTemp^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		HYOC32<-brm(bf(HYOC~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2) + lat +  MeanQ + I(MeanQ^2), zi~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2) + lat +  MeanQ + I(MeanQ^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		HYOC33<-brm(bf(HYOC~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2) + lat + Stage_cm + I(Stage_cm^2), zi~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2) + lat + Stage_cm + I(Stage_cm^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		summary(HYOC31); waic(HYOC31)
		summary(HYOC32); waic(HYOC32)
		summary(HYOC33); waic(HYOC33)
		###Save workspace
	save.image("U:/!Hydropsyche2018/brms/Rworkspaces/NewStage/HYOC/brms_HYOC_newSTAGE_31to33.RData")	


#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#

##HYOS

	#Stage Change
		HYOS19<-brm(bf(HYOS~Stage_cm, zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOS20<-brm(bf(HYOS~1, zi~Stage_cm + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOS21<-brm(bf(HYOS~Stage_cm, zi~Stage_cm + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		summary(HYOS19); waic(HYOS19)
		summary(HYOS20); waic(HYOS20)
		summary(HYOS21); waic(HYOS21)
		#StageChangeQuad 
		HYOS19q<-brm(bf(HYOS~Stage_cm + I(Stage_cm^2), zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOS20q<-brm(bf(HYOS~1, zi~Stage_cm + I(Stage_cm^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOS21q<-brm(bf(HYOS~Stage_cm + I(Stage_cm^2), zi~Stage_cm + I(Stage_cm^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOS19q); waic(HYOS19q)
		summary(HYOS20q); waic(HYOS20q)
		summary(HYOS21q); waic(HYOS21q)
		#StageChangeMixed
		HYOS19m<-brm(bf(HYOS~Stage_cm + I(Stage_cm^2), zi~Stage_cm + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOS20m<-brm(bf(HYOS~Stage_cm, zi~Stage_cm + I(Stage_cm^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		summary(HYOS19m); waic(HYOS19m)
		summary(HYOS20m); waic(HYOS20m)
		
	###Save workspace for models 19 to 21
	save.image("U:/!Hydropsyche2018/brms/Rworkspaces/NewStage/brms_HYOS_newSTAGE_19to21.RData")
	
	
	##HYOS MULTIVAR
		HYOS22<-brm(bf(HYOS~Stage_cm + RangeTemp + I(RangeTemp^2), zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))	
		HYOS23<-brm(bf(HYOS~Stage_cm + lat + I(lat^2), zi~lat + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOS24<-brm(bf(HYOS~Stage_cm + MeanQ, zi~MeanQ + I(MeanQ^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOS25<-brm(bf(HYOS~Stage_cm + MeanTemp + I(MeanTemp^2), zi~MeanTemp + I(MeanTemp^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOS26<-brm(bf(HYOS~Stage_cm + Jul + I(Jul^2), zi~Jul + I(Jul^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOS22); waic(HYOS22)
		summary(HYOS23); waic(HYOS23)
		summary(HYOS24); waic(HYOS24)
		summary(HYOS25); waic(HYOS25)
		summary(HYOS26); waic(HYOS26)
		summary(HYOS27); waic(HYOS27)
	###Save workspace for models 12 to 27
	save.image("U:/!Hydropsyche2018/brms/Rworkspaces/NewStage/brms_HYOS_newSTAGE_22to26.RData")
	
	HYOS27<-brm(bf(HYOS~Stage_cm + Jul + I(Jul^2) + RangeTemp + I(RangeTemp^2), zi~Jul + I(Jul^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
	HYOS28<-brm(bf(HYOS~Stage_cm + Jul + I(Jul^2) + lat + I(lat^2), zi~Jul + I(Jul^2) + lat + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
	HYOS29<-brm(bf(HYOS~Stage_cm + Jul + I(Jul^2) + MeanQ, zi~Jul + I(Jul^2) + MeanQ + I(MeanQ^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
	HYOS30<-brm(bf(HYOS~Stage_cm + Jul + I(Jul^2) + MeanTemp + I(MeanTemp^2), zi~Jul + I(Jul^2) + MeanTemp + I(MeanTemp^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		summary(HYOS27); waic(HYOS27)
		summary(HYOS28); waic(HYOS28)
		summary(HYOS29); waic(HYOS29)
		summary(HYOS30); waic(HYOS30)
	###Save workspace for models 19 to 21
	save.image("U:/!Hydropsyche2018/brms/Rworkspaces/NewStage/brms_HYOS_newSTAGE_27to30.RData")
	
	HYOS31<-brm(bf(HYOS~Stage_cm + Jul + I(Jul^2) + RangeTemp + I(RangeTemp^2) + lat + I(lat^2), zi~Jul + I(Jul^2) + lat + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
	HYOS32<-brm(bf(HYOS~Stage_cm + Jul + I(Jul^2) + RangeTemp + I(RangeTemp^2) + MeanQ, zi~Jul + I(Jul^2) + MeanQ + I(MeanQ^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
	HYOS33<-brm(bf(HYOS~Stage_cm + Jul + I(Jul^2) + RangeTemp + I(RangeTemp^2) + MeanTemp + I(MeanTemp^2), zi~Jul + I(Jul^2) + MeanTemp + I(MeanTemp^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		summary(HYOS31); waic(HYOS31)
		summary(HYOS32); waic(HYOS32)
		summary(HYOS33); waic(HYOS33)
	###Save workspace for models 31 to 33
	save.image("U:/!Hydropsyche2018/brms/Rworkspaces/NewStage/brms_HYOS_newSTAGE_31to33.RData")
	
		
		
		
		
		
#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#

##HYCO



	#Stage Change
		HYCO19<-brm(bf(HYCO~Stage_cm, zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYCO20<-brm(bf(HYCO~1, zi~Stage_cm + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYCO21<-brm(bf(HYCO~Stage_cm, zi~Stage_cm + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYCO19); waic(HYCO19)
		summary(HYCO20); waic(HYCO20)
		summary(HYCO21); waic(HYCO21)
		#StageChangeQuad 
		HYCO19q<-brm(bf(HYCO~Stage_cm + I(Stage_cm^2), zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYCO20q<-brm(bf(HYCO~1, zi~Stage_cm + I(Stage_cm^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYCO21q<-brm(bf(HYCO~Stage_cm + I(Stage_cm^2), zi~Stage_cm + I(Stage_cm^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		#summary(HYCO19q); waic(HYCO19q)
		#summary(HYCO20q); waic(HYCO20q)
		summary(HYCO21q); waic(HYCO21q)
		#StageChangeMixed
		#HYCO19m<-brm(bf(HYCO~Stage_cm + I(Stage_cm^2), zi~Stage_cm + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYCO20m<-brm(bf(HYCO~Stage_cm, zi~Stage_cm + I(Stage_cm^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		#summary(HYCO19m); waic(HYCO19m)
		summary(HYCO20m); waic(HYCO20m)
		
		HYCO24<-brm(bf(HYCO~MeanQ + I(MeanQ^2), zi~Jul + I(Jul^2) + MeanQ + I(MeanQ^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		summary(HYCO24); waic(HYCO24)

###Save workspace for models 19 to 21 and 24
save.image("U:/!Hydropsyche2018/brms/Rworkspaces/NewStage/brms_HYCO_newSTAGE_19to21and24.RData")

		

#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#

##HYCA

	#Stage Change
		#HYCA19<-brm(bf(HYCA~Stage_cm, zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYCA20<-brm(bf(HYCA~1, zi~Stage_cm + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		#HYCA21<-brm(bf(HYCA~Stage_cm, zi~Stage_cm + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		#summary(HYCA19); waic(HYCA19)
		summary(HYCA20); waic(HYCA20)
		#summary(HYCA21); waic(HYCA21)
		#StageChangeQuad 
		HYCA19q<-brm(bf(HYCA~Stage_cm + I(Stage_cm^2), zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		#HYCA20q<-brm(bf(HYCA~1, zi~Stage_cm + I(Stage_cm^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYCA21q<-brm(bf(HYCA~Stage_cm + I(Stage_cm^2), zi~Stage_cm + I(Stage_cm^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		summary(HYCA19q); waic(HYCA19q)
		#summary(HYCA20q); waic(HYCA20q)
		summary(HYCA21q); waic(HYCA21q)
		#StageChangeMixed
		HYCA19m<-brm(bf(HYCA~Stage_cm + I(Stage_cm^2), zi~Stage_cm + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		HYCA20m<-brm(bf(HYCA~Stage_cm, zi~Stage_cm + I(Stage_cm^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		summary(HYCA19m); waic(HYCA19m)
		summary(HYCA20m); waic(HYCA20m)
###Save workspace for models 19 to 21
save.image("U:/!Hydropsyche2018/brms/Rworkspaces/NewStage/brms_HYCA_newSTAGE_19to21.RData")

#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#

##HCOC

	#Stage Change
		#HCOC19<-brm(bf(HCOC~Stage_cm, zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HCOC20<-brm(bf(HCOC~1, zi~Stage_cm + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		HCOC21<-brm(bf(HCOC~Stage_cm, zi~Stage_cm + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		#summary(HCOC19); waic(HCOC19)
		summary(HCOC20); waic(HCOC20)
		summary(HCOC21); waic(HCOC21)
		#StageChangeQuad 
		#HCOC19q<-brm(bf(HCOC~Stage_cm + I(Stage_cm^2), zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HCOC20q<-brm(bf(HCOC~1, zi~Stage_cm + I(Stage_cm^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		HCOC21q<-brm(bf(HCOC~Stage_cm + I(Stage_cm^2), zi~Stage_cm + I(Stage_cm^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		#summary(HCOC19q); waic(HCOC19q)
		summary(HCOC20q); waic(HCOC20q)
		summary(HCOC21q); waic(HCOC21q)
		#StageChangeMixed
		HCOC19m<-brm(bf(HCOC~Stage_cm + I(Stage_cm^2), zi~Stage_cm + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		HCOC20m<-brm(bf(HCOC~Stage_cm, zi~Stage_cm + I(Stage_cm^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		summary(HCOC19m); waic(HCOC19m)
		summary(HCOC20m); waic(HCOC20m)
###Save workspace for models 19 to 21
save.image("U:/!Hydropsyche2018/brms/Rworkspaces/NewStage/brms_HCOC_newSTAGE_19to21.RData")
	
	
	
#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#

##HYQU

	#Stage Change
		#HYQU19<-brm(bf(HYQU~Stage_cm, zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYQU20<-brm(bf(HYQU~1, zi~Stage_cm + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		#HYQU21<-brm(bf(HYQU~Stage_cm, zi~Stage_cm + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		#summary(HYQU19); waic(HYQU19)
		summary(HYQU20); waic(HYQU20)
		#summary(HYQU21); waic(HYQU21)
		#StageChangeQuad 
		HYQU19q<-brm(bf(HYQU~Stage_cm + I(Stage_cm^2), zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		HYQU20q<-brm(bf(HYQU~1, zi~Stage_cm + I(Stage_cm^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		HYQU21q<-brm(bf(HYQU~Stage_cm + I(Stage_cm^2), zi~Stage_cm + I(Stage_cm^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		summary(HYQU19q); waic(HYQU19q)
		summary(HYQU20q); waic(HYQU20q)
		summary(HYQU21q); waic(HYQU21q)
		#StageChangeMixed
		HYQU19m<-brm(bf(HYQU~Stage_cm + I(Stage_cm^2), zi~Stage_cm + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		HYQU20m<-brm(bf(HYQU~Stage_cm, zi~Stage_cm + I(Stage_cm^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.95, max_treedepth=12))
		summary(HYQU19m); waic(HYQU19m)
		summary(HYQU20m); waic(HYQU20m)
		
		
		HYQU26<-brm(bf(HYQU~Jul + I(Jul^2) + MeanQ, zi~lat + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYQU27<-brm(bf(HYQU~Jul + I(Jul^2) + MeanTemp, zi~lat + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYQU28<-brm(bf(HYQU~Jul + I(Jul^2) + RangeTemp + I(RangeTemp^2), zi~lat + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
			#summary(HYQU26); waic(HYQU26)
			summary(HYQU27); waic(HYQU27)
			summary(HYQU28); waic(HYQU28)
			
			
			
###Save workspace for models 19 to 21
save.image("U:/!Hydropsyche2018/brms/Rworkspaces/NewStage/brms_HYQU_newSTAGE_19to21.RData")





##Figure comparing Jeffs model s6 and Hydropsyche data
 par(mfrow=c(2,2))
 hist(s7$MeanDeltaStage, main="Hydro model")
 hist(HYgc$Stage_cm, main="Hydropsyche data")
 plot(s7$RM,s7$MeanDeltaStage, main="Hydro model")
 plot(HYgc$RiverMile,HYgc$Stage_cm, main="Hydropsyche data")



