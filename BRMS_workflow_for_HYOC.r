##BRMS Bayesian regression models for data describing Hydropsyche captured in light traps in Colorado River Basin
##Example code for one species - Hydropsyche occidentalis
##loading brms package requires RTOOLS / C++ capability
library(brms)
library(rstan)
library(ggplot2)

#speed up processing, always run before running models
rstan_options (auto_write=TRUE)
options (mc.cores=parallel::detectCores ())

#load dataset
setwd("P:/BIOLOGICAL/Foodbase/LIGHT_TRAPS/Upper Basin/Hydropsyche")
HY5<-read.csv("HY5_Nov2018.csv")
head(HY5)

###Standardize certain columns of data ("Z scores") in a new dataframe ("HY")
HY<-HY5
#Z score ENV vars
	HY$Jul<-scale(as.numeric(HY$Jul))
	HY$DamDistKM<-scale(as.numeric(HY$DamDistKM))
	HY$MeanQ<-scale(HY$MeanQ)
	HY$SC_gradient<-scale(HY$SC_gradient)
	HY$lat<-scale(as.numeric(HY$lat))
	HY$MeanTemp<-scale(HY$MeanTemp)
	HY$RangeTemp<-scale(HY$RangeTemp)
head(HY)

#####Run models for Hydropsyche occidentalis using stepwise model selection#####
#I planned models in a seperate excel spreadsheet prior to model runs and kept track of wAIC and model parameters there
#Models that would not converge at adapt_delta=0.8 and max_treedepth=10 (default settings) were bumped up to a_d=0.9 and mtd=11, then a_d=0.95 and mtd=12.
#If model did not converge at maximum flexibility settings, it was excluded from further model selection
#Models take a lot of CPU. Running across multiple cores saves time. I ran models on a remote server with 10 cores.

#Set working directory for saving workspaces
#Save models intermittently in small batches to not overload workspaces
setwd("U:/!Hydropsyche2018/brms/Rworkspaces/HYOC")


##Base model for Hydropsyche occidentalis
	HYOC0<-brm(bf(HYOC~1, zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
	summary(HYOC0); waic(HYOC0)

#TIER ONE - SINGLE VARIABLE MODELS#
	##Latitude
		HYOC1<-brm(bf(HYOC~lat, zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC2<-brm(bf(HYOC~1, zi~lat + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC3<-brm(bf(HYOC~lat, zi~lat + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC1); waic(HYOC1)
		summary(HYOC2); waic(HYOC2)
		summary(HYOC3); waic(HYOC3)
		##LatQuad
		HYOC1q<-brm(bf(HYOC~lat + I(lat^2), zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC2q<-brm(bf(HYOC~1, zi~lat + I(lat^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC3q<-brm(bf(HYOC~lat + I(lat^2), zi~lat + I(lat^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC1q); waic(HYOC1q)
		summary(HYOC2q); waic(HYOC2q)
		summary(HYOC3q); waic(HYOC3q)
		##LatMixed
		HYOC1m<-brm(bf(HYOC~lat + I(lat^2), zi~lat + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC2m<-brm(bf(HYOC~lat, zi~lat + I(lat^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC1m); waic(HYOC1m)   
		summary(HYOC2m); waic(HYOC2m) 

	#Julian Day
		HYOC4<-brm(bf(HYOC~Jul, zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC5<-brm(bf(HYOC~1, zi~Jul + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC6<-brm(bf(HYOC~Jul, zi~Jul + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC4); waic(HYOC4)
		summary(HYOC5); waic(HYOC5)
		summary(HYOC6); waic(HYOC6)
		#JulQuad
		HYOC4q<-brm(bf(HYOC~Jul + I(Jul^2), zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC5q<-brm(bf(HYOC~1, zi~Jul + I(Jul^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC6q<-brm(bf(HYOC~Jul + I(Jul^2), zi~Jul + I(Jul^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC4q); waic(HYOC4q)
		summary(HYOC5q); waic(HYOC5q)
		summary(HYOC6q); waic(HYOC6q)
		#JulMixed
		HYOC4m<-brm(bf(HYOC~Jul + I(Jul^2), zi~Jul + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC5m<-brm(bf(HYOC~Jul, zi~Jul + I(Jul^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC4m); waic(HYOC4m)  
		summary(HYOC5m); waic(HYOC5m) 

###Save workspace for models 0 to 6
save.image("brms_HYOC_0to6.RData")


	#DamDist
		HYOC7<-brm(bf(HYOC~DamDistKM, zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8)) #25 divergent transitions with a_d=0.8 HYOC7.8<-HYOC7
		HYOC8<-brm(bf(HYOC~1, zi~DamDistKM + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC9<-brm(bf(HYOC~DamDistKM, zi~DamDistKM + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC7); waic(HYOC7)
		summary(HYOC8); waic(HYOC8)
		summary(HYOC9); waic(HYOC9)
		#DamDistQUAD
		HYOC7q<-brm(bf(HYOC~DamDistKM + I(DamDistKM^2), zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC8q<-brm(bf(HYOC~1, zi~DamDistKM + I(DamDistKM^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC9q<-brm(bf(HYOC~DamDistKM+ I(DamDistKM^2), zi~DamDistKM+ I(DamDistKM^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC7q); waic(HYOC7q)
		summary(HYOC8q); waic(HYOC8q)
		summary(HYOC9q); waic(HYOC9q)
		#DamDistMixed
		HYOC7m<-brm(bf(HYOC~DamDistKM + I(DamDistKM^2), zi~DamDistKM + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10)) 
		HYOC8m<-brm(bf(HYOC~DamDistKM, zi~DamDistKM + I(DamDistKM^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=11)) 
		summary(HYOC7m); waic(HYOC7m)
		summary(HYOC8m); waic(HYOC8m)

	#MeanTemp
		HYOC10<-brm(bf(HYOC~MeanTemp, zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC11<-brm(bf(HYOC~1, zi~MeanTemp + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC12<-brm(bf(HYOC~MeanTemp, zi~MeanTemp + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC10); waic(HYOC10)
		summary(HYOC11); waic(HYOC11)
		summary(HYOC12); waic(HYOC12)
		#MeanTempQuad
		HYOC10q<-brm(bf(HYOC~MeanTemp + I(MeanTemp^2), zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC11q<-brm(bf(HYOC~1, zi~MeanTemp + I(MeanTemp^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC12q<-brm(bf(HYOC~MeanTemp + I(MeanTemp^2), zi~MeanTemp + I(MeanTemp^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC10q); waic(HYOC10q)
		summary(HYOC11q); waic(HYOC11q)
		summary(HYOC12q); waic(HYOC12q)
		#MeanTempMixed
		HYOC10m<-brm(bf(HYOC~MeanTemp + I(MeanTemp^2), zi~MeanTemp + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC11m<-brm(bf(HYOC~MeanTemp, zi~MeanTemp + I(MeanTemp^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC10m); waic(HYOC10m)
		summary(HYOC11m); waic(HYOC11m)
		
###Save workspace for models 7 to 12
save.image("brms_HYOC_7to12.RData")

	#RangeTemp
		HYOC13<-brm(bf(HYOC~RangeTemp, zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC14<-brm(bf(HYOC~1, zi~RangeTemp + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC15<-brm(bf(HYOC~RangeTemp, zi~RangeTemp + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC13); waic(HYOC13)
		summary(HYOC14); waic(HYOC14)
		summary(HYOC15); waic(HYOC15)
		#RangeTempQUAD
		HYOC13q<-brm(bf(HYOC~RangeTemp + I(RangeTemp^2), zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC14q<-brm(bf(HYOC~1, zi~RangeTemp + I(RangeTemp^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC15q<-brm(bf(HYOC~RangeTemp + I(RangeTemp^2), zi~RangeTemp + I(RangeTemp^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC13q); waic(HYOC13q)
		summary(HYOC14q); waic(HYOC14q)
		summary(HYOC15q); waic(HYOC15q)
		#RangeTempMixed
		HYOC13m<-brm(bf(HYOC~RangeTemp + I(RangeTemp^2), zi~RangeTemp + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC14m<-brm(bf(HYOC~RangeTemp, zi~RangeTemp + I(RangeTemp^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=11))
		summary(HYOC13m); waic(HYOC13m)
		summary(HYOC14m); waic(HYOC14m)


	#MeanQ
		HYOC16<-brm(bf(HYOC~MeanQ, zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC17<-brm(bf(HYOC~1, zi~MeanQ + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC18<-brm(bf(HYOC~MeanQ, zi~MeanQ + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC16); waic(HYOC16)
		summary(HYOC17); waic(HYOC17)
		summary(HYOC18); waic(HYOC18)
		#MeanQ QUAD
		HYOC16q<-brm(bf(HYOC~MeanQ + I(MeanQ^2), zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC17q<-brm(bf(HYOC~1, zi~MeanQ + I(MeanQ^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8,  max_treedepth=10)) 
		HYOC18q<-brm(bf(HYOC~MeanQ + I(MeanQ^2), zi~MeanQ + I(MeanQ^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC16q); waic(HYOC16q)
		summary(HYOC17q); waic(HYOC17q)
		summary(HYOC18q); waic(HYOC18q)
		#MeanQMixed 
		HYOC16m<-brm(bf(HYOC~MeanQ + I(MeanQ^2), zi~MeanQ + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC17m<-brm(bf(HYOC~MeanQ, zi~MeanQ + I(MeanQ^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC16m); waic(HYOC16m)
		summary(HYOC17m); waic(HYOC17m)
		
###Save workspace for models 13 to 18
save.image("brms_HYOC_13to18.RData")

	#Stage Change
		HYOC19<-brm(bf(HYOC~SC_gradient, zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC20<-brm(bf(HYOC~1, zi~SC_gradient + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10, max_treedepth=10))
		HYOC21<-brm(bf(HYOC~SC_gradient, zi~SC_gradient + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC19); waic(HYOC19)
		summary(HYOC20); waic(HYOC20)
		summary(HYOC21); waic(HYOC21)
		#StageChangeQuad 
		HYOC19q<-brm(bf(HYOC~SC_gradient + I(SC_gradient^2), zi~1 + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC20qb<-HYOC20q
		HYOC21qb<-HYOC21q
		HYOC20q<-brm(bf(HYOC~1, zi~SC_gradient + I(SC_gradient^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC21q<-brm(bf(HYOC~SC_gradient + I(SC_gradient^2), zi~SC_gradient + I(SC_gradient^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC19q); waic(HYOC19q)
		summary(HYOC20q); waic(HYOC20q)
		summary(HYOC21q); waic(HYOC21q)
		#StageChangeMixed
		HYOC19m<-brm(bf(HYOC~SC_gradient + I(SC_gradient^2), zi~SC_gradient + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		HYOC20m<-brm(bf(HYOC~SC_gradient, zi~SC_gradient + I(SC_gradient^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
		summary(HYOC19m); waic(HYOC19m)
		summary(HYOC20m); waic(HYOC20m)

###Save workspace for models 19 to 21
save.image("brms_HYOC_19to21.RData")

#TIER TWO - TWO VARIABLE MODELS#
	HYOC22<-brm(bf(HYOC~MeanTemp + I(MeanTemp^2) + lat, zi~MeanTemp + I(MeanTemp^2) + lat + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
	HYOC23<-brm(bf(HYOC~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2), zi~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
	HYOC24<-brm(bf(HYOC~MeanTemp + I(MeanTemp^2) + RangeTemp + I(RangeTemp^2), zi~MeanTemp + I(MeanTemp^2) + RangeTemp + I(RangeTemp^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
	HYOC25<-brm(bf(HYOC~MeanTemp + I(MeanTemp^2) + MeanQ + I(MeanQ^2), zi~MeanTemp + I(MeanTemp^2) + MeanQ + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
	HYOC26<-brm(bf(HYOC~MeanTemp + I(MeanTemp^2) + SC_gradient, zi~MeanTemp + I(MeanTemp^2) + SC_gradient + I(SC_gradient^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10))
	summary(HYOC22); waic(HYOC22)
	summary(HYOC23); waic(HYOC23)
	summary(HYOC24); waic(HYOC24)
	summary(HYOC25); waic(HYOC25)
	summary(HYOC26); waic(HYOC26)
save.image("brms_HYOC_22to26.RData")


#TIER THREE - THREE VARIABLE MODELS#
	HYOC27<-brm(bf(HYOC~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2) + lat, zi~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2) + lat + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10, max_treedepth=10))
	HYOC28<-brm(bf(HYOC~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2) + RangeTemp + I(RangeTemp^2), zi~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2) + RangeTemp + I(RangeTemp^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10, max_treedepth=10))
	HYOC29<-brm(bf(HYOC~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2) + MeanQ + I(MeanQ^2), zi~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2) + MeanQ + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10, max_treedepth=10))
	HYOC30<-brm(bf(HYOC~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2) + SC_gradient, zi~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2) + SC_gradient + I(SC_gradient^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.8, max_treedepth=10, max_treedepth=10))
	summary(HYOC27); waic(HYOC27)
	summary(HYOC28); waic(HYOC28)
	summary(HYOC29); waic(HYOC29)
	summary(HYOC30); waic(HYOC30)
save.image("brms_HYOC_27to33.RData")


#TIER FOUR - FOUR VARIABLE MODELS#
	HYOC31<-brm(bf(HYOC~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2) + SC_gradient + lat, zi~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2) + SC_gradient + I(SC_gradient^2) + lat + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.85, max_treedepth=12))
	HYOC32<-brm(bf(HYOC~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2) + SC_gradient + RangeTemp + I(RangeTemp^2), zi~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2) + SC_gradient + I(SC_gradient^2) + RangeTemp + I(RangeTemp^2) + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.85, max_treedepth=12))
	HYOC33<-brm(bf(HYOC~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2) + SC_gradient + MeanQ + I(MeanQ^2), zi~MeanTemp + I(MeanTemp^2) + Jul + I(Jul^2) + SC_gradient + I(SC_gradient^2) + MeanQ + (1|SEGCODE)), data=HY, family=zero_inflated_negbinomial(link="log", link_shape="log", link_zi="logit"), control=list(adapt_delta=0.85, max_treedepth=12))
	summary(HYOC31); waic(HYOC31)
	summary(HYOC32); waic(HYOC32)
	summary(HYOC33); waic(HYOC33)
save.image("brms_HYOC_27to33.RData")




###PLOT THE BEST MODEL (HYOC30)
##Load Multiplot, a function for displaying multiple inset plots via ggplot2
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


##SCATTER PLOT OF RAW DATA
#Legend overlayed later in InkScape

#New dataframe for plotting
	HY3<-HY5
	
	##Set zeroes below 0 so they are visible
	HY3$HYOC<-ifelse(HY3$HYOC<1,-100, HY3$HYOC+5)
	HY3$HYOS<-ifelse(HY3$HYOS<1,-100, HY3$HYOS+5)
	HY3$HYCO<-ifelse(HY3$HYCO<1,-100, HY3$HYCO+5)
	HY3$HYCA<-ifelse(HY3$HYCA<1,-100, HY3$HYCA+5)
	HY3$HYQU<-ifelse(HY3$HYQU<1,-100, HY3$HYQU+5)
	HY3$HCOC<-ifelse(HY3$HCOC<1,-100, HY3$HCOC+5)
	HY3$black0<--3

	##Set stage change and Q to desired metrics
	HY3$Stage_cm<-HY2$SC_gradient*30.48 # change to cm
	HY3$MeanQ_m3s<-HY2$MeanQ * 0.028316846592 #change to m3s
	
	head(HY3)

##Scatterplots for all 6 species of Hydropsyche compared to environmental variables

	#Stage change
		ggStage<- ggplot(aes(Stage_cm, HYOC), data=HY3) + xlab("Stage change (cm)") + 
		ylab("Hydropsyche catch rate (#/hr)") + ylim(-12,600) +
		geom_point(aes(Stage_cm,black0),  col="grey", alpha=0.6, cex=1)+
		geom_point(aes(Stage_cm,HYOS),  col="springgreen1", alpha=0.6, cex=1)+
		geom_point(aes(Stage_cm,HYOC),  col="navyblue", alpha=0.6, cex=1)+
		geom_point(aes(Stage_cm,HYCO),  col="yellow", alpha=0.6, cex=1)+
		geom_point(aes(Stage_cm,HYQU),  col="magenta", alpha=0.6, cex=1)+
		geom_point(aes(Stage_cm,HYCA),  col="red", alpha=0.6, cex=1)+
		geom_point(aes(Stage_cm,HCOC),  col="orange", alpha=0.6, cex=1)+
		theme_bw()
		ggStage
		
	#MeanTemp
		ggmt<- ggplot(aes(MeanTemp, HYOC), data=HY3) + xlab("Mean Temp (°C)") + 
		ylab("Hydropsyche catch rate (#/hr)") + ylim(-12,600) +
		geom_point(aes(MeanTemp,black0),  col="grey", alpha=0.6, cex=1)+
		geom_point(aes(MeanTemp,HYOS),  col="springgreen1", alpha=0.6, cex=1)+
		geom_point(aes(MeanTemp,HYOC),  col="navyblue", alpha=0.6, cex=1)+
		geom_point(aes(MeanTemp,HYCO),  col="yellow", alpha=0.6, cex=1)+
		geom_point(aes(MeanTemp,HYQU),  col="magenta", alpha=0.6, cex=1)+
		geom_point(aes(MeanTemp,HYCA),  col="red", alpha=0.6, cex=1)+
		geom_point(aes(MeanTemp,HCOC),  col="orange", alpha=0.6, cex=1)+
		theme_bw()
		ggmt

	#RangeTemp
		ggrt<- ggplot(aes(RangeTemp, HYOC), data=HY3) + xlab("Range Temp (°C)") + 
		ylab("Hydropsyche catch rate (#/hr)") + ylim(-12,600) +
		geom_point(aes(RangeTemp,black0),  col="grey", alpha=0.6, cex=1)+
		geom_point(aes(RangeTemp,HYOS),  col="springgreen1", alpha=0.6, cex=1)+
		geom_point(aes(RangeTemp,HYOC),  col="navyblue", alpha=0.6, cex=1)+
		geom_point(aes(RangeTemp,HYCO),  col="yellow", alpha=0.6, cex=1)+
		geom_point(aes(RangeTemp,HYQU),  col="magenta", alpha=0.6, cex=1)+
		geom_point(aes(RangeTemp,HYCA),  col="red", alpha=0.6, cex=1)+
		geom_point(aes(RangeTemp,HCOC),  col="orange", alpha=0.6, cex=1)+
		theme_bw()
		ggrt

	#Mean discharge
		ggmq<- ggplot(aes(MeanQ_m3s, HYOC), data=HY3) + xlab(bquote("Mean discharge (m3s)")) + 
		ylab("Hydropsyche catch rate (#/hr)") + ylim(-12,600) +
		geom_point(aes(MeanQ_m3s,black0),  col="grey", alpha=0.6, cex=1)+
		geom_point(aes(MeanQ_m3s,HYOS),  col="springgreen1", alpha=0.6, cex=1)+
		geom_point(aes(MeanQ_m3s,HYOC),  col="navyblue", alpha=0.6, cex=1)+
		geom_point(aes(MeanQ_m3s,HYCO),  col="yellow", alpha=0.6, cex=1)+
		geom_point(aes(MeanQ_m3s,HYQU),  col="magenta", alpha=0.6, cex=1)+
		geom_point(aes(MeanQ_m3s,HYCA),  col="red", alpha=0.6, cex=1)+
		geom_point(aes(MeanQ_m3s,HCOC),  col="orange", alpha=0.6, cex=1)+
		theme_bw()
		ggmq
#View all 4 plots together
multiplot(ggStage,ggmq,ggmt,ggrt,cols=2)




###Plot the BRMS output
#Tutorial: http://www.flutterbys.com.au/stats/tut/tut7.2b.html

##Load R workspace that contains best model (HYOC30)
load("brms_HYOC_27to33.RData")


#Export data out of marginal_effects into a usable format for GGplot, one variable at a time
##Stage
	#Export data from model
		summary(marginal_effects(HYOC30))
		HYOCsc_man<-marginal_effects(HYOC30)$SC_gradient
		head(HYOCsc_man)
	#Plot exported model output with error
	#Manually Relabel the X axis with calculations so it is not Z scored (Z = (x-mean(x))/sd(x)
		meanSC<-mean(HY5$SC_gradient)
		sdSC<-sd(HY5$SC_gradient)
		#ie: What is the Z score of actual value 2 cm of stage change?
			(2-meanSC)/sdSC #1.104972

	STAGE<-ggplot() + 
		geom_line(data=HYOCsc_man, aes(y = estimate__, x = SC_gradient)) + 
		geom_ribbon(data=HYOCsc_man, aes(x = SC_gradient, ymin = lower__,ymax = upper__), fill = "springgreen1", alpha = 0.6) + 
		scale_y_continuous("Modeled abundance") + scale_x_continuous("Stage Change (cm)", breaks=c(-2.311031, -0.6030295,1.104972, 2.812973), labels=c("0", "1", "2","3")) +
		theme_bw()
	STAGE 


########MULTIPLOT SCATTER + MODELS#############
multiplot(ggStage,STAGE,cols=2)

















