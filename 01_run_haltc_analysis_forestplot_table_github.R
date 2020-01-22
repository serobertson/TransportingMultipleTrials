
#complete-case dataset
DF <- read.csv('haltc_mergebasedata_labs_NOTwithmissing.csv')

#renumber trials so 7 is target
DF$S<-ifelse(DF$S==7, 0, DF$S)
DF$S<-ifelse(DF$S==8, 7, DF$S)
DF$S<-ifelse(DF$S==9, 8, DF$S)
DF$S<-ifelse(DF$S==10, 9, DF$S)

table(DF$S)

#********************************************************

#Rename variables of interest:
DF$X1 = DF$PLATELETS
DF$X2 = DF$AGE_RAND
DF$X3 = DF$FEMALE
DF$X4 = DF$PEGINF 
DF$X5 = DF$WHITE
DF$X6 = DF$WBC
DF$X7 = DF$REC_DRUGS
DF$X8 = DF$TRANS_REC
DF$X9 = DF$BMI
DF$X10 = DF$CREATININE
DF$X11 = DF$SMOKESTAT

#create variable
DF$IN_TRIAL<-ifelse(DF$S!=0, 1,0) #variable called: `R' in paper
#-------sandwich variance------#
library("geex")
#--------------------------------------------------------------------
#No transportability analysis:
#Returns unadjusted estimate and adjusted estimates from confounding from OM and IPW in one dataset (one trial or OBS study)"
#--------------------------------------------------------------------
source('0_sourcecode_no_transport.R')

Sresults = list()
#trial 0
Sresults[[1]]<-no_transport_est(DF=DF,trialnum=0,str_covariates="X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11",
                                keep<-c("X1", "X2", "X3", "X4", "X5", "X6","X7","X8","X9","X10","X11"))

#trials 2-8
for(i in 3:9){
  Sresults[[i]]<-no_transport_est(DF=DF,trialnum=i-1,str_covariates="X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11",
                                  keep<-c("X1", "X2", "X3", "X4", "X5", "X6","X7","X8","X9","X10","X11"))
}

#trial 1
Sresults[[2]]<-no_transport_est(DF=DF,trialnum=1,str_covariates="X1+X2+X5+X6+X7+X8+X9+X11",
                                keep<-c("X1", "X2", "X5", "X6","X7","X8","X9","X11"))
#trial 9
Sresults[[10]]<-no_transport_est(DF=DF,trialnum=9,str_covariates="X1+X2+X6+X7+X8+X9+X11",
                                 keep<-c("X1", "X2", "X6","X7","X8","X9","X11"))

final_table<-do.call("rbind", Sresults)
write.csv(final_table, "results_no_transport.csv")

#--------------------------------------------------------------------
#Transporting one-trial-at-a-time to S=0
#--------------------------------------------------------------------

source('0_sourcecode_transport_one_trial.R')

Sresults_transport_one = list()
for(i in 2:8){
  Sresults_transport_one[[i]]<-transport_onetrial_est(DF=DF,trialnum=i,str_covariates="X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11",
                                        keep<-c("X1", "X2", "X3", "X4", "X5", "X6","X7","X8","X9","X10","X11"))
}

Sresults_transport_one[[1]]<-transport_onetrial_est(DF=DF,trialnum=1,str_covariates="X1+X2+X5+X6+X7+X8+X9+X11",
                                      keep<-c("X1", "X2", "X5", "X6","X7","X8","X9","X11"))
Sresults_transport_one[[9]]<-transport_onetrial_est(DF=DF,trialnum=9,str_covariates="X1+X2+X6+X7+X8+X9+X11",
                                      keep<-c("X1", "X2", "X6","X7","X8","X9","X11"))

final_table_transport_one<-do.call("rbind", Sresults_transport_one)

write.csv(final_table_transport_one, "results_transport_onetrial.csv")

#--------------------------------------------------------------------
#Transporting the entire collection of trials to S=0
#--------------------------------------------------------------------
source('0_sourcecode_transport_collection.R')

#Treatment effect (outcome-modeling) estimator
final_tableTE<-transport_collection_trial_TE(DF=DF,str_covariates="X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11",
                                              keep<-c("X1", "X2", "X3", "X4", "X5", "X6","X7","X8","X9","X10","X11"))

#Weighting estimator

final_tableW<-transport_collection_trial_W(DF=DF,str_covariates="X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11",
                                             keep<-c("X1", "X2", "X3", "X4", "X5", "X6","X7","X8","X9","X10","X11"))


final_table_pooled<-rbind(final_tableTE, final_tableW)
write.csv(final_table_pooled, "results_transport_collection.csv")









