
source('meta_analysis_sourcecode.R')

DF <- read.csv('sampledata.csv')

#rename variable
DF$IN_TRIAL=DF$in_trial

###################################################################
###################################################################
#WHEN TREATMENT DOESNT VARY
#*Working Models for when treatment assignment mechanisms does not vary across trials:

#This function only works for 3 trials and a source population, but can be modified for n trials. 
M<-calculate_estimators(DF)

#-------sandwich variance------#
library("geex")
#--------------------------------------------------------------------
#OM
#Calculate SE for the OM estimator

DF$A<-as.numeric(paste(DF$A)) #A can't be coded as a factor
(param_start_OM<-c(unlist(coef(M$Y_mod_treated)),unlist(coef(M$Y_mod_untreated)), unlist(M$ESTIMATES[1])))

OM_mest<-m_estimate(
  estFUN = OM_EE,
  data  = DF,
  root_control = setup_root_control(start = param_start_OM),
  compute_roots = T,
  compute_vcov = T
)

#save variance + SE
param_num_mu<-length(param_start_OM)
OM<-OM_mest@estimates[param_num_mu]
OM_sandwich_se <- diag(OM_mest@vcov)^0.5 
OM_SE<-OM_sandwich_se[param_num_mu]

#---------------------

#IPW

(param_start_IPW<-c(unlist(coef(M$IN_TRIAL_mod)) , unlist(coef(M$A_mod)), unlist(M$ESTIMATES[2]))) 

IPW_mest <-m_estimate(
  estFUN = IPW_EE,
  data  = DF,
  root_control = setup_root_control(start = param_start_IPW),
  compute_roots = T,
  compute_vcov = T
) 

#save variance + SE
param_num_mu<-length(param_start_IPW)
IPW<-IPW_mest@estimates[param_num_mu]
IPW_sandwich_se <- diag(IPW_mest@vcov)^0.5 
IPW_SE<-IPW_sandwich_se[param_num_mu]

#Compare point estimates
print(data.frame(OM, IPW))

#Compare SE from OM and IPW
print(data.frame(OM_SE, IPW_SE))






