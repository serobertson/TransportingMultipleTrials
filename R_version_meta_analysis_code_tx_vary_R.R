source('meta_analysis_sourcecode.R')

DF <- read.csv('sampledata.csv')
DF$IN_TRIAL=DF$in_trial
###################################################################
###################################################################


H<-calculate_estimators_tx_vary(DF)

####geex

library("geex")

mlogit.cf<-unlist(coef(H$mlogit.mod))
start_beta<-c(mlogit.cf[1],mlogit.cf[3],mlogit.cf[5],mlogit.cf[2],mlogit.cf[4],mlogit.cf[6])

#--------------------------------------------------------------------
#OM

DF$A<-as.numeric(paste(DF$A)) #A can't be coded as a factor

(param_start_OM<-c(unlist(coef(H$Y_mod_treated[[1]])), 
                   unlist(coef(H$Y_mod_treated[[2]])),
                   unlist(coef(H$Y_mod_treated[[3]])), 
                   unlist(coef(H$Y_mod_untreated[[1]])), 
                   unlist(coef(H$Y_mod_untreated[[2]])),
                   unlist(coef(H$Y_mod_untreated[[3]])),
                   start_beta, 
                   unlist(H$ESTIMATES$ESTIMATE_OM)))

OM_mest<-m_estimate(
  estFUN = OM_EE_vary,
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

(param_start_IPW<-c(unlist(coef(H$A_mods[[1]])), unlist(coef(H$A_mods[[2]]))
                     ,unlist(coef(H$A_mods[[3]])),
                     unlist(coef(H$IN_TRIAL_mod)), unlist(H$ESTIMATES$ESTIMATE_IPW)))

IPW_mest <-m_estimate(
  estFUN = IPW_EE_vary,
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
print(c(OM_SE, IPW_SE))

