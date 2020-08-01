transport_collection_trial_TE<-function(DF=DF,str_covariates,keep){
  
  
  OM<-OM_est_TE(data=DF,str_covariates=str_covariates) 
  DF$p<-OM$p
  
 
  #-------sandwich variance------#
  library("geex")
  
  b1<-DF[,c("Y", "A","S")]
  b2<-DF[,keep]
  
  DFsub<-data.frame(b1, b2)
  
  #OM
  
  (param_start_OM<-c(coef(OM$tx_mod), coef(OM$OMmod),ate=OM$OM))
  
  OM_mest<-m_estimate(
    estFUN =OM_EE_TE,
    data  = DFsub,
    root_control = setup_root_control(start = param_start_OM),
    compute_roots = T,
    compute_vcov = T 
  )
  
  # #return difference (ate) and corresponding standard errors
  OM_ate<-extractEST(geex_output=OM_mest, est_name="ate",param_start=param_start_OM)
  m <- matrix(c(OM_ate[1],OM_ate[2]), ncol = 2, nrow = 1)
  m <- data.frame(m)
  names(m) <- c("EST", "SE")
  m$trialnum<-"all"
  m$ESTNAME<-c("OM")
  rownames(m) <- NULL #print without rownames
  print(m,row.names = FALSE)
  return(m)
  
}



OM_est_TE<-function(data,str_covariates){
  

  #Pr[Z=z|X, S=s*]
  DFtrial<-subset(data, data$S!=0)
  Amod<-as.formula(paste("A ~ ",str_covariates,sep = ""))
  tx_mod<-glm(Amod, data = DFtrial, family = "binomial")
  pred_A1<-predict(tx_mod, type="response", newdata=data)
  data$pred_A1<-pred_A1
  data$W= ((data$A==1)*(pred_A1)^-1 - (1-data$A==1)*(1-pred_A1)^-1)*data$Y
  
  #model specification
  IN_TRIAL_data<-subset(data, S!=0)
  Ymod<-as.formula(paste("W ~ ",str_covariates,sep = ""))
  OMmod<-glm(Ymod, data=IN_TRIAL_data) 
  p<- predict(OMmod,newdata=data, type="response") 
  data$p<-p
  
  S0sub<-subset(data, S==0)
  
  OM<-mean(S0sub$p)
  
  list<-list(OM=OM, p=p,OMmod=OMmod,tx_mod=tx_mod)
  
  return(list)
  
}

OM_EE_TE<- function(data){
  
  A<-data$A
  S<- data$S
  Y <- data$Y
  trialnum<- data$trialnum
  
  
  Xdat<-data.frame(int=1,data) #matrix of only Xs
  Xdat$A<-NULL
  Xdat$S<-NULL
  Xdat$Y<-NULL
  Xdat$trialnum<-NULL
  
  X<-data.matrix(Xdat)
  matA <- cbind(1, data$A) 
  
  A[is.na(A)] <- 0 
  Y[is.na(Y)] <- 0
  
  function(theta){
    
    num_cov<- ncol(X)
    
    #TREATMENT MODEL PIECE
    lp  <- X %*% theta[1: num_cov]
    pa<- plogis(lp)
    score_eqns<-crossprod(X,(S!=0)*(A - pa) )
    
    W= ((A==1)*(pa)^-1 - (1-A==1)*(1-pa)^-1)*Y
    
    #OUTCOME MODEL PIECE 
    beta<-theta[(num_cov+1):(2* num_cov)]
    muate<-theta[(2*num_cov+1)]
    
    m <-X %*% beta
    
    #E[W|X, S=s*]
    ols <-crossprod(X, (S!=0)*(W - m))
    
    ate<-(S==0)*(m-muate) 
    
    c(score_eqns, ols,ate)
    
  }
}


transport_collection_trial_W<-function(DF=DF,str_covariates,keep){
  
  #Fit treatment models Pr[Z=z|X, S, R=1] ~ Pr[Z=z|X,R=1] for our application
  DFtrial<-subset(DF, DF$S!=0)
  Amod<-as.formula(paste("A ~ ",str_covariates,sep = ""))
  tx_mod<-glm(Amod, data = DFtrial, family = "binomial")
  pred_A1<-predict(tx_mod, type="response", newdata=DF)
  DF$pred_A1<-pred_A1
  
  #Pr[S=0|X] 
  S0modform<-as.formula(paste("IN_TRIAL ~ ",str_covariates,sep = ""))
  S0_mod<-glm(S0modform, data = DF, family = "binomial")
  pred_in_trial<-predict(S0_mod, type="response", newdata=DF)
  pred_S0<-1-pred_in_trial
  DF$pred_S0<-pred_S0
  
  w= (DF$IN_TRIAL==1)* pred_S0 * (1-pred_S0)^-1 
  
  IOW1_1 <-(sum(DF$S==0)^-1)* sum(DF$A*(DF$pred_A1)^-1*(DF$IN_TRIAL==1)*w*DF$Y)
  IOW1_0 <-(sum(DF$S==0)^-1)* sum((1-DF$A)*(1-DF$pred_A1)^-1*(DF$IN_TRIAL==1)*w*DF$Y)
  
  IOW1<-((IOW1_1-IOW1_0)) 
  print(IOW1)
  
  (param_start_IPWnew<-c(coef(tx_mod),-1*coef(S0_mod),
                         ate=IOW1))

  b1<-DF[,c("Y", "A","S")]
  b2<-DF[,keep]
  DFsub<-data.frame(b1, b2)
  
  #-------sandwich variance------#
  
  IPW_NEW_mest<-m_estimate(
    estFUN = IPW_EE_W,
    data  = DFsub,
    root_control = setup_root_control(start = param_start_IPWnew),
    compute_roots = T,
    compute_vcov = T
  )
  
  
  IPW_ate<-extractEST(geex_output=IPW_NEW_mest, est_name="ate",param_start=param_start_IPWnew)
  
  m <- matrix(c(IPW_ate[1],IPW_ate[2]), ncol = 2, nrow = 1)
  m <- data.frame(m)
  names(m) <- c("EST", "SE")
  m$trialnum<-"all"
  m$ESTNAME<-c("IPW")
  rownames(m) <- NULL #print without rownames
  print(m,row.names = FALSE)
  return(m)
  
  
}


IPW_EE_W <- function(data){
  
  A<-data$A
  S<- data$S
  Y <- data$Y
  IN_TRIAL<-data$IN_TRIAL
  
  Xdat<-data.frame(int=1,data) #matrix of only Xs
  Xdat$A<-NULL
  Xdat$S<-NULL
  Xdat$Y<-NULL
  X<-data.matrix(Xdat)
  matA <- cbind(1, data$A) 
  
  A[is.na(A)] <- 0 
  Y[is.na(Y)] <- 0
  
  
  function(theta){
    
    num_cov<- ncol(X)
    
    #TREATMENT MODEL PIECE
    lp2  <- X %*% theta[1: num_cov]
    pred_A1<- plogis(lp2)
    score_eqnstx<-crossprod(X,(S!=0)*(A -  pred_A1) )
    
    #TREATMENT MODEL PIECE
    lp  <- X %*% theta[(num_cov+1):(2* num_cov)]
    pred_S0<- plogis(lp)
    score_eqns2<-crossprod(X,((S==0)-  pred_S0) )
    
    mu<-theta[(2* num_cov+1)]
    
    w= (S!=0)* pred_S0 * (1-pred_S0)^-1
    
    IOW1_1 <-A*(pred_A1)^-1*(S!=0)*w*Y
    IOW1_0 <-(1-A)*(1-pred_A1)^-1*(S!=0)*w*Y
    
    
    m<-(IOW1_1-IOW1_0)
    ate<-(m) - (S==0)*mu
    
    c(score_eqnstx, score_eqns2, ate)

    
  }
}



#Function to extract point estimate and SE from geex output
extractEST<-function(geex_output=OM_mest, est_name="m1",param_start=param_start_OM){
  param_num_EST<-match(est_name,names(param_start))
  EST<-geex_output@estimates[param_num_EST]
  
  sandwich_se <- diag(geex_output@vcov)^0.5 
  SE<-sandwich_se[param_num_EST]
  return(c(EST, SE=SE))
}
