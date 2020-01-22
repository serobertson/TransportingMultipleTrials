

transport_onetrial_est<-function(DF=DF, trialnum=1,str_covariates,keep){
  
  DFsub<-subset(DF, S==0 | S==trialnum) 
  DFsub$S<-ifelse(DFsub$S!=0, 1, 0)
  
  Smod<-as.formula(paste("S ~ ",str_covariates,sep = ""))
  Amod<-as.formula(paste("A ~ ",str_covariates,sep = ""))
  
  weights<-generate_weights(Smod=Smod,
                            Amod=Amod, 
                            data=DFsub) 
  
  DF2<-weights$dat

  
  OM<-OM_est(data=DF2,str_covariates=str_covariates)
  DF2$p1<-OM$p1
  DF2$p0<-OM$p0
  
  IOW2<-IOW2_est(data=DF2)
  
  #Estimates of ATE
  results<-rbind(OM=OM$OM, IOW2=IOW2$IOW2)
  print(results)
  
  #-------sandwich variance------#
  
  library("geex")
  
  b1<-DFsub[,c("Y", "A","S")]
  b2<-DFsub[,keep]
  
  DFsub<-data.frame(b1, b2)
  
  #OM
  (param_start_OM<-c(coef(OM$OM1mod), coef(OM$OM0mod), 
                     m1=OM$OM_1, m0=OM$OM_0,ate=OM$OM))
  
  OM_mest<-m_estimate(
    estFUN = OM_EE,
    data  = DFsub,
    root_control = setup_root_control(start = param_start_OM),
    compute_roots = T,
    compute_vcov = T
  )
  
  #return difference (ate) and corresponding standard errors
  OM_ate<-extractEST(geex_output=OM_mest, est_name="ate",param_start=param_start_OM)
  
  
  #IOW2 sandwich variance
  
  (param_start_IOW2<-c(coef(weights$Smod),coef(weights$Amod),
                       int1=coef(IOW2$IOW1mod)["(Intercept)"],int0=coef(IOW2$IOW0mod)["(Intercept)"],
                       m1=IOW2$IOW2_1, m0=IOW2$IOW2_0, ate=IOW2$IOW2))
  
  IOW2_mest <-m_estimate(
    estFUN = IOW2_EE,
    data  = DFsub,
    root_control = setup_root_control(start = param_start_IOW2),
    compute_roots = T,
    compute_vcov = T
  ) 
  
  #save variance + SE
  
  IOW2_ate<-extractEST(geex_output=IOW2_mest, est_name="ate",param_start=param_start_IOW2)
  
  print("Returns transportability estimates from OM and IPW one trial `n'` at a time to POP `0' (i.e., trial 1 to target 0, 
      trial 2 to target 0, etc)")
  
  m <- matrix(c(OM_ate[1], IOW2_ate[1], OM_ate[2],IOW2_ate[2]), ncol = 2, nrow = 2)
  m <- data.frame(m)
  names(m) <- c("EST", "SE")
  m$trialnum<-trialnum
  m$ESTNAME<-c("OM", "IPW")
  rownames(m) <- NULL #print without rownames
  print(m,row.names = FALSE)
  return(m)
  
  
}



generate_weights<-function(Smod,Amod, data){
  
  S1data<-subset(data, S==1)
  
  w_reg<-glm(Smod, family="binomial", data=data)
  ps<- predict(w_reg,newdata=data, type="response") 
  w_reg2<-glm(Amod, family="binomial", data=S1data)
  pa<- predict(w_reg2,newdata=data, type="response") 
  
  w= (data$A*data$S*(1-ps) )/(ps*pa) + ((1 -data$A) *data$S*(1-ps) ) /(ps*(1-pa))
  data$w<-w
  list<-list(dat=data, Smod=w_reg, Amod=w_reg2)
  return(list)
}

OM_est<-function(data,str_covariates){
  
  #model specification
  Ymod<-as.formula(paste("Y ~ ",str_covariates,sep = ""))

  S1data_A1<-subset(data, S==1 & A==1)
  OM1mod<-glm(Ymod, data=S1data_A1) #X4 deleted
  p1<- predict(OM1mod,newdata=data, type="response") 
  data$p1<-p1
  
  S1data_A0<-subset(data, S==1 & A==0)
  OM0mod<-glm( Ymod, data=S1data_A0)
  p0<- predict(OM0mod,newdata=data, type="response") 
  data$p0<-p0
  
  S0sub<-subset(data, S==0)
  
  OM_1<-mean(S0sub$p1)
  OM_0<-mean(S0sub$p0)
  OM<-mean(S0sub$p1)-mean(S0sub$p0)
  list<-list(OM_1=OM_1, OM_0=OM_0, OM=OM, p1=p1, p0=p0, OM1mod=OM1mod, OM0mod=OM0mod)
  return(list)
  
}


IOW1_est<-function(data){
  
  A<-data$A
  S<-data$S
  w<-data$w
  Y<-data$Y
  
  IOW1_1 <-(sum((1-S))^-1)* sum(A*S*w*Y)
  IOW1_0 <-(sum((1-S))^-1)* sum((1-A)*S*w*Y)
  IOW1 = IOW1_1 - IOW1_0
  
  return(list(IOW1_1=IOW1_1,IOW1_0=IOW1_0, IOW1=IOW1))
  
}

IOW2_est<-function(data){
  S0data<-subset(data, S==0)
  
  S1data_A1<-subset(data, S==1 & A==1)
  IOW1mod<-glm(formula=Y~1, data=S1data_A1, weights=w)
  p1<- predict(IOW1mod,newdata=S0data, type="response") 
  
  S1data_A0<-subset(data, S==1 & A==0)
  IOW0mod<-glm(formula=Y~1, data=S1data_A0, weights=w)
  p0<- predict(IOW0mod,newdata=S0data, type="response") 
  
  IOW2_1<-mean(p1)
  IOW2_0<-mean(p0)
  IOW2<-mean(p1)-mean(p0)
  
  list<-list(IOW2_1=IOW2_1,IOW2_0=IOW2_0, IOW2=IOW2,IOW1mod=IOW1mod,IOW0mod=IOW0mod)
  return(list)
}


#functions for geex 

OM_EE <- function(data){
  
  A<-data$A
  S<- data$S
  Y <- data$Y
  
  Xdat<-data.frame(int=1,data) #matrix of only Xs
  Xdat$A<-NULL
  Xdat$S<-NULL
  Xdat$Y<-NULL
  
  X<-data.matrix(Xdat)
  matA <- cbind(1, data$A) 
  
  A[is.na(A)] <- 0 
  Y[is.na(Y)] <- 0
  
  function(theta){
    
    #OUTCOME MODEL PIECE 
    num_cov<- ncol(X)
    
    beta<-theta[1: num_cov]
    alpha<-theta[(num_cov+1):(2* num_cov)]
    mu1<-theta[(2* num_cov+1)]
    mu0<-theta[(2* num_cov+2)]
    muate<-theta[(2* num_cov+3)]
    
    m_A1 <-X %*% beta
    m_A0<-X %*% alpha
    
    #E[Y|X, A=1]
    ols_A1 <-crossprod(X, (S*A)*(Y - m_A1))
    
    #E[Y|X, A=0]
    ols_A0 <-crossprod(X, (S*(1-A))*(Y - m_A0))
    
    mean1<-(1-S)*(m_A1-mu1) 
    mean0 <- (1-S)*(m_A0-mu0) 
    ate<-(1-S)*(m_A1-m_A0-muate) 
    
    c(ols_A1,ols_A0,mean1, mean0,ate)
    
    
  }
}




IOW2_EE <- function(data){
  
  A<-data$A
  S<- data$S
  Y <- data$Y
  
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
    
    #PARTICIPATION MODEL PIECE
    lp  <- X %*% theta[1: num_cov]
    ps <- plogis(lp)
    score_eqns<-crossprod(X, S-ps)
    
    #TREATMENT MODEL PIECE
    lp2  <- X %*% theta[(num_cov+1):(2* num_cov)]
    pa<- plogis(lp2)
    score_eqns2<-crossprod(X,S*(A - pa) )
    
    w = (A * S*(1-ps))/(ps*pa) + ((1 - A)*S*(1-ps))/(ps*(1-pa)) 
    
    #OUTCOME MODEL PIECE for weighted E[Y|X, A=1]
    m_A1<-1 %*% theta[(2* num_cov+1)]
    m_A0<-1 %*% theta[(2* num_cov+2)]
    
    linear_eqns1<-crossprod(1, (S*A*w)*(Y -  m_A1) )
    
    #OUTCOME MODEL PIECE for weighted E[Y|X, A=0]
    linear_eqns0<-crossprod(1, (S*(1-A)*w)*(Y - m_A0) )
    
    mu1<-theta[(2* num_cov+3)]
    mu0<-theta[(2* num_cov+4)]
    muate<-theta[(2* num_cov+5)]
    
    mean1 <- (1-S)*( m_A1 -mu1) 
    mean0 <- (1-S)*( m_A0 -mu0)
    ate <- (1-S)*(m_A1- m_A0 - muate)
    
    c(score_eqns,score_eqns2, linear_eqns1,linear_eqns0,  mean1,  mean0,ate)
    
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
