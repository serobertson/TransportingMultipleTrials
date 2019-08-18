##Source code for non-varying treatment estimators:

calculate_estimators<-function(DF){
  

  #1) Calculate working models
  
  #a) E[Y|X, IN_TRIAL=1, A=a]
  
  #OUTCOME MODEL FOR TREATED
  DFtrialA1 <-subset(DF, DF$S!=0 & DF$A==1)
  Y_mod_treated <- glm(Y ~ X1 + X2, data=DFtrialA1, family="gaussian")
  g_A1 <-predict(Y_mod_treated, newdata=DF)  #/*g_a(X)*/
  #OUTCOME MODEL FOR UNTREATED
  DFtrialA0 <-subset(DF, DF$S!=0 & DF$A==0)
  Y_mod_untreated <- glm(Y ~ X1 + X2, data=DFtrialA0, family="gaussian")
  g_A0 <-predict(Y_mod_untreated, newdata=DF, type="response")  #/*g_a(X)*/
  
  # b) Pr[A=a|X, IN_TRIAL=1]
  DFtrial<-subset(DF, DF$S!=0)
  A_mod <- glm(A ~ X1 + X2, data=DFtrial, family="binomial")
  e_a1 <- predict(A_mod, newdata=DF, type="response")
  e_a0<-1-e_a1
  
  # c) Pr[IN_TRIAL==1|X]
  IN_TRIAL_mod <- glm(IN_TRIAL ~ X1 + X2, data=DF, family="binomial")
  p <- predict(IN_TRIAL_mod, newdata=DF, type="response")
  
  #2) Calculations for estimators
  
  
  S0_size <- sum(DF$S==0)
  
  #-----------------------------
  #estimators for treated
  #OM1
  term1 = (DF$S == 0) * (g_A1)
  ESTIMATE_OM1=sum(term1)/S0_size
  
  #IPW1
  w_1= ((1-p)/(p*e_a1))
  term2_mod=(DF$S != 0) *(DF$A==1) * w_1 * (DF$Y)
  ESTIMATE_IPW1=sum(term2_mod)/S0_size
  #-----------------------------
  #estimators for untreated
  #OM0
  term1 = (DF$S == 0) * (g_A0)
  ESTIMATE_OM0=sum(term1)/S0_size
  
  #IPW0
  w_0= ((1-p)/(p*e_a0))
  term2_mod=(DF$S != 0) *(DF$A==0) * w_0 * (DF$Y)
  ESTIMATE_IPW0=sum(term2_mod)/S0_size
  #-----------------------------
  #ATE
  IPW=ESTIMATE_IPW1-ESTIMATE_IPW0
  OM=ESTIMATE_OM1-ESTIMATE_OM0
  
  ESTIMATES<-data.frame(OM, IPW)
  
  #return working models (useful for starting values for geex) and estimates
  newList <- list("Y_mod_treated"=Y_mod_treated,"Y_mod_untreated"=Y_mod_untreated,
                  "A_mod"=A_mod, "IN_TRIAL_mod"=IN_TRIAL_mod, "ESTIMATES"=ESTIMATES  )
  
  return(newList)
  
}



#estimating equation functions 

OM_EE <- function(data){
  
  A<-data$A
  S<- data$S
  Y <- data$Y
  IN_TRIAL<-data$IN_TRIAL
  
  X <- cbind(1, data$X1, data$X2) 
  matA <- cbind(1, data$A) 
  
  A[is.na(A)] <- 0 
  Y[is.na(Y)] <- 0
  
  function(theta){
    
    #OUTCOME MODEL PIECE 
    beta1<-theta[1:3]
    beta0<-theta[4:6]
    mu<-theta[7]
    
    m_A1 <-X %*% beta1

    #E[Y|X,A=1]
    ols_A1 <-crossprod(X, (IN_TRIAL*A)*(Y - m_A1))
    
    #mean1 <- (1-IN_TRIAL)*(m_A1-mu1) 
    
    #E[Y|X,A=0]
    m_A0 <-X %*% beta0
    ols_A0 <-crossprod(X, (IN_TRIAL*(1-A))*(Y - m_A0))
    
    #mean0 <- (1-IN_TRIAL)*(m_A0-mu0)
    
    mean<-(1-IN_TRIAL)*(m_A1-m_A0-mu)
    
    c(ols_A1,ols_A0, mean)
  }
}


IPW_EE <- function(data){
  
  A<-data$A
  S<- data$S
  Y <- data$Y
  IN_TRIAL<-data$IN_TRIAL
  
  X <- cbind(1, data$X1, data$X2) 
  matA <- cbind(1, data$A) 
  
  A[is.na(A)] <- 0 
  Y[is.na(Y)] <- 0 
  
  
  function(theta){
    
    mu_mean<-theta[7]
    
    #PARTICIPATION MODEL PIECE
    lp  <- X %*% theta[1:3]
    ps <- plogis(lp)
    score_eqns<-crossprod(X, IN_TRIAL-ps)
    
    #TREATMENT MODEL PIECE
    lp2  <- X %*% theta[4:6] 
    pa<- plogis(lp2)
    score_eqns2<-crossprod(X,IN_TRIAL*(A - pa) )
    
    w = (A * IN_TRIAL*(1-ps))/(ps*pa) + ((1 - A)*IN_TRIAL*(1-ps))/(ps*(1-pa)) 
    
    #means
    summand1<-(w*A*IN_TRIAL*Y)  #term2
    summand0<-(w*(1-A)*IN_TRIAL*Y)  
    mean<-(summand1-summand0)- (1-IN_TRIAL)*mu_mean
    c(score_eqns,score_eqns2, mean)
    
    
  }
}



#########################################
###Estimators for when treatment varies



calculate_estimators_tx_vary<-function(DF){
  
  DFtrial<-subset(DF, S!=0)
  
  #TREATMENT MODEL 
  A_mods <- lapply(split(DFtrial, DFtrial$S), 'glm', formula = A ~X1 + X2, family="binomial")
  A1_S <- lapply(A_mods,'predict', newdata=DF, type="response")
  
  #OUTCOME MODEL FOR TREATED
  DFtrialA1 <-subset(DF, DF$S!=0 & DF$A==1)
  Y_mod_treated <- lapply(split(DFtrialA1, DFtrialA1$S), 'glm', formula = Y ~ X1 + X2, family="gaussian")
  Y_A1 <- lapply( Y_mod_treated,'predict', newdata=DF)
  #OUTCOME MODEL FOR UNTREATED
  DFtrialA0 <-subset(DF, DF$S!=0 & DF$A==0)
  Y_mod_untreated <- lapply(split(DFtrialA0, DFtrialA0$S), 'glm', formula = Y ~ X1 + X2, family="gaussian")
  Y_A0 <- lapply( Y_mod_untreated,'predict', newdata=DF)
  
  
  #S regression
  library(mlogit)
  df2 <- mlogit.data(DFtrial, choice="S", shape="wide")
  mlogit.mod <- mlogit(S ~ 1 | X1 + X2, data=df2)
  (mlogit.cf <- coef(mlogit.mod))

  #Calculate probability of being in each trial over all individuals:
  beta<-matrix(mlogit.cf,nrow = 2,ncol = 3)
  beta<-t(beta)
  X<-cbind(1,DF$X1, DF$X2, DF$X3)
  m_A1 <-X %*% beta
  denom<-1+exp(m_A1[,1]) + exp(m_A1[,2])
  p2<-exp(m_A1[,1])/denom
  p3<-exp(m_A1[,2])/denom
  p1<-1-(p2+p3)  ##reference category is base level 
  probabilties<-cbind(p1, p2,p3) #Strial
  
  
  ####
  #Model for in_trial
  #Pr[IN_TRIAL==1|X]
  IN_TRIAL_mod <- glm(IN_TRIAL ~ X1 + X2, data=DF, family="binomial")
  p <- predict(IN_TRIAL_mod, newdata=DF, type="response")
  
  #final_prod = prod((S==i)*Y_A1[[i]])
  prod1=(DF$S==1)*Y_A1[[1]]
  prod2=(DF$S==2)*Y_A1[[2]]
  prod3=(DF$S==3)*Y_A1[[3]]
  
  DF$prod1<-prod1
  DF$prod2<-prod2
  DF$prod3<-prod3
  
  g_aS <- rowSums(DF[c("prod1", "prod2", "prod3")], na.rm = TRUE)
  
  #newprod
  DF$newprod1=p1*Y_A1[[1]]
  DF$newprod2=p2*Y_A1[[2]]
  DF$newprod3=p3*Y_A1[[3]]
  
  summand_g_as_ps <- rowSums(DF[c("newprod1", "newprod2","newprod3")], na.rm = TRUE)
  
  #term3summand
  DF$term3summand1= (DF$S==1)*A1_S[[1]] 
  DF$term3summand2= (DF$S==2)*A1_S[[2]] 
  DF$term3summand3= (DF$S==3)*A1_S[[3]] 
  #calculations for estimators:
  
  S0_size <- sum(DF$S==0)
  
  #OM1
  term1<-(DF$S == 0) * (summand_g_as_ps)
  ESTIMATE_OM1<-sum(term1)/S0_size
  
  #IPW1
  
  #/*calculate Pr[A=a|X, S] == e_a*/
  e_a<-rowSums(DF[c("term3summand1", "term3summand2","term3summand3")], na.rm = TRUE)
  inv_e_a<-ifelse(e_a==0, 0, (e_a)^-1)
  
  term3_mod<-(DF$S != 0) * (DF$A==1) * ((1-p)/(p))*(inv_e_a)*(DF$Y)
  ESTIMATE_IPW1<-sum(term3_mod)/S0_size
  
  #**********
  #UNTREATED calculations:
  #final_prod = prod((S==i)*Y_A1[[i]])
  prod1=(DF$S==1)*Y_A0[[1]]
  prod2=(DF$S==2)*Y_A0[[2]]
  prod3=(DF$S==3)*Y_A0[[3]]
  
  DF$prod1<-prod1
  DF$prod2<-prod2
  DF$prod3<-prod3
  
  g_aS <- rowSums(DF[c("prod1", "prod2", "prod3")], na.rm = TRUE)
  
  #newprod
  DF$newprod1=p1*Y_A0[[1]]
  DF$newprod2=p2*Y_A0[[2]]
  DF$newprod3=p3*Y_A0[[3]]
  
  summand_g_as_ps <- rowSums(DF[c("newprod1", "newprod2","newprod3")], na.rm = TRUE)
  
  #term3summand
  DF$term3summand1= (DF$S==1)*(1-A1_S[[1]])
  DF$term3summand2= (DF$S==2)*(1-A1_S[[2]]) 
  DF$term3summand3= (DF$S==3)*(1-A1_S[[3]]) 
  
  #untreated
  #OM0
  term1<-(DF$S == 0) * (summand_g_as_ps)
  ESTIMATE_OM0<-sum(term1)/S0_size
  
  #IPW0
  #/*calculate Pr[A=a|X, S] == e_a*/
  e_a<-rowSums(DF[c("term3summand1", "term3summand2","term3summand3")], na.rm = TRUE)
  inv_e_a<-ifelse(e_a==0, 0, (e_a)^-1)
  
  term3_mod<-(DF$S != 0) * (DF$A==0) * ((1-p)/(p))*(inv_e_a)*(DF$Y)
  ESTIMATE_IPW0<-sum(term3_mod)/S0_size
  
  ESTIMATE_OM<-ESTIMATE_OM1-ESTIMATE_OM0
  ESTIMATE_IPW<-ESTIMATE_IPW1-ESTIMATE_IPW0
  ESTIMATES<-data.frame(ESTIMATE_OM, ESTIMATE_IPW)
  
  #return working models (useful for starting values for geex) and estimates
  newList <- list("Y_mod_treated"=Y_mod_treated,"Y_mod_untreated"=Y_mod_untreated,
                  "A_mods"=A_mods, "IN_TRIAL_mod"=IN_TRIAL_mod, "mlogit.mod"=mlogit.mod, "ESTIMATES"=ESTIMATES  )
  
  return(newList)
  
}



#estimating equation functions 

OM_EE_vary <- function(data){
  
  A<-data$A
  S<- data$S
  Y <- data$Y
  IN_TRIAL<-data$IN_TRIAL
  
  X <- cbind(1, data$X1, data$X2) 
  matA <- cbind(1, data$A) 
  
  A[is.na(A)] <- 0 
  Y[is.na(Y)] <- 0
  
  dummy2 <- as.numeric(data$S== 2)
  dummy3 <- as.numeric(data$S== 3)
  dummies<-cbind(dummy2, dummy3)
  
  function(theta){
    
    #OUTCOME MODEL PIECE FOR TREATED
    OM1_beta1<-theta[1:3]
    OM2_beta1<-theta[4:6]
    OM3_beta1<-theta[7:9]
    
    g1_A1 <-X %*% OM1_beta1
    g1_A2 <-X %*% OM2_beta1
    g1_A3 <-X %*% OM3_beta1

    #E[Y|X,A=1, S=s]
    ols1_A1 <-crossprod(X, ((S==1)*A)*(Y - g1_A1))
    ols1_A2 <-crossprod(X, ((S==2)*A)*(Y - g1_A2))
    ols1_A3 <-crossprod(X, ((S==3)*A)*(Y - g1_A3))
    
    #OUTCOME MODEL PIECE FOR UNTREATED
    OM1_beta0<-theta[10:12]
    OM2_beta0<-theta[13:15]
    OM3_beta0<-theta[16:18]
    
    g0_A1 <-X %*% OM1_beta0
    g0_A2 <-X %*% OM2_beta0
    g0_A3 <-X %*% OM3_beta0
    
    #E[Y|X,A=0, S=s]
    ols0_A1 <-crossprod(X, ((S==1)*(1-A))*(Y - g0_A1))
    ols0_A2 <-crossprod(X, ((S==2)*(1-A))*(Y - g0_A2))
    ols0_A3 <-crossprod(X, ((S==3)*(1-A))*(Y - g0_A3))
    
    #SELECTION MODELS (multinomial model)
    
    gamma<-theta[19:21] #can beta be a matrix instead of column vector...
    m_A1 <-X %*% gamma  #get predictions
    
    alpha<-theta[22:24]
    m_A2 <-X %*% alpha 
    
    denom<-1+exp(m_A1) + exp(m_A2)
    p2<-exp(m_A1)/denom
    p3<-exp(m_A2)/denom
    # p1<-1-(p2+p3)  ##reference category is base level 
    
    probs<-cbind(p2, p3)
    
    score_eqns4<- crossprod(X,IN_TRIAL*(dummy2 - p2) )
    score_eqns5<- crossprod(X,IN_TRIAL*(dummy3 - p3) )
    
    ##calculations for estimator for TREATED
    p1<- 1-(p2+p3)
    newprod1=p1*g1_A1
    newprod2=p2*g1_A2
    newprod3=p3*g1_A3
    
    summand_g1_as_ps<- newprod1+newprod2+newprod3
    
    #untreated
    newprod1=p1*g0_A1
    newprod2=p2*g0_A2
    newprod3=p3*g0_A3
    
    summand_g0_as_ps<- newprod1+newprod2+newprod3
    
    
    mu<-theta[25]
    
    mean<-(1-IN_TRIAL)*((summand_g1_as_ps-summand_g0_as_ps)- mu) 
    
    c(ols1_A1,ols1_A2, ols1_A3,ols0_A1,ols0_A2, ols0_A3, score_eqns4, score_eqns5, mean)
    
    
    
    
  }
}



IPW_EE_vary <- function(data){
  
  A<-data$A
  S<- data$S
  Y <- data$Y
  IN_TRIAL<-data$IN_TRIAL
  
  X <- cbind(1, data$X1, data$X2) 
  matA <- cbind(1, data$A) 
  
  A[is.na(A)] <- 0 
  Y[is.na(Y)] <- 0 
  
  
  function(theta){
    
    A1_beta<-theta[1:3]
    A2_beta<-theta[4:6]
    A3_beta<-theta[7:9]
    
    g_A1 <-plogis(X %*% A1_beta)
    g_A2 <- plogis(X %*% A2_beta)
    g_A3 <- plogis(X %*% A3_beta)

    #E[Y|S=s]
    ols_A1 <-crossprod(X, ((S==1))*(A - g_A1))
    ols_A2 <-crossprod(X, ((S==2))*(A - g_A2))
    ols_A3 <-crossprod(X, ((S==3))*(A - g_A3))
    
    #term3summand treated
    #probability of being in each trial * probability of being untreated
    term3summand1_tx= (S==1)*g_A1
    term3summand2_tx= (S==2)*g_A2 
    term3summand3_tx= (S==3)*g_A3 
    
    #Model for being in the trial
    IN_TRIAL_beta<-theta[10:12]
    g_IN_TRIAL <-plogis(X %*% IN_TRIAL_beta)
    ols_IN_TRIAL <-crossprod(X, (IN_TRIAL - g_IN_TRIAL))
    
    ##calculate estimand
    mu<-theta[13]
    #summand for treated:
    e_a1<-term3summand1_tx + term3summand2_tx + term3summand3_tx
    inv_e_a1<-ifelse(e_a1==0, 0, (e_a1)^-1)
    
    p<-g_IN_TRIAL
    term3_mod_tx<-(S != 0) * (A==1) * ((1-p)/(p))*(inv_e_a1)*(Y)
    
    #term3summand untreated
    term3summand1_untx= (S==1)*(1-g_A1)
    term3summand2_untx= (S==2)*(1-g_A2) 
    term3summand3_untx= (S==3)*(1-g_A3) 
    
    e_a0<-term3summand1_untx + term3summand2_untx + term3summand3_untx
    inv_e_a0<-ifelse(e_a0==0, 0, (e_a0)^-1)
    term3_mod_untx<-(S != 0) * (A==0) * ((1-p)/(p))*(inv_e_a0)*(Y)
    
    
    mean<-(term3_mod_tx-term3_mod_untx)- (1-IN_TRIAL)*mu

    c(ols_A1, ols_A2, ols_A3, ols_IN_TRIAL, mean)
    
    
    
    
  }
}


