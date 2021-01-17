####
# This simulation is report the performance of the true response


#### 0 Simulation set up####
library(parallel)
library(scales)

## 0.1 Generate the seed for the simulation ####

set.seed(2018)
seed_i <- sample(1000000,1000)

## 0.2 Global Parameters ####

ncores <- 30
ProjectName <- "Simulation_INS8"
nsample <- 1000

## 0.3 Functions  ####
GEE_UI <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                  gamma1, gamma, alpha1, alpha0, sigma_e){
  # cat(theta, " \n")
  return(GEE_UfuncIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                      theta[1:3], theta[4:6], sigma = theta[7], xi = theta[8], 
                      gamma1, gamma, alpha1=alpha1, alpha0=alpha0,
                      sigma_e))
}



GEE_GAMMA_IV0 <- function(theta, Y1star, Y2star, Y1, Y2, ValidationMatrix1, ValidationMatrix2, CovMis1, CovMis2, 
                          gamma1, gamma, alpha1, alpha0, sigma_e){
  return(GEE_GAMMAInsIV0(Y1star, Y2star, Y1, Y2,
                         CovMis1, CovMis2, ValidationMatrix1, ValidationMatrix2,
                         beta1=theta[1:3], beta2=theta[4:6], xi=theta[8], sigma=theta[7],
                         gamma1, gamma, alpha1, alpha0, sigma_e, 
                         fixgamma1=0, fixgamma=c(0,0), fixsigma_e=0, 
                         fixalpha1=0, fixalpha0=0)
  )
}


GEE_SIGMA <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, 
                      gamma1, gamma, alpha1, alpha0, sigma_e){
  return(GEE_SIGMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                      theta[1:3], theta[4:6], sigma = theta[7], xi = theta[8], 
                      gamma1, gamma, alpha1, alpha0, sigma_e))
}

GEE_GAMMAIV <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2){
  GAMMA <- GEE_GAMMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, beta1=theta[1:3], beta2=theta[4:6],
                        xi=theta[8], sigma = theta[7])
  return(GAMMA)
}

GEE_GAMMA.inv <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2){
  GAMMA <- GEE_GAMMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, beta1=theta[1:3], beta2=theta[4:6],
                        xi=theta[8], sigma = theta[7])
  GAMMA.inv <- solve(GAMMA)
  return(GAMMA.inv)
}

GEE_cov <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                    gamma1, gamma, alpha1, alpha0, sigma_e){
  GAMMA.inv <- GEE_GAMMA.inv(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
  SIGMA <- GEE_SIGMA(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                     gamma1, gamma, alpha1, alpha0, sigma_e)
  covmatrix <- GAMMA.inv %*% SIGMA %*% t(as.matrix(GAMMA.inv))
  return(covmatrix)
}

##### 1 Implementation Function ####

### A example in deguging
i <- 1
sigma_e <- gamma <- 0.1
alpha <- -1.39

INS_int <- function(i, sigma_e, alpha, gamma){
  ## 1.1 Set up ####
  set.seed(2019)
  seed_i <- sample(1000000,1000)
  set.seed(seed_i[i])
  library(GeneErrorMis)
  library(nleqslv)
  
  ## 1.2 Data Generation ####
  
  ## true parameters
  
  beta1 <- c(0.7,1.5,-1)
  beta2 <- c(0.7,-1.5,1)
  
  # theta <- c(beta1, beta2, sigma, xi)
  sigma <- 1
  rho <- 0
  # gamma <- 0.8
  
  ### Generate the true data sets
  nsample <- 1500
  nvalidation <- 1000
  X <- runif(nsample,-3,4)
  W <- rnorm(nsample,0,sd=1)
  
  
  mu1 <- beta1[1] + beta1[2] * X + beta1[3] * W 
  mu2 <- beta2[1] + beta2[2] * X + beta2[3] * W 
  
  
  expit <- function(x){
    value <- exp(x)/(1+exp(x))
    ifelse(is.na(value),1,value) 
  }
  
  ## Response
  epsilon <- rnorm(nsample,0,1)
  U <- runif(nsample,0,1)
  mu2expit <- expit(mu2)
  
  Y1 <- mu1 +  epsilon
  Y2 <- ifelse(U < mu2expit,1,0)
  
  ## obtain the actuall correlation
  rho <- cor(Y1-mu1,Y2-mu2expit)
  
  ## measurement error and misclassification
  e <- rnorm(nsample,0,sigma_e)
  U2 <- runif(nsample,0,1)
  
  Y1star <- Y1 + gamma * Y2 + e
  Y2star <- ifelse(U2 > expit(alpha),Y2,1-Y2)
  
  ## Naive model
  naive.model1 <- lm(Y1star ~ X + W)
  true.model1 <- lm(Y1 ~ X + W)
  naive.model2 <- glm(Y2star ~ X + W, family = binomial(link = logit))
  true.model2 <- glm(Y2 ~ X + W, family = binomial(link = logit))
  
  
  ## 1.3 Implementation Generation ###
  ## 1.3.1 Preperation ###
  
  DesignMatrix1 <-  DesignMatrix2 <- cbind(rep(1,length(X)),X,W)
  CovMis1 <- cbind(rep(0,length(X)),rep(1,length(X)))
  CovMis2 <- c(rep(1,length(X)))
  
  DesignMatrix1 <- as.matrix(DesignMatrix1)
  DesignMatrix2 <- as.matrix(DesignMatrix2)
  CovMis1 <- as.matrix(CovMis1)
  CovMis2 <- as.matrix (CovMis2)
  
  ## Create the mismeasured data and the validation data
  
  data.mismeasure <- data.frame(Y1star=Y1star[1:(nsample - nvalidation)],Y2star=Y2star[1:(nsample - nvalidation)], DesignMatrix1[1:(nsample - nvalidation),]) 
  data.validation <- data.frame(Y1=Y1[(nsample - nvalidation+1):nsample],Y2=Y2[(nsample - nvalidation+1):nsample], Y1star=Y1star[(nsample - nvalidation+1):nsample],Y2star=Y2star[(nsample - nvalidation+1):nsample],
                                DesignMatrix1[(nsample - nvalidation+1):nsample,],CovMis1[(nsample - nvalidation+1):nsample,],CovMis2[(nsample - nvalidation+1):nsample,]) 

  ## 1.3.2 Prepare different choices of initial variables ###

  beta_Y1_0 <- mean(Y1star)
  beta_Y2_0 <- log(mean(Y2star)/(1-mean(Y2star)))

  intial3 <- c(beta_Y1_0,0,0,beta_Y2_0,0,0,1,0.001)
  intial4 <- c(naive.model1$coefficients,naive.model2$coefficients,1,0)
  
  return(list(seed = seed_i[i],
              naive1coef = naive.model1$coefficients,
              naive1vcov = vcov(naive.model1),
              naive2coef = naive.model2$coefficients,
              naive2vcov = vcov(naive.model2),
              true1coef = true.model1$coefficients,
              true1vcov = vcov(true.model1),
              true2coef = true.model2$coefficients,
              true2vcov = vcov(true.model2)))
  
}


##### 3 Simulation and Results Evaluation ####
### 3.1 Simulation ####
# results_simp <- lapply(1:1000, FUN = INS_noint)
# results_comp <- lapply(1:1000, FUN = INS_int,  sigma_e=0.5, alpha=-1.38, gamma=0)


#### evaluaton: make sure correct  ####
# ### Results Evaluation ###
# 
# beta1 <- c(0.7,1.5,-1,0.7,-1.5,1,1,0)
# 
# results_simp_M<-matrix(unlist(results_simp), ncol = 8, byrow = T)
# 
# bias1 <- colMeans(na.omit(results_simp_M),na.rm = T) - beta1
# 
# results_comp_M<-matrix(unlist(results_comp), ncol = 8, byrow = T)
# 
# bias2 <- colMeans(na.omit(results_comp_M),na.rm = T) - beta1
# 
# xtable(data.frame(bias1,bias2),digits = 3)
# 
# ### 4. Variance Estimation
# 
# try1 <- GEE_SIGMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, 
#              beta1, beta2, xi=0, sigma, gamma1 = 1, gamma=c(0,0.8), alpha1=-1.386294, alpha0=-1.386294,
#              sigma_e = 0.5)
# GEE_SIGMA(NR$x,Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
# M2 <- GEE_GAMMA(NR$x,Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
# GEE_GAMMA.inv(NR$x,Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
# GEE_cov(NR$x,Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
# try2 <- GEE_GAMMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, beta1, beta2, xi, sigma)
# 
# (NR$jac - M2)/NR$jac

### 4.1 Simulation 1: under different degree of measurement error ####
results_1 <- lapply(c(0.1,0.5,0.7), FUN= function(x){
  results_x <- lapply(1:1000, FUN = INS_int,
                         sigma_e = x, alpha = -1.386294, gamma = 0.8)
  return(results_x)
})


truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0)
sigma_e_range <- c(0.1,0.5,0.7)



Results <- NULL

for (k in 1:3) {
  results <- results_1[[k]]
  naive1coef <- NULL
  naive1sd <- NULL
  naive2coef <- NULL
  naive2sd <- NULL
  CI1naive <- NULL
  CI2naive <- NULL
  
  true1coef <- NULL
  true1sd <- NULL
  true2coef <- NULL
  true2sd <- NULL
  CI1true <- NULL
  CI2true <- NULL


  for (i in 1:1000){
    if (is.null(results[[i]])) {
        next}
    naive1coef <- rbind(naive1coef, results[[i]]$naive1coef)
    naive1sd <- rbind(naive1sd,sqrt(diag( results[[i]]$naive1vcov)))
    naive2coef <- rbind(naive2coef, results[[i]]$naive2coef)
    naive2sd <- rbind(naive2sd,sqrt(diag( results[[i]]$naive2vcov)))

    CILBnaive1 <- results[[i]]$naive1coef - 1.96 *(sqrt(diag(results[[i]]$naive1vcov)))
    CIUBnaive1 <- results[[i]]$naive1coef + 1.96 *(sqrt(diag(results[[i]]$naive1vcov)))
    CI1naive <- rbind(CI1naive,ifelse((truebeta[1:3]<CIUBnaive1) & (truebeta[1:3]>CILBnaive1),1,0))
    CILBnaive2 <- results[[i]]$naive2coef - 1.96 *(sqrt(diag(results[[i]]$naive2vcov)))
    CIUBnaive2 <- results[[i]]$naive2coef + 1.96 *(sqrt(diag(results[[i]]$naive2vcov)))
    CI2naive <- rbind(CI2naive,ifelse((truebeta[4:6]<CIUBnaive2) & (truebeta[4:6]>CILBnaive2),1,0))
    
    
    true1coef <- rbind(true1coef, results[[i]]$true1coef)
    true1sd <- rbind(true1sd,sqrt(diag( results[[i]]$true1vcov)))
    true2coef <- rbind(true2coef, results[[i]]$true2coef)
    true2sd <- rbind(true2sd,sqrt(diag( results[[i]]$true2vcov)))
    
    CILBtrue1 <- results[[i]]$true1coef - 1.96 *(sqrt(diag(results[[i]]$true1vcov)))
    CIUBtrue1 <- results[[i]]$true1coef + 1.96 *(sqrt(diag(results[[i]]$true1vcov)))
    CI1true <- rbind(CI1true,ifelse((truebeta[1:3]<CIUBtrue1) & (truebeta[1:3]>CILBtrue1),1,0))
    CILBtrue2 <- results[[i]]$true2coef - 1.96 *(sqrt(diag(results[[i]]$true2vcov)))
    CIUBtrue2 <- results[[i]]$true2coef + 1.96 *(sqrt(diag(results[[i]]$true2vcov)))
    CI2true <- rbind(CI2true,ifelse((truebeta[4:6]<CIUBtrue2) & (truebeta[4:6]>CILBtrue2),1,0))
  }

  biasnaive1 <- colMeans(naive1coef)-truebeta[1:3]
  biasnaive2 <- colMeans(naive2coef)-truebeta[4:6]
  naive_esd <- apply(cbind(naive1coef,naive2coef), MARGIN = 2 , FUN=sd, na.rm=T)
  sdnaive1 <- colMeans(naive1sd)
  sdnaive2 <- colMeans(naive2sd)
  CInaive1 <- colMeans(CI1naive,na.rm=T)
  CInaive2 <- colMeans(CI2naive,na.rm=T)

  naivebias <- c(biasnaive1,biasnaive2,0,0)
  naive_esd <- c(naive_esd,0,0)
  naivesd <- c(sdnaive1,sdnaive2,0,0)
  naiveCI <- c(CInaive1,CInaive2,0,0)

  biastrue1 <- colMeans(true1coef)-truebeta[1:3]
  biastrue2 <- colMeans(true2coef)-truebeta[4:6]
  true_esd <- apply(cbind(true1coef,true2coef), MARGIN = 2 , FUN=sd, na.rm=T)
  sdtrue1 <- colMeans(true1sd)
  sdtrue2 <- colMeans(true2sd)
  CItrue1 <- colMeans(CI1true,na.rm=T)
  CItrue2 <- colMeans(CI2true,na.rm=T)
  
  truebias <- c(biastrue1,biastrue2,0,0)
  true_esd <- c(true_esd,0,0)
  truesd <- c(sdtrue1,sdtrue2,0,0)
  trueCI <- c(CItrue1,CItrue2,0,0)

  Results0 <- data.frame(sigma_e = sigma_e_range[k],
                naivebias=round(naivebias,3),naive_esd=round(naive_esd,3),naivesd=round(naivesd,3),naiveCI=percent(round(naiveCI,3)),
                truebias=round(truebias,3),true_esd=round(true_esd,3),truesd=round(truesd,3),trueCI=percent(round(trueCI,3)))

  Results <- rbind(Results,Results0)
}


library(xtable)
xtable(Results,digits = 3)

### 4.2 Simulation 2: under different degree of misclassification rates ####
results_2 <- lapply(c(-4.595,-2.197,-1.386), FUN= function(x){
  results_x <- lapply(1:1000, FUN = INS_int, 
                      sigma_e = 0.5, alpha = x, gamma = 0.8)
  return(results_x)
})

truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0)
alpha_range <- c(-4.595,-2.197,-1.386)


Results <- NULL

for (k in 1:3) {
  results <- results_1[[k]]
  naive1coef <- NULL
  naive1sd <- NULL
  naive2coef <- NULL
  naive2sd <- NULL
  CI1naive <- NULL
  CI2naive <- NULL
  
  true1coef <- NULL
  true1sd <- NULL
  true2coef <- NULL
  true2sd <- NULL
  CI1true <- NULL
  CI2true <- NULL
  
  
  for (i in 1:1000){
    if (is.null(results[[i]])) {
      next}
    naive1coef <- rbind(naive1coef, results[[i]]$naive1coef)
    naive1sd <- rbind(naive1sd,sqrt(diag( results[[i]]$naive1vcov)))
    naive2coef <- rbind(naive2coef, results[[i]]$naive2coef)
    naive2sd <- rbind(naive2sd,sqrt(diag( results[[i]]$naive2vcov)))
    
    CILBnaive1 <- results[[i]]$naive1coef - 1.96 *(sqrt(diag(results[[i]]$naive1vcov)))
    CIUBnaive1 <- results[[i]]$naive1coef + 1.96 *(sqrt(diag(results[[i]]$naive1vcov)))
    CI1naive <- rbind(CI1naive,ifelse((truebeta[1:3]<CIUBnaive1) & (truebeta[1:3]>CILBnaive1),1,0))
    CILBnaive2 <- results[[i]]$naive2coef - 1.96 *(sqrt(diag(results[[i]]$naive2vcov)))
    CIUBnaive2 <- results[[i]]$naive2coef + 1.96 *(sqrt(diag(results[[i]]$naive2vcov)))
    CI2naive <- rbind(CI2naive,ifelse((truebeta[4:6]<CIUBnaive2) & (truebeta[4:6]>CILBnaive2),1,0))
    
    
    true1coef <- rbind(true1coef, results[[i]]$true1coef)
    true1sd <- rbind(true1sd,sqrt(diag( results[[i]]$true1vcov)))
    true2coef <- rbind(true2coef, results[[i]]$true2coef)
    true2sd <- rbind(true2sd,sqrt(diag( results[[i]]$true2vcov)))
    
    CILBtrue1 <- results[[i]]$true1coef - 1.96 *(sqrt(diag(results[[i]]$true1vcov)))
    CIUBtrue1 <- results[[i]]$true1coef + 1.96 *(sqrt(diag(results[[i]]$true1vcov)))
    CI1true <- rbind(CI1true,ifelse((truebeta[1:3]<CIUBtrue1) & (truebeta[1:3]>CILBtrue1),1,0))
    CILBtrue2 <- results[[i]]$true2coef - 1.96 *(sqrt(diag(results[[i]]$true2vcov)))
    CIUBtrue2 <- results[[i]]$true2coef + 1.96 *(sqrt(diag(results[[i]]$true2vcov)))
    CI2true <- rbind(CI2true,ifelse((truebeta[4:6]<CIUBtrue2) & (truebeta[4:6]>CILBtrue2),1,0))
  }
  
  biasnaive1 <- colMeans(naive1coef)-truebeta[1:3]
  biasnaive2 <- colMeans(naive2coef)-truebeta[4:6]
  naive_esd <- apply(cbind(naive1coef,naive2coef), MARGIN = 2 , FUN=sd, na.rm=T)
  sdnaive1 <- colMeans(naive1sd)
  sdnaive2 <- colMeans(naive2sd)
  CInaive1 <- colMeans(CI1naive,na.rm=T)
  CInaive2 <- colMeans(CI2naive,na.rm=T)
  
  naivebias <- c(biasnaive1,biasnaive2,0,0)
  naive_esd <- c(naive_esd,0,0)
  naivesd <- c(sdnaive1,sdnaive2,0,0)
  naiveCI <- c(CInaive1,CInaive2,0,0)
  
  biastrue1 <- colMeans(true1coef)-truebeta[1:3]
  biastrue2 <- colMeans(true2coef)-truebeta[4:6]
  true_esd <- apply(cbind(true1coef,true2coef), MARGIN = 2 , FUN=sd, na.rm=T)
  sdtrue1 <- colMeans(true1sd)
  sdtrue2 <- colMeans(true2sd)
  CItrue1 <- colMeans(CI1true,na.rm=T)
  CItrue2 <- colMeans(CI2true,na.rm=T)
  
  truebias <- c(biastrue1,biastrue2,0,0)
  true_esd <- c(true_esd,0,0)
  truesd <- c(sdtrue1,sdtrue2,0,0)
  trueCI <- c(CItrue1,CItrue2,0,0)
  
  Results0 <- data.frame(alpha = alpha_range[k],
                         naivebias=round(naivebias,3),naive_esd=round(naive_esd,3),naivesd=round(naivesd,3),naiveCI=percent(round(naiveCI,3)),
                         truebias=round(truebias,3),true_esd=round(true_esd,3),truesd=round(truesd,3),trueCI=percent(round(trueCI,3)))
  
  Results <- rbind(Results,Results0)
}


library(xtable)
xtable(Results,digits = 3)

### 4.3 Simulation 3: under different level of interation ####
results_3 <- lapply(c(0,-0.8,0.8), FUN= function(x){
  results_x <- lapply(1:1000, FUN = INS_int, 
                      sigma_e = 0.5, alpha = -1.386, gamma = x)
  return(results_x)
})

truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0)
gamma_range <- c(0,-0.8,0.8)

Results <- NULL

for (k in 1:3) {
  results <- results_1[[k]]
  naive1coef <- NULL
  naive1sd <- NULL
  naive2coef <- NULL
  naive2sd <- NULL
  CI1naive <- NULL
  CI2naive <- NULL
  
  true1coef <- NULL
  true1sd <- NULL
  true2coef <- NULL
  true2sd <- NULL
  CI1true <- NULL
  CI2true <- NULL
  
  
  for (i in 1:1000){
    if (is.null(results[[i]])) {
      next}
    naive1coef <- rbind(naive1coef, results[[i]]$naive1coef)
    naive1sd <- rbind(naive1sd,sqrt(diag( results[[i]]$naive1vcov)))
    naive2coef <- rbind(naive2coef, results[[i]]$naive2coef)
    naive2sd <- rbind(naive2sd,sqrt(diag( results[[i]]$naive2vcov)))
    
    CILBnaive1 <- results[[i]]$naive1coef - 1.96 *(sqrt(diag(results[[i]]$naive1vcov)))
    CIUBnaive1 <- results[[i]]$naive1coef + 1.96 *(sqrt(diag(results[[i]]$naive1vcov)))
    CI1naive <- rbind(CI1naive,ifelse((truebeta[1:3]<CIUBnaive1) & (truebeta[1:3]>CILBnaive1),1,0))
    CILBnaive2 <- results[[i]]$naive2coef - 1.96 *(sqrt(diag(results[[i]]$naive2vcov)))
    CIUBnaive2 <- results[[i]]$naive2coef + 1.96 *(sqrt(diag(results[[i]]$naive2vcov)))
    CI2naive <- rbind(CI2naive,ifelse((truebeta[4:6]<CIUBnaive2) & (truebeta[4:6]>CILBnaive2),1,0))
    
    
    true1coef <- rbind(true1coef, results[[i]]$true1coef)
    true1sd <- rbind(true1sd,sqrt(diag( results[[i]]$true1vcov)))
    true2coef <- rbind(true2coef, results[[i]]$true2coef)
    true2sd <- rbind(true2sd,sqrt(diag( results[[i]]$true2vcov)))
    
    CILBtrue1 <- results[[i]]$true1coef - 1.96 *(sqrt(diag(results[[i]]$true1vcov)))
    CIUBtrue1 <- results[[i]]$true1coef + 1.96 *(sqrt(diag(results[[i]]$true1vcov)))
    CI1true <- rbind(CI1true,ifelse((truebeta[1:3]<CIUBtrue1) & (truebeta[1:3]>CILBtrue1),1,0))
    CILBtrue2 <- results[[i]]$true2coef - 1.96 *(sqrt(diag(results[[i]]$true2vcov)))
    CIUBtrue2 <- results[[i]]$true2coef + 1.96 *(sqrt(diag(results[[i]]$true2vcov)))
    CI2true <- rbind(CI2true,ifelse((truebeta[4:6]<CIUBtrue2) & (truebeta[4:6]>CILBtrue2),1,0))
  }
  
  biasnaive1 <- colMeans(naive1coef)-truebeta[1:3]
  biasnaive2 <- colMeans(naive2coef)-truebeta[4:6]
  naive_esd <- apply(cbind(naive1coef,naive2coef), MARGIN = 2 , FUN=sd, na.rm=T)
  sdnaive1 <- colMeans(naive1sd)
  sdnaive2 <- colMeans(naive2sd)
  CInaive1 <- colMeans(CI1naive,na.rm=T)
  CInaive2 <- colMeans(CI2naive,na.rm=T)
  
  naivebias <- c(biasnaive1,biasnaive2,0,0)
  naive_esd <- c(naive_esd,0,0)
  naivesd <- c(sdnaive1,sdnaive2,0,0)
  naiveCI <- c(CInaive1,CInaive2,0,0)
  
  biastrue1 <- colMeans(true1coef)-truebeta[1:3]
  biastrue2 <- colMeans(true2coef)-truebeta[4:6]
  true_esd <- apply(cbind(true1coef,true2coef), MARGIN = 2 , FUN=sd, na.rm=T)
  sdtrue1 <- colMeans(true1sd)
  sdtrue2 <- colMeans(true2sd)
  CItrue1 <- colMeans(CI1true,na.rm=T)
  CItrue2 <- colMeans(CI2true,na.rm=T)
  
  truebias <- c(biastrue1,biastrue2,0,0)
  true_esd <- c(true_esd,0,0)
  truesd <- c(sdtrue1,sdtrue2,0,0)
  trueCI <- c(CItrue1,CItrue2,0,0)
  
  Results0 <- data.frame(gamma = gamma_range[k],
                         naivebias=round(naivebias,3),naive_esd=round(naive_esd,3),naivesd=round(naivesd,3),naiveCI=percent(round(naiveCI,3)),
                         truebias=round(truebias,3),true_esd=round(true_esd,3),truesd=round(truesd,3),trueCI=percent(round(trueCI,3)))
  
  Results <- rbind(Results,Results0)
}



library(xtable)
xtable(Results,digits = 3)