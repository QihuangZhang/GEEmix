####
# This simulation is testing whether the first version of Insertion method has unbiased property.


#### 0 Simulation set up####
library(parallel)

## 0.1 Generate the seed for the simulation ####

set.seed(2018)
seed_i <- sample(1000000,1000)

## 0.2 Global Parameters ####

ncores <- 30
ProjectName <- "Simulation_INS1"
nsample <- 1000

## 0.3 Functions  ####
GEE_UI <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2){
  # cat(theta, " \n")
  return(GEE_UfuncIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                      theta[1:3], theta[4:6], phi = theta[7], xi = theta[8], 
                      gamma1 = 1, gamma=c(0,0), alpha1=-1.386294, alpha0=-1.386294,
                      sigma_e = 0.5))
}

GEE_UI2 <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2){
  # cat(theta, " \n")
  # if ((theta[8]>1)|(theta[8]<-1)) {theta[8] <- 0}
  return(GEE_UfuncIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                      theta[1:3], theta[4:6], phi = theta[7], xi = theta[8], 
                      gamma1 = 1, gamma=c(0,0.8), alpha1=-1.386294, alpha0=-1.386294,
                      sigma_e = 0.5))
}

GEE_SIGMA <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,gamma){
  return(GEE_SIGMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                      theta[1:3], theta[4:6], phi = theta[7], xi = theta[8], 
                      gamma1 = 1, gamma, alpha1=-1.386294, alpha0=-1.386294,
                      sigma_e = 0.5))
}

GEE_GAMMA <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2){
  GAMMA <- GEE_GAMMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, beta1=theta[1:3], beta2=theta[4:6],
                        xi=theta[8], phi = theta[7])
  return(GAMMA)
}

GEE_GAMMA.inv <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2){
  GAMMA <- GEE_GAMMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, beta1=theta[1:3], beta2=theta[4:6],
                        xi=theta[8], phi = theta[7])
  GAMMA.inv <- solve(GAMMA)
  return(GAMMA.inv)
}

GEE_cov <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,gamma){
  GAMMA.inv <- GEE_GAMMA.inv(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
  SIGMA <- GEE_SIGMA(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,gamma)
  covmatrix <- GAMMA.inv %*% SIGMA %*% t(as.matrix(GAMMA.inv))
  return(covmatrix)
}

#### 1 Implementation Function: No Interatction in Measurement Error ####

INS_noint <- function(i){
  
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
  
  # theta <- c(beta1, beta2, phi, xi)
  phi <- 1
  rho <- 0
  
  ### Generate the true data sets
  nsample <- 1000
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
  e <- rnorm(nsample,0,0.5)
  U2 <- runif(nsample,0,1)
  
  Y1star <- Y1 + e
  Y2star <- ifelse(U2>0.2,Y2,1-Y2)
  
  ## Naive model
  naive.model1 <- lm(Y1star ~ X + W)
  true.model1 <- lm(Y1 ~ X + W)
  naive.model2 <- glm(Y2star ~ X + W, family = binomial(link = logit))
  true.model2 <- glm(Y2 ~ X + W, family = binomial(link = logit))
  
  ## 1.2 Implementation of the methods #####
  ## 1.2.1 Preparation ###
  
  DesignMatrix1 <-  DesignMatrix2 <- cbind(rep(1,length(X)),X,W)
  CovMis1 <- cbind(rep(0,length(X)),rep(0,length(X)))
  CovMis2 <- c(rep(1,length(X)))
  
  DesignMatrix1 <- as.matrix(DesignMatrix1)
  DesignMatrix2 <- as.matrix(DesignMatrix2)
  CovMis1 <- as.matrix(CovMis1)
  CovMis2 <- as.matrix (CovMis2)
  
  ## 1.2.2 Implementation ###

  beta_Y1_0 <- mean(Y1star)
  beta_Y2_0 <- log(mean(Y2star)/(1-mean(Y2star)))
  
  intial3 <- c(beta_Y1_0,0,0,beta_Y2_0,0,0,1,0.001)
  intial4 <- c(naive.model1$coefficients,naive.model2$coefficients,1,0)
  
  NR <- nleqslv(intial4,GEE_UI,Y1star=Y1star, Y2star=Y2star, DesignMatrix1=DesignMatrix1, DesignMatrix2=DesignMatrix2,
          CovMis1=CovMis1, CovMis2=CovMis2,jacobian=T,control=list(maxit=1000))
  
  betahat <- ifelse(abs(NR$x)<10,NR$x,NA)
  
  if (!any(is.na(betahat))) {
    cov <- GEE_cov(betahat,Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,c(0,0))
    sd <- sqrt(diag(cov))}
  else {
    sd <- rep(NA,length(betahat))
  }
  
  # if (!any(is.na(betahat))) {
  #   SIGMA <- GEE_SIGMA(betahat,Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
  #   GAMMA <- NR$jac
  #   GAMMA.inv <- solve(GAMMA,tol=1e-250)
  #   cov <-  GAMMA.inv %*% SIGMA %*% t(as.matrix(GAMMA.inv))
  #   sd <- sqrt(diag(cov))}
  # else {
  #   sd <- rep(NA,length(betahat))
  # }
  
  return(list(seed = seed_i[i],
              naive1coef = naive.model1$coefficients,
              naive1vcov = vcov(naive.model1),
              naive2coef = naive.model2$coefficients,
              naive2vcov = vcov(naive.model2),
              betahat=betahat,
              sd=sd))
}


##### 2 Implementation Function: The correlated mismeasruement ####

INS_int <- function(i){
  ## 2.1 Set up ####
  set.seed(2019)
  seed_i <- sample(1000000,1000)
  set.seed(seed_i[i])
  library(GeneErrorMis)
  library(nleqslv)
  
  ## 2.2 Data Generation ####
  
  ## true parameters
  
  beta1 <- c(0.7,1.5,-1)
  beta2 <- c(0.7,-1.5,1)
  
  # theta <- c(beta1, beta2, phi, xi)
  phi <- 1
  rho <- 0
  gamma <- 0.8
  
  ### Generate the true data sets
  nsample <- 1000
  X <- runif(nsample,-3,4)
  W <- rnorm(nsample,0,sd=1)
  
  
  mu1 <- beta1[1] + beta1[2] * X + beta1[3] * W 
  mu2 <- beta2[1] + beta2[2] * X + beta2[3] * W 
  
  
  expit <- function(x){
    value <- exp(x)/(1+exp(x))
    ifelse(is.na(value),1,value) 
  }
  
  logit <- function(x){
    value <- log(x/(1-x))
    return(value)
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
  e <- rnorm(nsample,0,0.5)
  U2 <- runif(nsample,0,1)
  
  Y1star <- Y1 + gamma * Y2 + e
  Y2star <- ifelse(U2>0.2,Y2,1-Y2)
  
  ## Naive model
  naive.model1 <- lm(Y1star ~ X + W)
  true.model1 <- lm(Y1 ~ X + W)
  naive.model2 <- glm(Y2star ~ X + W, family = binomial(link = logit))
  true.model2 <- glm(Y2 ~ X + W, family = binomial(link = logit))
  
  
  ## 2.3 Implementation Generation ###
  ## 2.3.1 Preperation ###
  
  DesignMatrix1 <-  DesignMatrix2 <- cbind(rep(1,length(X)),X,W)
  CovMis1 <- cbind(rep(0,length(X)),rep(1,length(X)))
  CovMis2 <- c(rep(1,length(X)))
  
  DesignMatrix1 <- as.matrix(DesignMatrix1)
  DesignMatrix2 <- as.matrix(DesignMatrix2)
  CovMis1 <- as.matrix(CovMis1)
  CovMis2 <- as.matrix (CovMis2)
  
  ## 2.3.2 Implementation ###
  intial4 <- c(naive.model1$coefficients,naive.model2$coefficients,1,0)
  
  NR <- nleqslv(intial4,GEE_UI2,Y1star=Y1star, Y2star=Y2star, DesignMatrix1=DesignMatrix1, DesignMatrix2=DesignMatrix2,
                CovMis1=CovMis1, CovMis2=CovMis2 ,control=list(maxit=2000))
  
  betahat <- ifelse(abs(NR$x)<10,NR$x,NA)
  
  if (!any(is.na(betahat))) {
    cov <- GEE_cov(betahat,Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,c(0,0.8))
    sd <- sqrt(diag(cov))}
  else {
    sd <- rep(NA,length(betahat))
  }
  
  return(list(seed = seed_i[i],
              naive1coef = naive.model1$coefficients,
              naive1vcov = vcov(naive.model1),
              naive2coef = naive.model2$coefficients,
              naive2vcov = vcov(naive.model2),
              betahat=betahat,
              sd=sd))
}


##### 3 Simulation and Results Evaluation ####
### 3.1 Simulation ###
results_simp <- lapply(1:1000,FUN = INS_noint)
results_comp <- lapply(1:1000,FUN = INS_int)

# ### Results Evaluation ####
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
#              beta1, beta2, xi=0, phi, gamma1 = 1, gamma=c(0,0.8), alpha1=-1.386294, alpha0=-1.386294,
#              sigma_e = 0.5)
# GEE_SIGMA(NR$x,Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
# M2 <- GEE_GAMMA(NR$x,Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
# GEE_GAMMA.inv(NR$x,Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
# GEE_cov(NR$x,Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
# try2 <- GEE_GAMMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, beta1, beta2, xi, phi)
# 
# (NR$jac - M2)/NR$jac

### 4.1 Simulation 1: No interaction ####
results_simp <- lapply(1:1000,FUN = INS_noint)

truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0)

naive1coef <- NULL
naive1sd <- NULL
naive2coef <- NULL
naive2sd <- NULL
CI1naive <- NULL
CI2naive <- NULL
betas <- NULL
sds <- NULL
CIs <- NULL
results <- results_simp

for (i in 1:1000){
  naive1coef <- rbind(naive1coef, results[[i]]$naive1coef)
  naive1sd <- rbind(naive1sd,sqrt(diag( results[[i]]$naive1vcov)))
  naive2coef <- rbind(naive2coef, results[[i]]$naive2coef)
  naive2sd <- rbind(naive2sd,sqrt(diag( results[[i]]$naive2vcov)))
  
  betahat0 <- results[[i]]$betahat
  sd0 <- results[[i]]$sd
  sd0 <- ifelse(abs(sd0)<10,sd0,NA)
  betas <- rbind(betas, betahat0)
  sds <- rbind(sds, sd0)
  
  CILBnaive1 <- results[[i]]$naive1coef - 1.96 *(sqrt(diag(results[[i]]$naive1vcov)))
  CIUBnaive1 <- results[[i]]$naive1coef + 1.96 *(sqrt(diag(results[[i]]$naive1vcov)))
  CI1naive <- rbind(CI1naive,ifelse((truebeta[1:3]<CIUBnaive1) & (truebeta[1:3]>CILBnaive1),1,0))
  CILBnaive2 <- results[[i]]$naive2coef - 1.96 *(sqrt(diag(results[[i]]$naive2vcov)))
  CIUBnaive2 <- results[[i]]$naive2coef + 1.96 *(sqrt(diag(results[[i]]$naive2vcov)))
  CI2naive <- rbind(CI2naive,ifelse((truebeta[4:6]<CIUBnaive2) & (truebeta[4:6]>CILBnaive2),1,0))
  CILB <- betahat0 - 1.96 *(sd0)
  CIUB <- betahat0 + 1.96 *(sd0)
  CIs <- rbind(CIs,ifelse((truebeta<as.vector(CIUB)) & (truebeta>as.vector(CILB)),1,0))
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

bias1 <- colMeans(na.omit(betas),na.rm = T) - truebeta
sd_emp <- apply(na.omit(betas),MARGIN = 2, FUN = sd)
sd_mod <- colMeans(na.omit(sds),na.rm = T)
CIrate <- colMeans(na.omit(CIs),na.rm = T)

Results <- data.frame(naivebias,naive_esd,naivesd,naiveCI,bias1,sd_emp,sd_mod,CIrate)

library(xtable)
xtable(Results)

### 4.2 Simulation 2: interaction 
results_comp <- lapply(1:1000,FUN = INS_int)

truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0)

betas <- NULL
sds <- NULL
CIs <- NULL
results <- results_comp

for (i in 1:1000){
  betahat0 <- results[[i]]$betahat
  sd0 <- results[[i]]$sd
  # sd0 <- ifelse(abs(sd0)<10,sd0,NA)
  betas <- rbind(betas, betahat0)
  sds <- rbind(sds, sd0)
  CILB <- betahat0 - 1.96 *(sd0)
  CIUB <- betahat0 + 1.96 *(sd0)
  CIs <- rbind(CIs,ifelse((truebeta<as.vector(CIUB)) & (truebeta>as.vector(CILB)),1,0))
}

bias1 <- colMeans(na.omit(betas),na.rm = T) - truebeta
sd_emp <- apply(na.omit(betas),MARGIN = 2, FUN = sd)
sd_mod <- colMeans(na.omit(sds),na.rm = T)
CIrate <- colMeans(na.omit(CIs),na.rm = T)

