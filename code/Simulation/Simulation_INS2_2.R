####
# This simulation is forked from Version 2 with following changes:

#   1:  Adding the covariates to the misclassification model



#### 0 Simulation set up ####
library(parallel)
library(scales)

## 0.1 Generate the seed for the simulation ####

set.seed(2019)
seed_i <- sample(1000000,1000)

## 0.2 Global Parameters ####

ncores <- 30
ProjectName <- "Simulation_INS2"
nsample <- 1000

## 0.3 Functions  ####
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

GEE_UI <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                  gamma1, gamma, alpha1, alpha0, sigma_e){
  # cat(theta, " \n")
  return(GEE_UfuncIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                      theta[1:3], theta[4:6], sigma = theta[7], xi = theta[8], 
                      gamma1, gamma, alpha1=alpha1, alpha0=alpha0,
                      sigma_e))
}


GEE_SIGMA <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, 
                      gamma1, gamma, alpha1, alpha0, sigma_e){
  return(GEE_SIGMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                      theta[1:3], theta[4:6], sigma = theta[7], xi = theta[8], 
                      gamma1, gamma, alpha1, alpha0, sigma_e))
}

GEE_GAMMA <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2){
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

#### 1 Implementation Function: No Interatction in Measurement Error ####
i <- 1
sigma_e <- 0.1


# INS_noint <- function(i,sigma_e,alpha){
#   
#   ## 1.1 Set up ####
#   set.seed(2019)
#   seed_i <- sample(1000000,1000)
#   set.seed(seed_i[i])
#   library(GeneErrorMis)
#   library(nleqslv)
#   
#   ## 1.2 Data Generation ####
#   
#   ## true parameters
#   
#   beta1 <- c(0.7,1.5,-1)
#   beta2 <- c(0.7,-1.5,1)
#   
#   # theta <- c(beta1, beta2, sigma, xi)
#   sigma <- 1
#   rho <- 0
#   
#   ### Generate the true data sets
#   nsample <- 1000
#   X <- runif(nsample,-3,4)
#   W <- rnorm(nsample,0,sd=1)
#   
#   
#   mu1 <- beta1[1] + beta1[2] * X + beta1[3] * W 
#   mu2 <- beta2[1] + beta2[2] * X + beta2[3] * W 
#   
#   
#   expit <- function(x){
#     value <- exp(x)/(1+exp(x))
#     ifelse(is.na(value),1,value) 
#   }
#   
#   ## Response
#   epsilon <- rnorm(nsample,0,1)
#   U <- runif(nsample,0,1)
#   mu2expit <- expit(mu2)
#   
#   Y1 <- mu1 +  epsilon
#   Y2 <- ifelse(U < mu2expit,1,0)
#   
#   ## obtain the actuall correlation
#   rho <- cor(Y1-mu1,Y2-mu2expit)
#   
#   ## measurement error and misclassification
#   e <- rnorm(nsample,0,sigma_e)
#   U2 <- runif(nsample,0,1)
#   
#   Y1star <- Y1 + e
#   Y2star <- ifelse(U2>expit(alpha),Y2,1-Y2)
#   
#   ## Naive model
#   naive.model1 <- lm(Y1star ~ X + W)
#   true.model1 <- lm(Y1 ~ X + W)
#   naive.model2 <- glm(Y2star ~ X + W, family = binomial(link = logit))
#   true.model2 <- glm(Y2star ~ X + W, family = binomial(link = logit))
#   
#   ## 1.2 Implementation of the methods #####
#   
#   ## 1.2.1 Preparation ###
#   
#   DesignMatrix1 <-  DesignMatrix2 <- cbind(rep(1,length(X)),X,W)
#   CovMis1 <- cbind(rep(0,length(X)),rep(0,length(X)))
#   CovMis2 <- c(rep(1,length(X)))
#   
#   DesignMatrix1 <- as.matrix(DesignMatrix1)
#   DesignMatrix2 <- as.matrix(DesignMatrix2)
#   CovMis1 <- as.matrix(CovMis1)
#   CovMis2 <- as.matrix (CovMis2)
#   
#   ## 1.2.2 Prepare different choices of initial variables ###
# 
#   beta_Y1_0 <- mean(Y1star)
#   beta_Y2_0 <- log(mean(Y2star)/(1-mean(Y2star)))
#   
#   intial3 <- c(beta_Y1_0,0,0,beta_Y2_0,0,0,1,0.001)
#   intial4 <- c(naive.model1$coefficients,naive.model2$coefficients,1,0)
#   
#   ### < uncomment the following part to study a broader list ###
#   # 1.2.5 The proposed method ###
#   NR <- nleqslv(intial4, GEE_UI,Y1star=Y1star, Y2star=Y2star, DesignMatrix1=DesignMatrix1, DesignMatrix2=DesignMatrix2,
#                       CovMis1=CovMis1, CovMis2=CovMis2, jacobian=T, control=list(maxit=2000),
#                       gamma1 = 1, gamma=c(0,0), alpha1= alpha, alpha0= alpha, sigma_e = sigma_e)
#   
#   betahat <- ifelse(abs(NR$x)<10,NR$x,NA)
#   
#   # GEE_SIGMA(betahat, Y1star=Y1star, Y2star=Y2star, DesignMatrix1=DesignMatrix1, DesignMatrix2=DesignMatrix2, 
#   #           CovMis1=CovMis1, CovMis2=CovMis2, gamma1 = 1, gamma=c(0,0), alpha1= alpha, alpha0= alpha, sigma_e = sigma_e)
#   
#   if (!any(is.na(betahat))) {
#     cov <- GEE_cov(betahat,Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
#                    gamma1 = 1, gamma=c(0,0), alpha1= alpha, alpha0= alpha, sigma_e = sigma_e)
#     sd <- sqrt(diag(cov))} else {
#     sd <- rep(NA,length(betahat))
#   }
#   
#   # 1.2.3 Naive Model of only consider the measurement error ###
#   measonly <- nleqslv(intial4,GEE_UI,Y1star=Y1star, Y2star=Y2star, DesignMatrix1=DesignMatrix1, DesignMatrix2=DesignMatrix2,
#                       CovMis1=CovMis1, CovMis2=CovMis2, jacobian=T, control=list(maxit=2000),
#                       gamma1 = 1, gamma=c(0,0), alpha1= -Inf, alpha0= -Inf, sigma_e = sigma_e)
#   
#   betahat_measonly <- ifelse(abs(measonly$x)<10,measonly$x,NA)
#   
#   if (!any(is.na(betahat_measonly))) {
#     cov <- GEE_cov(betahat_measonly,Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
#                    gamma1 = 1, gamma=c(0,0), alpha1= -Inf, alpha0= -Inf, sigma_e = sigma_e)
#     sd_measonly  <- sqrt(diag(cov))} else {
#       sd_measonly  <- rep(NA,length(betahat_measonly))
#     }
#   
#   # 1.2.4 Naive Model of only consider the misclassification error ###
#   misconly <- nleqslv(intial4, GEE_UI,Y1star=Y1star, Y2star=Y2star, DesignMatrix1=DesignMatrix1, DesignMatrix2=DesignMatrix2,
#                       CovMis1=CovMis1, CovMis2=CovMis2, jacobian=T, control=list(maxit=2000),
#                       gamma1 = 1, gamma=c(0,0), alpha1= alpha, alpha0= alpha, sigma_e = 0)
#   
#   betahat_misconly <- ifelse(abs(misconly$x)<10,misconly$x,NA)
#   
#   if (!any(is.na(betahat_misconly))) {
#     cov <- GEE_cov(betahat_misconly,Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
#                    gamma1 = 1, gamma=c(0,0), alpha1= -Inf, alpha0= -Inf, sigma_e = sigma_e)
#     sd_misconly <- sqrt(diag(cov))} else {
#       sd_misconly <- rep(NA,length(betahat_misconly))
#     }
#   
#   return(list(seed = seed_i[i],
#               naive1coef = naive.model1$coefficients,
#               naive1vcov = vcov(naive.model1),
#               naive2coef = naive.model2$coefficients,
#               naive2vcov = vcov(naive.model2),
#               betameasonly = betahat_measonly,
#               sdmeasonly = sd_measonly,
#               betamisconly = betahat_misconly,
#               sdmisconly = sd_misconly,
#               betahat = betahat,
#               sd = sd))
# }

# alpha <- c(-3,-0.5)   #5%
# alpha <- c(-2,-0.5)   #10%
# alpha <- c(-1.5,-0.5)   #20%
##### 2 Implementation Function: The correlated mismeasruement ####
gamma <- 0.5

INS_int <- function(i, sigma_e, alpha, gamma){
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
  
  # theta <- c(beta1, beta2, sigma, xi)
  sigma <- 1
  rho <- 0
  # gamma <- 0.8
  
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
  e <- rnorm(nsample,0,sigma_e)
  U2 <- runif(nsample,0,1)
  
  Z <- X
  CovMis2 <- data.frame(intercept= rep(1,nsample), Z=Z)
  CovMis2 <- as.matrix(CovMis2)
  
  alphastcov <- as.matrix(CovMis2) %*% t(t(alpha))
  
  Y1star <- Y1 + gamma * Y2 + e
  Y2star <- ifelse(U2 > expit(alphastcov),Y2,1-Y2)
  
  mean(U2 > expit(alphastcov))
  ## Naive model
  naive.model1 <- lm(Y1star ~ X + W)
  true.model1 <- lm(Y1 ~ X + W)
  naive.model2 <- glm(Y2star ~ X + W, family = binomial(link = logit))
  true.model2 <- glm(Y2 ~ X + W, family = binomial(link = logit))
  
  
  ## 2.3 Implementation Generation ###
  ## 2.3.1 Preperation ###
  
  DesignMatrix1 <-  DesignMatrix2 <- cbind(rep(1,length(X)),X,W)
  CovMis1 <- cbind(rep(0,length(X)),rep(1,length(X)))
  CovMis2 <- cbind(rep(1,length(X)),Z)
  
  DesignMatrix1 <- as.matrix(DesignMatrix1)
  DesignMatrix2 <- as.matrix(DesignMatrix2)
  CovMis1 <- as.matrix(CovMis1)
  # CovMis2 <- as.matrix (CovMis2)
  
  ## 2.3.2 Prepare different choices of initial variables ###

  beta_Y1_0 <- mean(Y1star)
  beta_Y2_0 <- log(mean(Y2star)/(1-mean(Y2star)))

  intial3 <- c(beta_Y1_0,0,0,beta_Y2_0,0,0,1,0.001)
  intial4 <- c(naive.model1$coefficients,naive.model2$coefficients,1,0)
  
  tryCatch({
  # 2.3.5 The proposed method ###
  NR <- nleqslv(intial4, GEE_UI,Y1star=Y1star, Y2star=Y2star, DesignMatrix1=DesignMatrix1, DesignMatrix2=DesignMatrix2,
                CovMis1=CovMis1, CovMis2=CovMis2, jacobian=T, control=list(maxit=2000),
                gamma1 = 1, gamma=c(0,gamma), alpha1= alpha, alpha0= alpha, sigma_e = sigma_e)
  
  betahat <- ifelse(abs(NR$x)<10,NR$x,NA)
  
  if (!any(is.na(betahat))) {
    cov <- GEE_cov(betahat,Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                   gamma1 = 1, gamma=c(0,gamma), alpha1= alpha, alpha0= alpha, sigma_e = sigma_e)
    sd <- sqrt(diag(cov))} else {
      sd <- rep(NA,length(betahat))
    }

  # 2.3.3 Naive Model of only consider the measurement error ###
  measonly <- nleqslv(intial4,GEE_UI,Y1star=Y1star, Y2star=Y2star, DesignMatrix1=DesignMatrix1, DesignMatrix2=DesignMatrix2,
                      CovMis1=CovMis1, CovMis2=CovMis2, jacobian=T, control=list(maxit=2000),
                      gamma1 = 1, gamma=c(0,gamma), alpha1= c(-Inf,0), alpha0= c(-Inf,0), sigma_e = sigma_e)

  betahat_measonly <- ifelse(abs(measonly$x)<10,measonly$x,NA)

  if (!any(is.na(betahat_measonly))) {
    cov <- GEE_cov(betahat_measonly,Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                   gamma1 = 1, gamma=c(0,gamma), alpha1= c(-Inf,0), alpha0= c(-Inf,0), sigma_e = sigma_e)
    sd_measonly  <- sqrt(diag(cov))} else {
      sd_measonly  <- rep(NA,length(betahat_measonly))
    }

  # 2.3.4 Naive Model of only consider the misclassification error ###
  misconly <- nleqslv(intial4, GEE_UI,Y1star=Y1star, Y2star=Y2star, DesignMatrix1=DesignMatrix1, DesignMatrix2=DesignMatrix2,
                      CovMis1=CovMis1, CovMis2=CovMis2, jacobian=T, control=list(maxit=2000),
                      gamma1 = 1, gamma=c(0,0), alpha1= alpha, alpha0= alpha, sigma_e = 0)

  betahat_misconly <- ifelse(abs(misconly$x)<10,misconly$x,NA)

  if (!any(is.na(betahat_misconly))) {
    cov <- GEE_cov(betahat_misconly,Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                   gamma1 = 1, gamma=c(0,0), alpha1= alpha, alpha0= alpha, sigma_e = 0)
    sd_misconly <- sqrt(diag(cov))} else {
      sd_misconly <- rep(NA,length(betahat_misconly))
    }
  
  
  
  return(list(seed = seed_i[i],
              naive1coef = naive.model1$coefficients,
              naive1vcov = vcov(naive.model1),
              naive2coef = naive.model2$coefficients,
              naive2vcov = vcov(naive.model2),
              betameasonly = betahat_measonly,
              sdmeasonly = sd_measonly,
              betamisconly = betahat_misconly,
              sdmisconly = sd_misconly,
              betahat = betahat,
              sd = sd))
  
  }, error = function(e) return(NULL))
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
# results_1 <- lapply(c(0.1,0.5,0.7), FUN= function(x){
#   results_x <- lapply(1:1000, FUN = INS_int,
#                          sigma_e = x, alpha = -2.19722, gamma = 0.8)
#   return(results_x)
# })
# 
# save(results_1, file="KNOWN_R1.RData")
# 
# truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0)
# sigma_e_range <- c(0.1,0.5,0.7)
# 
# 
# 
# Results <- NULL
# 
# for (k in 1:3) {
#   results <- results_1[[k]]
#   naive1coef <- NULL
#   naive1sd <- NULL
#   naive2coef <- NULL
#   naive2sd <- NULL
#   CI1naive <- NULL
#   CI2naive <- NULL
#   measonlycoef <- NULL
#   measonlysd <- NULL
#   CImeasonly <- NULL
#   misconlycoef <- NULL
#   misconlysd <- NULL
#   CImisconly <- NULL
# 
#   betas <- NULL
#   sds <- NULL
#   CIs <- NULL
#   for (i in 1:1000){
#     if (is.null(results[[i]])) {
#         next}
#     naive1coef <- rbind(naive1coef, results[[i]]$naive1coef)
#     naive1sd <- rbind(naive1sd,sqrt(diag( results[[i]]$naive1vcov)))
#     naive2coef <- rbind(naive2coef, results[[i]]$naive2coef)
#     naive2sd <- rbind(naive2sd,sqrt(diag( results[[i]]$naive2vcov)))
# 
#     if ((!is.null(results[[i]]$betameasonly)) & (!is.null(results[[i]]$sdmeasonly))) {
#       measonlycoef <- rbind(measonlycoef,as.vector(results[[i]]$betameasonly))
#       measonlysd_i <- results[[i]]$sdmeasonly
#       measonlysd_i <- ifelse(abs(measonlysd_i)<100,measonlysd_i,NA)
#       measonlysd <- rbind(measonlysd,measonlysd_i)
# 
#       CILBmeasonly <- results[[i]]$betameasonly - 1.96 *(measonlysd_i)
#       CIUBmeasonly <- results[[i]]$betameasonly + 1.96 *(measonlysd_i)
#       CImeasonly <- rbind(CImeasonly,ifelse((truebeta<as.vector(CIUBmeasonly)) & (truebeta>as.vector(CILBmeasonly)),1,0))
#     }
# 
#     if ((!is.null(results[[i]]$betamisconly)) & (!is.null(results[[i]]$sdmisconly))) {
#       misconlycoef <- rbind(misconlycoef,as.vector(results[[i]]$betamisconly))
#       misconlysd_i <- (results[[i]]$sdmisconly)
#       misconlysd_i <- ifelse(abs(misconlysd_i)<100,misconlysd_i,NA)
#       misconlysd <- rbind(misconlysd, misconlysd_i)
# 
#       CILBmisconly <- results[[i]]$betamisconly - 1.96 *(misconlysd_i)
#       CIUBmisconly <- results[[i]]$betamisconly + 1.96 *(misconlysd_i)
#       CImisconly <- rbind(CImisconly,ifelse((truebeta<as.vector(CIUBmisconly)) & (truebeta>as.vector(CILBmisconly)),1,0))
#     }
# 
#     betahat0 <- results[[i]]$betahat
#     sd0 <- results[[i]]$sd
#     # sd0 <- ifelse(abs(sd0)<100,sd0,NA)
#     betas <- rbind(betas, betahat0)
#     sds <- rbind(sds, sd0)
# 
#     CILBnaive1 <- results[[i]]$naive1coef - 1.96 *(sqrt(diag(results[[i]]$naive1vcov)))
#     CIUBnaive1 <- results[[i]]$naive1coef + 1.96 *(sqrt(diag(results[[i]]$naive1vcov)))
#     CI1naive <- rbind(CI1naive,ifelse((truebeta[1:3]<CIUBnaive1) & (truebeta[1:3]>CILBnaive1),1,0))
#     CILBnaive2 <- results[[i]]$naive2coef - 1.96 *(sqrt(diag(results[[i]]$naive2vcov)))
#     CIUBnaive2 <- results[[i]]$naive2coef + 1.96 *(sqrt(diag(results[[i]]$naive2vcov)))
#     CI2naive <- rbind(CI2naive,ifelse((truebeta[4:6]<CIUBnaive2) & (truebeta[4:6]>CILBnaive2),1,0))
#     CILB <- betahat0 - 1.96 *(sd0)
#     CIUB <- betahat0 + 1.96 *(sd0)
#     CIs <- rbind(CIs,ifelse((truebeta<as.vector(CIUB)) & (truebeta>as.vector(CILB)),1,0))
#   }
# 
#   biasnaive1 <- colMeans(naive1coef)-truebeta[1:3]
#   biasnaive2 <- colMeans(naive2coef)-truebeta[4:6]
#   naive_esd <- apply(cbind(naive1coef,naive2coef), MARGIN = 2 , FUN=sd, na.rm=T)
#   sdnaive1 <- colMeans(naive1sd)
#   sdnaive2 <- colMeans(naive2sd)
#   CInaive1 <- colMeans(CI1naive,na.rm=T)
#   CInaive2 <- colMeans(CI2naive,na.rm=T)
# 
#   naivebias <- c(biasnaive1,biasnaive2,0,0)
#   naive_esd <- c(naive_esd,0,0)
#   naivesd <- c(sdnaive1,sdnaive2,0,0)
#   naiveCI <- c(CInaive1,CInaive2,0,0)
# 
#   bias_measonly <- colMeans(na.omit(measonlycoef),na.rm = T) - truebeta
#   sd_emp_measonly <- apply(na.omit(measonlycoef),MARGIN = 2, FUN = function(x){
#     x.noout <-  remove_outliers(x)
#     return( sd(x.noout,na.rm = T ) )
#   })
#   sd_mod_measonly <- apply(na.omit(measonlysd),MARGIN = 2, FUN = function(x){
#     x.noout <-  remove_outliers(x)
#     return( mean(x.noout,na.rm = T ))
#   })
#   CI_measonly <- colMeans(na.omit(CImeasonly),na.rm = T)
# 
#   bias_misconly <- colMeans(na.omit(misconlycoef),na.rm = T) - truebeta
#   sd_emp_misconly <- apply(na.omit(misconlycoef),MARGIN = 2, FUN = function(x){
#     x.noout <-  remove_outliers(x)
#     return( sd(x.noout,na.rm = T ) )
#   })
#   sd_mod_misconly <- apply(na.omit(misconlysd),MARGIN = 2, FUN = function(x){
#     x.noout <-  remove_outliers(x)
#     return( mean(x.noout,na.rm = T ))
#   })
#   CI_misconly <- colMeans(na.omit(CImisconly),na.rm = T)
# 
#   bias1 <- colMeans(na.omit(betas),na.rm = T) - truebeta
#   sd_emp <- apply(na.omit(betas),MARGIN = 2, FUN = function(x){
#    x.noout <-  remove_outliers(x)
#    return( sd(x.noout,na.rm = T ) )
#   })
#   
#   sd_mod <- sd_mod <- apply(na.omit(sds),MARGIN = 2, FUN = function(x){
#     x.noout <-  remove_outliers(x)
#     return( mean(x.noout,na.rm = T ))
#   })
#   CIrate <- colMeans(na.omit(CIs),na.rm = T)
# 
#   Results0 <- data.frame(sigma_e = sigma_e_range[k],
#                 naivebias=round(naivebias,3),naive_esd=round(naive_esd,3),naivesd=round(naivesd,3),naiveCI=percent(round(naiveCI,3)),
#                 biasmeas=round(bias_measonly,3),esdmeas=round(sd_emp_measonly,3),sdmeas=round(sd_mod_measonly,3),CImeas=percent(round(CI_measonly,3)),
#                 biasmisc=round(bias_misconly,3),esdmisc=round(sd_emp_misconly,3),sdmisc=round(sd_mod_misconly,3),CImisc=percent(round(CI_misconly,3)),
#                 biasprop=round(bias1,3),propose_esd=round(sd_emp,3),sdpropose=round(sd_mod,3),CI_propose=percent(round(CIrate,3)))
# 
#   Results <- rbind(Results,Results0)
# }
# 
# 
# 
# library(xtable)
# xtable(Results,digits = 3)

### 4.2 Simulation 2: under different degree of misclassification rates ####

# alpha <- c(-3,-0.5)   #5%
# alpha <- c(-2,-0.5)   #10%
# alpha <- c(-1.5,-0.5)   #20%

alphatable <- matrix(c(-3,-0.5,-2,-0.5,-1.5,-0.5),nrow=2)


results_2 <- lapply(1:3, FUN= function(x){
  results_x <- lapply(1:1000, FUN = INS_int,
                      sigma_e = 0.1, alpha = alphatable[,x], gamma = 0.8)
  return(results_x)
})

save(results_2, file="KNOWN_R2_2.RData")
truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0)
alpha_range <- c(-4.595,-2.19722,-1.38629)
# alpha_range <- c(-2.197,-2.2,-2.197225)
Results <- NULL

for (k in 1:3) {
  results <- results_2[[k]]
  naive1coef <- NULL
  naive1sd <- NULL
  naive2coef <- NULL
  naive2sd <- NULL
  CI1naive <- NULL
  CI2naive <- NULL
  measonlycoef <- NULL
  measonlysd <- NULL
  CImeasonly <- NULL
  misconlycoef <- NULL
  misconlysd <- NULL
  CImisconly <- NULL

  betas <- NULL
  sds <- NULL
  CIs <- NULL
  for (i in 1:1000){
    if (is.null(results[[i]])) {
        next}
    naive1coef <- rbind(naive1coef, results[[i]]$naive1coef)
    naive1sd <- rbind(naive1sd,sqrt(diag( results[[i]]$naive1vcov)))
    naive2coef <- rbind(naive2coef, results[[i]]$naive2coef)
    naive2sd <- rbind(naive2sd,sqrt(diag( results[[i]]$naive2vcov)))

    if ((!is.null(results[[i]]$betameasonly)) & (!is.null(results[[i]]$sdmeasonly))) {
      measonlycoef <- rbind(measonlycoef,as.vector(results[[i]]$betameasonly))
      measonlysd_i <- results[[i]]$sdmeasonly
      measonlysd_i <- ifelse(abs(measonlysd_i)<100,measonlysd_i,NA)
      measonlysd <- rbind(measonlysd,measonlysd_i)

      CILBmeasonly <- results[[i]]$betameasonly - 1.96 *(measonlysd_i)
      CIUBmeasonly <- results[[i]]$betameasonly + 1.96 *(measonlysd_i)
      CImeasonly <- rbind(CImeasonly,ifelse((truebeta<as.vector(CIUBmeasonly)) & (truebeta>as.vector(CILBmeasonly)),1,0))
    }

    if ((!is.null(results[[i]]$betamisconly)) & (!is.null(results[[i]]$sdmisconly))) {
      misconlycoef <- rbind(misconlycoef,as.vector(results[[i]]$betamisconly))
      misconlysd_i <- (results[[i]]$sdmisconly)
      misconlysd_i <- ifelse(abs(misconlysd_i)<100,misconlysd_i,NA)
      misconlysd <- rbind(misconlysd, misconlysd_i)

      CILBmisconly <- results[[i]]$betamisconly - 1.96 *(misconlysd_i)
      CIUBmisconly <- results[[i]]$betamisconly + 1.96 *(misconlysd_i)
      CImisconly <- rbind(CImisconly,ifelse((truebeta<as.vector(CIUBmisconly)) & (truebeta>as.vector(CILBmisconly)),1,0))
    }

    betahat0 <- results[[i]]$betahat
    sd0 <- results[[i]]$sd
    # sd0 <- ifelse(abs(sd0)<100,sd0,NA)
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

  bias_measonly <- colMeans(na.omit(measonlycoef),na.rm = T) - truebeta
  sd_emp_measonly <- apply(na.omit(measonlycoef),MARGIN = 2, FUN = function(x){
    x.noout <-  remove_outliers(x)
    return( sd(x.noout,na.rm = T ))
  })
  
  sd_mod_measonly <- apply(na.omit(measonlysd),MARGIN = 2, FUN = function(x){
    x.noout <-  remove_outliers(x)
    return( mean(x.noout,na.rm = T ))
  })
  CI_measonly <- colMeans(na.omit(CImeasonly),na.rm = T)

  bias_misconly <- colMeans(na.omit(misconlycoef),na.rm = T) - truebeta
  sd_emp_misconly <- apply(na.omit(misconlycoef),MARGIN = 2, FUN = function(x){
    x.noout <-  remove_outliers(x)
    return( sd(x.noout,na.rm = T ))
  })
  sd_mod_misconly <- apply(na.omit(misconlysd),MARGIN = 2, FUN = function(x){
    x.noout <-  remove_outliers(x)
    return( mean(x.noout,na.rm = T ))
  })
  CI_misconly <- colMeans(na.omit(CImisconly),na.rm = T)

  
  
  bias1 <- apply(na.omit(betas),MARGIN = 2, FUN = function(x){
    x.noout <-  remove_outliers(x)
    return( mean(x.noout,na.rm = T ))
  }) - truebeta
  
  sd_emp <- apply(na.omit(betas),MARGIN = 2, FUN = function(x){
    x.noout <-  remove_outliers(x)
    return( sd(x.noout,na.rm = T ))
  })
  
  sd_mod <- apply(na.omit(sds),MARGIN = 2, FUN = function(x){
    x.noout <-  remove_outliers(x)
    return( mean(x.noout,na.rm = T ))
  })
  
  CIrate <- colMeans(na.omit(CIs),na.rm = T)

  Results0 <- data.frame(alpha0 = alphatable[1,k], alpha1 = alphatable[2,k],
                         naivebias=round(naivebias,3),naive_esd=round(naive_esd,3),naivesd=round(naivesd,3),naiveCI=percent(naiveCI,accuracy=.1),
                         biasmeas=round(bias_measonly,3),esdmeas=round(sd_emp_measonly,3),sdmeas=round(sd_mod_measonly,3),CImeas=percent(CI_measonly,accuracy=.1),
                         biasmisc=round(bias_misconly,3),esdmisc=round(sd_emp_misconly,3),sdmisc=round(sd_mod_misconly,3),CImisc=percent(CI_misconly,accuracy=.1),
                         biasprop=round(bias1,3),propose_esd=round(sd_emp,3),sdpropose=round(sd_mod,3),CI_propose=percent(CIrate,accuracy=.1))

  Results <- rbind(Results,Results0)
}

betas0new <- remove_outliers(betas[,1])

library(xtable)
xtable(Results,digits = 3)

### 4.3 Simulation 3: under different level of interation ####
# results_3 <- lapply(c(0,-0.8,0.8), FUN= function(x){
#   results_x <- lapply(1:1000, FUN = INS_int,
#                       sigma_e = 0.1, alpha = -2.19722, gamma = x)
#   return(results_x)
# })
# 
# save(results_3, file="KNOWN_R3.RData")
# 
# truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0)
# gamma_range <- c(0,-0.8,0.8)
# 
# Results <- NULL
# 
# for (k in 1:3) {
#   results <- results_3[[k]]
#   naive1coef <- NULL
#   naive1sd <- NULL
#   naive2coef <- NULL
#   naive2sd <- NULL
#   CI1naive <- NULL
#   CI2naive <- NULL
#   measonlycoef <- NULL
#   measonlysd <- NULL
#   CImeasonly <- NULL
#   misconlycoef <- NULL
#   misconlysd <- NULL
#   CImisconly <- NULL
# 
#   betas <- NULL
#   sds <- NULL
#   CIs <- NULL
#   for (i in 1:1000){
#     naive1coef <- rbind(naive1coef, results[[i]]$naive1coef)
#     naive1sd <- rbind(naive1sd,sqrt(diag( results[[i]]$naive1vcov)))
#     naive2coef <- rbind(naive2coef, results[[i]]$naive2coef)
#     naive2sd <- rbind(naive2sd,sqrt(diag( results[[i]]$naive2vcov)))
# 
#     if ((!is.null(results[[i]]$betameasonly)) & (!is.null(results[[i]]$sdmeasonly))) {
#       measonlycoef <- rbind(measonlycoef,as.vector(results[[i]]$betameasonly))
#       measonlysd_i <- results[[i]]$sdmeasonly
#       measonlysd_i <- ifelse(abs(measonlysd_i)<100,measonlysd_i,NA)
#       measonlysd <- rbind(measonlysd,measonlysd_i)
# 
#       CILBmeasonly <- results[[i]]$betameasonly - 1.96 *(measonlysd_i)
#       CIUBmeasonly <- results[[i]]$betameasonly + 1.96 *(measonlysd_i)
#       CImeasonly <- rbind(CImeasonly,ifelse((truebeta<as.vector(CIUBmeasonly)) & (truebeta>as.vector(CILBmeasonly)),1,0))
#     }
# 
#     if ((!is.null(results[[i]]$betamisconly)) & (!is.null(results[[i]]$sdmisconly))) {
#       misconlycoef <- rbind(misconlycoef,as.vector(results[[i]]$betamisconly))
#       misconlysd_i <- (results[[i]]$sdmisconly)
#       misconlysd_i <- ifelse(abs(misconlysd_i)<100,misconlysd_i,NA)
#       misconlysd <- rbind(misconlysd, misconlysd_i)
# 
#       CILBmisconly <- results[[i]]$betamisconly - 1.96 *(misconlysd_i)
#       CIUBmisconly <- results[[i]]$betamisconly + 1.96 *(misconlysd_i)
#       CImisconly <- rbind(CImisconly,ifelse((truebeta<as.vector(CIUBmisconly)) & (truebeta>as.vector(CILBmisconly)),1,0))
#     }
# 
#     betahat0 <- results[[i]]$betahat
#     sd0 <- results[[i]]$sd
#     # sd0 <- ifelse(abs(sd0)<100,sd0,NA)
#     betas <- rbind(betas, betahat0)
#     sds <- rbind(sds, sd0)
# 
#     CILBnaive1 <- results[[i]]$naive1coef - 1.96 *(sqrt(diag(results[[i]]$naive1vcov)))
#     CIUBnaive1 <- results[[i]]$naive1coef + 1.96 *(sqrt(diag(results[[i]]$naive1vcov)))
#     CI1naive <- rbind(CI1naive,ifelse((truebeta[1:3]<CIUBnaive1) & (truebeta[1:3]>CILBnaive1),1,0))
#     CILBnaive2 <- results[[i]]$naive2coef - 1.96 *(sqrt(diag(results[[i]]$naive2vcov)))
#     CIUBnaive2 <- results[[i]]$naive2coef + 1.96 *(sqrt(diag(results[[i]]$naive2vcov)))
#     CI2naive <- rbind(CI2naive,ifelse((truebeta[4:6]<CIUBnaive2) & (truebeta[4:6]>CILBnaive2),1,0))
#     CILB <- betahat0 - 1.96 *(sd0)
#     CIUB <- betahat0 + 1.96 *(sd0)
#     CIs <- rbind(CIs,ifelse((truebeta<as.vector(CIUB)) & (truebeta>as.vector(CILB)),1,0))
#   }
# 
#   biasnaive1 <- colMeans(naive1coef)-truebeta[1:3]
#   biasnaive2 <- colMeans(naive2coef)-truebeta[4:6]
#   naive_esd <- apply(cbind(naive1coef,naive2coef), MARGIN = 2 , FUN=sd, na.rm=T)
#   sdnaive1 <- colMeans(naive1sd)
#   sdnaive2 <- colMeans(naive2sd)
#   CInaive1 <- colMeans(CI1naive,na.rm=T)
#   CInaive2 <- colMeans(CI2naive,na.rm=T)
# 
#   naivebias <- c(biasnaive1,biasnaive2,0,0)
#   naive_esd <- c(naive_esd,0,0)
#   naivesd <- c(sdnaive1,sdnaive2,0,0)
#   naiveCI <- c(CInaive1,CInaive2,0,0)
# 
#   bias_measonly <- colMeans(na.omit(measonlycoef),na.rm = T) - truebeta
#   sd_emp_measonly <- apply(na.omit(measonlycoef),MARGIN = 2, FUN = function(x){
#     x.noout <-  remove_outliers(x)
#     return( sd(x.noout,na.rm = T ))
#   })
#   sd_mod_measonly <- apply(na.omit(measonlysd),MARGIN = 2, FUN = function(x){
#     x.noout <-  remove_outliers(x)
#     return( mean(x.noout,na.rm = T ))
#   })
#   CI_measonly <- colMeans(na.omit(CImeasonly),na.rm = T)
# 
#   bias_misconly <- colMeans(na.omit(misconlycoef),na.rm = T) - truebeta
#   sd_emp_misconly <- apply(na.omit(misconlycoef),MARGIN = 2, FUN = function(x){
#     x.noout <-  remove_outliers(x)
#     return( sd(x.noout,na.rm = T ))
#   })
#   sd_mod_misconly <- apply(na.omit(misconlysd),MARGIN = 2, FUN = function(x){
#     x.noout <-  remove_outliers(x)
#     return( mean(x.noout,na.rm = T ))
#   })
#   CI_misconly <- colMeans(na.omit(CImisconly),na.rm = T)
# 
#   bias1 <- colMeans(na.omit(betas),na.rm = T) - truebeta
#   sd_emp <- apply(na.omit(betas),MARGIN = 2, FUN = function(x){
#     x.noout <-  remove_outliers(x)
#     return( sd(x.noout,na.rm = T ))
#   })
#   sd_mod <- apply(na.omit(sds),MARGIN = 2, FUN = function(x){
#     x.noout <-  remove_outliers(x)
#     return( mean(x.noout,na.rm = T ))
#   })
#   CIrate <- colMeans(na.omit(CIs),na.rm = T)
# 
#   Results0 <- data.frame( gamma = gamma_range[k],
#     naivebias=round(naivebias,3),naive_esd=round(naive_esd,3),naivesd=round(naivesd,3),naiveCI=percent(round(naiveCI,3)),
#     biasmeas=round(bias_measonly,3),esdmeas=round(sd_emp_measonly,3),sdmeas=round(sd_mod_measonly,3),CImeas=percent(round(CI_measonly,3)),
#     biasmisc=round(bias_misconly,3),esdmisc=round(sd_emp_misconly,3),sdmisc=round(sd_mod_misconly,3),CImisc=percent(round(CI_misconly,3)),
#     biasprop=round(bias1,3),propose_esd=round(sd_emp,3),sdpropose=round(sd_mod,3),CI_propose=percent(round(CIrate,3)))
# 
#   Results <- rbind(Results,Results0)
# }
# 
# 
# 
# library(xtable)
# xtable(Results,digits = 3)