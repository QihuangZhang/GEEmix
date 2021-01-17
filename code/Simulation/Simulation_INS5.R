####
# This simulation is testing whether the first version of Insertion method has unbiased property.


#### 0 Simulation set up####
library(parallel)
library(scales)

## 0.1 Generate the seed for the simulation ####

set.seed(2018)
seed_i <- sample(1000000,1000)

## 0.2 Global Parameters ####

ncores <- 30
ProjectName <- "Simulation_INS5"
nsample <- 1000

## 0.3 Functions  ####
## 0.3.1 Original Functions  ####
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
  GAMMA.inv <- solve(GAMMA,tol=1e-200)
  return(GAMMA.inv)
}

GEE_cov <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                    gamma1, gamma, alpha1, alpha0, sigma_e){
  GAMMA <- GEE_GAMMA(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
  GAMMA.inv <- GEE_GAMMA.inv(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
  SIGMA <- GEE_SIGMA(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                     gamma1, gamma, alpha1, alpha0, sigma_e)
  covmatrix <- GAMMA.inv %*% SIGMA %*% t(as.matrix(GAMMA.inv))
   return(covmatrix)
  # return(list(GAMMA=GAMMA,SIGMA=SIGMA,covmatrix))
}


## 0.3.1 Functions with Internal Validation  ####

GEE_UI_IV <- function(theta,  data.validation, data.mismeasure,
                    gamma1, gamma, alpha1, alpha0, sigma_e, Weight){
  # cat(theta, " \n")
  return(GEE_UfuncInsIVWeight(Y1star=data.mismeasure$Y1star,
                        Y2star=data.mismeasure$Y2star,
                        Y1 = data.validation$Y1,
                        Y2 = data.validation$Y2,
                        DesignMatrix1 = as.matrix(data.mismeasure[,3:5]), 
                        DesignMatrix2 = as.matrix(data.mismeasure[,3:5]), 
                        ValidationMatrix1 = as.matrix(data.validation[,5:7]),
                        ValidationMatrix2 = as.matrix(data.validation[,5:7]),
                        CovMis1 = as.matrix(data.mismeasure[,6:7]),
                        CovMis2 = as.matrix(data.mismeasure[,8]),
                        Weight = Weight,
                        beta1=theta[1:3], beta2=theta[4:6], sigma = theta[7], xi = theta[8], 
                        gamma1=gamma1, gamma=gamma, alpha1=alpha1, alpha0=alpha0, sigma_e=sigma_e
                      ))
}



GEE_GAMMA_IV0 <- function(theta, Y1star, Y2star, Y1, Y2, ValidationMatrix1, ValidationMatrix2, CovMis1, CovMis2, 
                          gamma1, gamma, alpha1, alpha0, sigma_e,
                          fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0){
  return(GEE_GAMMAInsIV0(Y1star, Y2star, Y1, Y2,
                         CovMis1, CovMis2, ValidationMatrix1, ValidationMatrix2,
                         beta1=theta[1:3], beta2=theta[4:6], xi=theta[8], sigma=theta[7],
                         gamma1, gamma, alpha1, alpha0, sigma_e, 
                         fixgamma1=fixgamma1, fixgamma=fixgamma, fixsigma_e=fixsigma_e, fixalpha1=fixalpha1, fixalpha0=fixalpha0)
  )
}


GEE_GAMMA_IVI <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, 
                          gamma1, gamma, alpha1, alpha0, sigma_e,
                          fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0){
  return(GEE_GAMMAInsIVI(Y1star, Y2star, DesignMatrix1, DesignMatrix2,
                         CovMis1, CovMis2,
                         beta1=theta[1:3], beta2=theta[4:6], xi=theta[8], sigma=theta[7],
                         gamma1, gamma, alpha1, alpha0, sigma_e,
                         fixgamma1=fixgamma1, fixgamma=fixgamma, fixsigma_e=fixsigma_e, fixalpha1=fixalpha1, fixalpha0=fixalpha0)
         )
}

GEE_SIGMA_IV0 <- function(theta, Y1star, Y2star, Y1, Y2, ValidationMatrix1, ValidationMatrix2, CovMis1, CovMis2, 
                          gamma1, gamma, alpha1, alpha0, sigma_e,
                          fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0){
  return(GEE_GAMMAInsIV0(Y1star, Y2star, Y1, Y2,
                         CovMis1, CovMis2, ValidationMatrix1, ValidationMatrix2,
                         beta1=theta[1:3], beta2=theta[4:6], xi=theta[8], sigma=theta[7],
                         gamma1, gamma, alpha1, alpha0, sigma_e, 
                         fixgamma1=fixgamma1, fixgamma=fixgamma, fixsigma_e=fixsigma_e, fixalpha1=fixalpha1, fixalpha0=fixalpha0)
  )
}

GEE_SIGMA_IV0 <- function(theta, Y1star, Y2star, Y1, Y2, ValidationMatrix1, ValidationMatrix2, CovMis1, CovMis2, 
                          gamma1, gamma, alpha1, alpha0, sigma_e,
                          fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0){
  return(GEE_SIGMAInsIV0(Y1star, Y2star, Y1, Y2, ValidationMatrix1, ValidationMatrix2,
                         CovMis1, CovMis2,
                         beta1=theta[1:3], beta2=theta[4:6], xi=theta[8], sigma=theta[7],
                         gamma1, gamma, alpha1, alpha0, sigma_e,
                         fixgamma1=fixgamma1, fixgamma=fixgamma, fixsigma_e=fixsigma_e, fixalpha1=fixalpha1, fixalpha0=fixalpha0)
  )
}


GEE_SIGMA_IVI <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, 
                          gamma1, gamma, alpha1, alpha0, sigma_e,
                          fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0){
  return(GEE_SIGMAInsIVI(Y1star, Y2star, DesignMatrix1, DesignMatrix2,
                         CovMis1, CovMis2,
                         beta1=theta[1:3], beta2=theta[4:6], xi=theta[8], sigma=theta[7],
                         gamma1, gamma, alpha1, alpha0, sigma_e,
                         fixgamma1=fixgamma1, fixgamma=fixgamma, fixsigma_e=fixsigma_e, fixalpha1=fixalpha1, fixalpha0=fixalpha0)
  )
}

GEE_covIV <- function(theta, data.validation, data.mismeasure, Weight,
                    gamma1, gamma, alpha1, alpha0, sigma_e,
                    fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0){
  nvalidation <- dim(data.validation)[1]
  nsample <- dim(data.mismeasure)[1] + nvalidation
  Omega <- diag(Weight)
  I_Omega <- diag(1-Weight)
  
  M0 <- GEE_GAMMA_IV0(theta, 
                      Y1star=data.validation$Y1star, 
                      Y2star=data.validation$Y2star, 
                      Y1 = data.validation$Y1, 
                      Y2 = data.validation$Y2, 
                      ValidationMatrix1 = as.matrix(data.validation[,5:7]),
                      ValidationMatrix2 = as.matrix(data.validation[,5:7]), 
                      CovMis1 = as.matrix(data.validation[,8:9]), 
                      CovMis2 = as.matrix(data.validation[,10]), 
                      gamma1, gamma, alpha1, alpha0, sigma_e,
                      fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0)
  
  M1 <- GEE_GAMMA_IVI(theta, 
                      Y1star=data.mismeasure$Y1star, 
                      Y2star=data.mismeasure$Y2star, 
                      DesignMatrix1 = as.matrix(data.mismeasure[,3:5]),
                      DesignMatrix2 = as.matrix(data.mismeasure[,3:5]), 
                      CovMis1 = as.matrix(data.mismeasure[,6:7]), 
                      CovMis2 = as.matrix(data.mismeasure[,8]), 
                      gamma1, gamma, alpha1, alpha0, sigma_e,
                      fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0)
  
  GAMMA_IV <- Omega %*% M1 + I_Omega %*% M0
  
  B0 <- GEE_SIGMA_IV0(theta, 
                      Y1star=data.validation$Y1star, 
                      Y2star=data.validation$Y2star, 
                      Y1 = data.validation$Y1, 
                      Y2 = data.validation$Y2, 
                      ValidationMatrix1 = as.matrix(data.validation[,5:7]),
                      ValidationMatrix2 = as.matrix(data.validation[,5:7]), 
                      CovMis1 = as.matrix(data.validation[,8:9]), 
                      CovMis2 = as.matrix(data.validation[,10]), 
                      gamma1, gamma, alpha1, alpha0, sigma_e,
                      fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0)
  
  B1 <- GEE_SIGMA_IVI(theta, 
                      Y1star=data.mismeasure$Y1star, 
                      Y2star=data.mismeasure$Y2star, 
                      DesignMatrix1 = as.matrix(data.mismeasure[,3:5]),
                      DesignMatrix2 = as.matrix(data.mismeasure[,3:5]), 
                      CovMis1 = as.matrix(data.mismeasure[,6:7]), 
                      CovMis2 = as.matrix(data.mismeasure[,8]), 
                      gamma1, gamma, alpha1, alpha0, sigma_e,
                      fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0)
  
  SIGMA_IV <- Omega %*% B1 %*% Omega + I_Omega %*% B0 %*% I_Omega 
  
  GAMMA.inv <- solve(GAMMA_IV,tol=1e-200)
  covmatrix <- GAMMA.inv %*% SIGMA_IV %*% t(as.matrix(GAMMA.inv))
  return(covmatrix)
  # return(list(M1=M1,M0=M0,B1=B1,B0=B0,GAMMA_IV=GAMMA_IV,SIGMA_IV=SIGMA_IV,covmatrix))
}





##### 1 Implementation Function ####

### A example in deguging
i <- 1
nsample <- 1500
nvalidation <- 500
omega_j <- 1

INS_int <- function(i, nsample, nvalidation, omega_j){
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
  
  omegas <- rep(omega_j,6+2)
  
  # theta <- c(beta1, beta2, sigma, xi)
  sigma <- 1
  rho <- 0
  
  sigma_e <- 0.1 
  gamma <- 0.8
  alpha <- -2.197
  # gamma <- 0.8
  
  ### Generate the true data sets
  # nsample <- 1000
  # nvalidation <- 500
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
  
  data.mismeasure <- data.frame(Y1star=Y1star[1:(nsample - nvalidation)],Y2star=Y2star[1:(nsample - nvalidation)], DesignMatrix1[1:(nsample - nvalidation),],CovMis1[1:(nsample - nvalidation),],CovMis2[1:(nsample - nvalidation),]) 
  data.validation <- data.frame(Y1=Y1[(nsample - nvalidation+1):nsample],Y2=Y2[(nsample - nvalidation+1):nsample], Y1star=Y1star[(nsample - nvalidation+1):nsample],Y2star=Y2star[(nsample - nvalidation+1):nsample],
                                DesignMatrix1[(nsample - nvalidation+1):nsample,],CovMis1[(nsample - nvalidation+1):nsample,],CovMis2[(nsample - nvalidation+1):nsample,]) 

  ## 1.3.2 Prepare different choices of initial variables ###

  beta_Y1_0 <- mean(Y1star)
  beta_Y2_0 <- log(mean(Y2star)/(1-mean(Y2star)))

  intial3 <- c(beta_Y1_0,0,0,beta_Y2_0,0,0,1,0.001)
  intial4 <- c(naive.model1$coefficients,naive.model2$coefficients,1,0)
  
  
  ## 1.4 Estimating Procedure
  
  ## 1.4.1 Measurement Error and Misclassification Parameters
  model.measure <- lm(Y1star ~ -1 + offset(Y1) + Y2,data = data.validation) 
  model.class1 <- glm((1-Y2star) ~ 1, data = data.validation[data.validation$Y2==1,],family = binomial(link="logit")) 
  model.class0 <- glm(Y2star ~ 1, data = data.validation[data.validation$Y2==0,],family = binomial(link="logit")) 
   
  gamma2 <- model.measure$coefficients
  sigma_e <- sigma(model.measure)
  alpha1 <- model.class1$coefficients
  alpha0 <- model.class0$coefficients
  
  tryCatch({
  # 1.4.2 The proposed method ####
  NR <- nleqslv(intial4, GEE_UI_IV, jacobian=T, control=list(maxit=2000),
                data.validation = data.validation, data.mismeasure = data.mismeasure, 
                Weight = c(omegas,rep(0,4)),
                gamma1 = 1, gamma=c(0,gamma2), alpha1= alpha1, alpha0= alpha0, sigma_e = sigma_e)
  
  betahat <- ifelse(abs(NR$x)<10,NR$x,NA)
  
  ### variance estimation with validation data
  if (!any(is.na(betahat))) {
    cov <- GEE_covIV (betahat, data.validation, data.mismeasure, Weight = c(omegas,rep(0,4)),
                                 gamma1=1, gamma = c(0,gamma2), alpha1, alpha0, sigma_e,
                      fixgamma1=1, fixgamma=c(1,0), fixsigma_e=0, fixalpha1=0, fixalpha0=0)
    sd <- sqrt(diag(cov))} else {
      sd <- rep(NA,length(betahat))
    }

  # 2.3.3 Naive Model of only consider the measurement error ###
  measonly <- nleqslv(intial4,GEE_UI_IV,jacobian=T, control=list(maxit=2000),
                      data.validation = data.validation, data.mismeasure = data.mismeasure, 
                      Weight = c(omegas,rep(0,4)),
                      gamma1 = 1, gamma=c(0,gamma2), alpha1= -Inf, alpha0= -Inf, sigma_e = sigma_e)

  betahat_measonly <- ifelse(abs(measonly$x)<10,measonly$x,NA)

  if (!any(is.na(betahat_measonly))) {
    cov <- GEE_covIV(betahat_measonly,  data.validation, data.mismeasure, Weight = c(omegas,rep(0,2)),
                   gamma1 = 1, gamma = c(0,gamma2), alpha1= -Inf, alpha0= -Inf, sigma_e = sigma_e,
                   fixgamma1=1, fixgamma=c(1,0), fixsigma_e=0, fixalpha1=1, fixalpha0=1)
    sd_measonly  <- sqrt(diag(cov))} else {
      sd_measonly  <- rep(NA,length(betahat_measonly))
    }

  # 2.3.4 Naive Model of only consider the misclassification error ###
  misconly <- nleqslv(intial4,GEE_UI_IV,jacobian=T, control=list(maxit=2000),
                      Weight = c(omegas,rep(0,4)),
                      data.validation = data.validation, data.mismeasure = data.mismeasure,
                      gamma1 = 1, gamma=c(0,0), alpha1= alpha1, alpha0= alpha0, sigma_e = 0)

  betahat_misconly <- ifelse(abs(misconly$x)<10,misconly$x,NA)

  if (!any(is.na(betahat_misconly))) {
    cov <- GEE_covIV(betahat_measonly,  data.validation, data.mismeasure, Weight = c(omegas,rep(0,2)),
                   gamma1 = 1, gamma=c(0,0), alpha1= alpha1, alpha0= alpha0, sigma_e = 0,
                   fixgamma1=1, fixgamma=c(1,1), fixsigma_e=1, fixalpha1=0, fixalpha0=0)
    sd_misconly <- sqrt(diag(cov))} else {
      sd_misconly <- rep(NA,length(betahat_misconly))
    }
  
  
  
  return(list(seed = seed_i[i],
              naive1coef = naive.model1$coefficients,
              naive1vcov = vcov(naive.model1),
              naive2coef = naive.model2$coefficients,
              naive2vcov = vcov(naive.model2),
              betameasonly = c(betahat_measonly,gamma2,sigma_e,0,0),
              sdmeasonly = c(sd_measonly,0,0),
              betamisconly = c(betahat_misconly,0,0,alpha1,alpha0),
              sdmisconly = c(sd_misconly[1:8],0,0,sd_misconly[9:10]),
              betahat = c(betahat,gamma2,sigma_e,alpha1,alpha0),
              sd = sd))
  
  }, error = function(e) return(NULL))
}


# ### 4.1 Simulation 1: under different degree of measurement error ####
# results_1 <- lapply(c(750,1500,3000), FUN= function(x){
#   results_x <- lapply(1:1000, FUN = INS_int,
#                       nsample = x, nvalidation = x*1/3, omega_j = 0.5)
#   return(results_x)
# })
# 
# # re1 <- INS_int(1,nsample=nsample, nvalidation=nvalidation, omega_j= omega_j)
# 
# nsample <- c(750,1500,3000)
# 
# # truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0,1,0.1,0.1,-1.39,-1.39)
# 
# save(results_1,file="WINV_R1.RData")
# 
# Results <- NULL
# 
# for (k in 1:3) {
#   results <- results_1[[k]]
#   truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0,0.8,0.1,-2.197,-2.197)
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
#       measonlysd_i <- ifelse(abs(measonlysd_i)<10,measonlysd_i,NA)
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
#       misconlysd_i <- ifelse(abs(misconlysd_i)<10,misconlysd_i,NA)
#       misconlysd <- rbind(misconlysd, misconlysd_i)
# 
#       CILBmisconly <- results[[i]]$betamisconly - 1.96 *(misconlysd_i)
#       CIUBmisconly <- results[[i]]$betamisconly + 1.96 *(misconlysd_i)
#       CImisconly <- rbind(CImisconly,ifelse((truebeta<as.vector(CIUBmisconly)) & (truebeta>as.vector(CILBmisconly)),1,0))
#     }
# 
#     betahat0 <- results[[i]]$betahat
#     sd0 <- results[[i]]$sd
#     sd0 <- ifelse(abs(sd0)<10,sd0,NA)
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
#   biasnaive1 <- colMeans(naive1coef,na.rm=T)-truebeta[1:3]
#   biasnaive2 <- colMeans(naive2coef,na.rm=T)-truebeta[4:6]
#   naive_esd <- apply(cbind(naive1coef,naive2coef), MARGIN = 2 , FUN=sd, na.rm=T)
#   sdnaive1 <- colMeans(naive1sd,na.rm=T)
#   sdnaive2 <- colMeans(naive2sd,na.rm=T)
#   CInaive1 <- colMeans(CI1naive,na.rm=T)
#   CInaive2 <- colMeans(CI2naive,na.rm=T)
# 
#   naivebias <- c(biasnaive1,biasnaive2,rep(0,6))
#   naive_esd <- c(naive_esd,rep(0,6))
#   naivesd <- c(sdnaive1,sdnaive2,rep(0,6))
#   naiveCI <- c(CInaive1,CInaive2,rep(0,6))
# 
#   bias_measonly <- colMeans(na.omit(measonlycoef),na.rm = T) - truebeta
#   sd_emp_measonly <- apply(na.omit(measonlycoef),MARGIN = 2, FUN = sd)
#   sd_mod_measonly <- colMeans(na.omit(measonlysd),na.rm = T)
#   CI_measonly <- colMeans(na.omit(CImeasonly),na.rm = T)
# 
#   bias_misconly <- colMeans(na.omit(misconlycoef),na.rm = T) - truebeta
#   sd_emp_misconly <- apply(na.omit(misconlycoef),MARGIN = 2, FUN = sd)
#   sd_mod_misconly <- colMeans(na.omit(misconlysd),na.rm = T)
#   CI_misconly <- colMeans(na.omit(CImisconly),na.rm = T)
# 
#   bias1 <- colMeans(na.omit(betas),na.rm = T) - truebeta
#   sd_emp <- apply(na.omit(betas),MARGIN = 2, FUN = sd)
#   sd_mod <- colMeans(na.omit(sds),na.rm = T)
#   CIrate <- colMeans(na.omit(CIs),na.rm = T)
# 
#   Results0 <- data.frame(nsample = nsample[k],
#                 naivebias=round(naivebias,3),naive_esd=round(naive_esd,3),naivesd=round(naivesd,3),naiveCI=percent(round(naiveCI,3)),
#                 biasmeas=round(bias_measonly,3),esdmeas=round(sd_emp_measonly,3),sdmeas=round(sd_mod_measonly,3),CImeas=percent(round(CI_measonly,3)),
#                 biasmisc=round(bias_misconly,3),esdmisc=round(sd_emp_misconly,3),sdmisc=round(sd_mod_misconly,3),CImisc=percent(round(CI_misconly,3)),
#                 biasprop=round(bias1,3),propose_esd=round(sd_emp,3),sdpropose=round(sd_mod,3),CI_propose=percent(round(CIrate,3)))
# 
#   Results <- rbind(Results,Results0)
# }
# 
# save(Results,file="WINV_RTable1.RData")
# 
# library(xtable)
# xtable(Results,digits = 3)

### 4.2 Simulation 2: under different degree of misclassification rates ####
results_2 <- lapply(c(1/3,1/2,2/3), FUN= function(x){
  results_x <- lapply(1:1000, FUN = INS_int,
                      nsample = 1500, nvalidation = round(x*1500), omega_j = 0.5)
  return(results_x)
})

# re1 <- INS_int(1,nsample=nsample, nvalidation=nvalidation, omega_j= omega_j)

nrate <- round(c(1/3,1/2,2/3),3)


truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0)
alpha_range <- c(-4.595,-2.197,-1.386)

save(results_2,file="WINV_R2.RData")


  Results <- NULL
  
  
for (k in 1:3) {
  results <- results_2[[k]]
  truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0,0.8,0.1,-2.197,-2.197)
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
      measonlysd_i <- ifelse(abs(measonlysd_i)<10,measonlysd_i,NA)
      measonlysd <- rbind(measonlysd,measonlysd_i)
      
      CILBmeasonly <- results[[i]]$betameasonly - 1.96 *(measonlysd_i)
      CIUBmeasonly <- results[[i]]$betameasonly + 1.96 *(measonlysd_i)
      CImeasonly <- rbind(CImeasonly,ifelse((truebeta<as.vector(CIUBmeasonly)) & (truebeta>as.vector(CILBmeasonly)),1,0))
    } 
    
    if ((!is.null(results[[i]]$betamisconly)) & (!is.null(results[[i]]$sdmisconly))) {
      misconlycoef <- rbind(misconlycoef,as.vector(results[[i]]$betamisconly))
      misconlysd_i <- (results[[i]]$sdmisconly)
      misconlysd_i <- ifelse(abs(misconlysd_i)<10,misconlysd_i,NA)
      misconlysd <- rbind(misconlysd, misconlysd_i)
      
      CILBmisconly <- results[[i]]$betamisconly - 1.96 *(misconlysd_i)
      CIUBmisconly <- results[[i]]$betamisconly + 1.96 *(misconlysd_i)
      CImisconly <- rbind(CImisconly,ifelse((truebeta<as.vector(CIUBmisconly)) & (truebeta>as.vector(CILBmisconly)),1,0))
    }
    
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
  
  biasnaive1 <- colMeans(naive1coef,na.rm=T)-truebeta[1:3]
  biasnaive2 <- colMeans(naive2coef,na.rm=T)-truebeta[4:6]
  naive_esd <- apply(cbind(naive1coef,naive2coef), MARGIN = 2 , FUN=sd, na.rm=T)
  sdnaive1 <- colMeans(naive1sd,na.rm=T)
  sdnaive2 <- colMeans(naive2sd,na.rm=T)
  CInaive1 <- colMeans(CI1naive,na.rm=T) 
  CInaive2 <- colMeans(CI2naive,na.rm=T) 
  
  naivebias <- c(biasnaive1,biasnaive2,rep(0,6))
  naive_esd <- c(naive_esd,rep(0,6))
  naivesd <- c(sdnaive1,sdnaive2,rep(0,6))
  naiveCI <- c(CInaive1,CInaive2,rep(0,6))
  
  bias_measonly <- colMeans(na.omit(measonlycoef),na.rm = T) - truebeta
  sd_emp_measonly <- apply(na.omit(measonlycoef),MARGIN = 2, FUN = function(x){
    x.noout <-  remove_outliers(x)
    return( sd(x.noout,na.rm = T ))
  })
  sd_mod_measonly <- colMeans(na.omit(measonlysd),na.rm = T)
  CI_measonly <- colMeans(na.omit(CImeasonly),na.rm = T)
  
  bias_misconly <- colMeans(na.omit(misconlycoef),na.rm = T) - truebeta
  sd_emp_misconly <- apply(na.omit(misconlycoef),MARGIN = 2, FUN = function(x){
    x.noout <-  remove_outliers(x)
    return( sd(x.noout,na.rm = T ))
  })
  sd_mod_misconly <- colMeans(na.omit(misconlysd),na.rm = T)
  CI_misconly <- colMeans(na.omit(CImisconly),na.rm = T)
  
  bias1 <- colMeans(na.omit(betas),na.rm = T) - truebeta
  sd_emp <- apply(na.omit(betas),MARGIN = 2, FUN = function(x){
    x.noout <-  remove_outliers(x)
    return( sd(x.noout,na.rm = T ))
  })
  sd_mod <- colMeans(na.omit(sds),na.rm = T)
  CIrate <- colMeans(na.omit(CIs),na.rm = T)
  
  Results0 <- data.frame(nrate = nrate[k], 
                         naivebias=round(naivebias,3),naive_esd=round(naive_esd,3),naivesd=round(naivesd,3),naiveCI=percent(round(naiveCI,3)),
                         biasmeas=round(bias_measonly,3),esdmeas=round(sd_emp_measonly,3),sdmeas=round(sd_mod_measonly,3),CImeas=percent(round(CI_measonly,3)),
                         biasmisc=round(bias_misconly,3),esdmisc=round(sd_emp_misconly,3),sdmisc=round(sd_mod_misconly,3),CImisc=percent(round(CI_misconly,3)),
                         biasprop=round(bias1,3),propose_esd=round(sd_emp,3),sdpropose=round(sd_mod,3),CI_propose=percent(round(CIrate,3)))
  
  Results <- rbind(Results,Results0)
}


save(Results,file="WINV_RTable2small.RData")
library(xtable)
xtable(Results,digits = 3)

### 4.2 Simulation 2: under different degree of misclassification rates ####
results_2 <- lapply(c(1/3,1/2,2/3), FUN= function(x){
  results_x <- lapply(1:1000, FUN = INS_int,
                      nsample = 3000, nvalidation = round(x*3000), omega_j = 0.5)
  return(results_x)
})

# re1 <- INS_int(1,nsample=nsample, nvalidation=nvalidation, omega_j= omega_j)

nrate <- round(c(1/3,1/2,2/3),3)


truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0)
alpha_range <- c(-4.595,-2.197,-1.386)

save(results_2,file="WINV_R22.RData")


Results <- NULL


for (k in 1:3) {
  results <- results_2[[k]]
  truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0,0.8,0.1,-2.197,-2.197)
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
      measonlysd_i <- ifelse(abs(measonlysd_i)<10,measonlysd_i,NA)
      measonlysd <- rbind(measonlysd,measonlysd_i)
      
      CILBmeasonly <- results[[i]]$betameasonly - 1.96 *(measonlysd_i)
      CIUBmeasonly <- results[[i]]$betameasonly + 1.96 *(measonlysd_i)
      CImeasonly <- rbind(CImeasonly,ifelse((truebeta<as.vector(CIUBmeasonly)) & (truebeta>as.vector(CILBmeasonly)),1,0))
    } 
    
    if ((!is.null(results[[i]]$betamisconly)) & (!is.null(results[[i]]$sdmisconly))) {
      misconlycoef <- rbind(misconlycoef,as.vector(results[[i]]$betamisconly))
      misconlysd_i <- (results[[i]]$sdmisconly)
      misconlysd_i <- ifelse(abs(misconlysd_i)<10,misconlysd_i,NA)
      misconlysd <- rbind(misconlysd, misconlysd_i)
      
      CILBmisconly <- results[[i]]$betamisconly - 1.96 *(misconlysd_i)
      CIUBmisconly <- results[[i]]$betamisconly + 1.96 *(misconlysd_i)
      CImisconly <- rbind(CImisconly,ifelse((truebeta<as.vector(CIUBmisconly)) & (truebeta>as.vector(CILBmisconly)),1,0))
    }
    
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
  
  biasnaive1 <- colMeans(naive1coef,na.rm=T)-truebeta[1:3]
  biasnaive2 <- colMeans(naive2coef,na.rm=T)-truebeta[4:6]
  naive_esd <- apply(cbind(naive1coef,naive2coef), MARGIN = 2 , FUN=sd, na.rm=T)
  sdnaive1 <- colMeans(naive1sd,na.rm=T)
  sdnaive2 <- colMeans(naive2sd,na.rm=T)
  CInaive1 <- colMeans(CI1naive,na.rm=T) 
  CInaive2 <- colMeans(CI2naive,na.rm=T) 
  
  naivebias <- c(biasnaive1,biasnaive2,rep(0,6))
  naive_esd <- c(naive_esd,rep(0,6))
  naivesd <- c(sdnaive1,sdnaive2,rep(0,6))
  naiveCI <- c(CInaive1,CInaive2,rep(0,6))
  
  bias_measonly <- colMeans(na.omit(measonlycoef),na.rm = T) - truebeta
  sd_emp_measonly <- apply(na.omit(measonlycoef),MARGIN = 2, FUN = function(x){
    x.noout <-  remove_outliers(x)
    return( sd(x.noout,na.rm = T ))
  })
  sd_mod_measonly <- colMeans(na.omit(measonlysd),na.rm = T)
  CI_measonly <- colMeans(na.omit(CImeasonly),na.rm = T)
  
  bias_misconly <- colMeans(na.omit(misconlycoef),na.rm = T) - truebeta
  sd_emp_misconly <- apply(na.omit(misconlycoef),MARGIN = 2, FUN = function(x){
    x.noout <-  remove_outliers(x)
    return( sd(x.noout,na.rm = T ))
  })
  sd_mod_misconly <- colMeans(na.omit(misconlysd),na.rm = T)
  CI_misconly <- colMeans(na.omit(CImisconly),na.rm = T)
  
  bias1 <- colMeans(na.omit(betas),na.rm = T) - truebeta
  sd_emp <- apply(na.omit(betas),MARGIN = 2, FUN = function(x){
    x.noout <-  remove_outliers(x)
    return( sd(x.noout,na.rm = T ))
  })
  sd_mod <- colMeans(na.omit(sds),na.rm = T)
  CIrate <- colMeans(na.omit(CIs),na.rm = T)
  
  Results0 <- data.frame(nrate = nrate[k], 
                         naivebias=round(naivebias,3),naive_esd=round(naive_esd,3),naivesd=round(naivesd,3),naiveCI=percent(round(naiveCI,3)),
                         biasmeas=round(bias_measonly,3),esdmeas=round(sd_emp_measonly,3),sdmeas=round(sd_mod_measonly,3),CImeas=percent(round(CI_measonly,3)),
                         biasmisc=round(bias_misconly,3),esdmisc=round(sd_emp_misconly,3),sdmisc=round(sd_mod_misconly,3),CImisc=percent(round(CI_misconly,3)),
                         biasprop=round(bias1,3),propose_esd=round(sd_emp,3),sdpropose=round(sd_mod,3),CI_propose=percent(round(CIrate,3)))
  
  Results <- rbind(Results,Results0)
}


save(Results,file="WINV_RTable2large.RData")
library(xtable)
xtable(Results,digits = 3)

### 4.3 Simulation 3: under different level of interation ####
# results_3 <- lapply(c(1/3,0.5,2/3), FUN= function(x){
#   results_x <- lapply(1:1000, FUN = INS_int,
#                       nsample = 1500, nvalidation = 500, omega_j = x)
#   return(results_x)
# })
# 
# # re1 <- INS_int(1,nsample=nsample, nvalidation=nvalidation, omega_j= omega_j)
# 
# omega <- round(c(1/3,1/2,2/3),3)
# 
# 
# truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0)
# gamma_range <- c(0,-0.8,0.8)
# 
# save(results_3,file="WINV_R3.RData")
# 
# Results <- NULL
# 
# 
# for (k in 1:3) {
#   results <- results_3[[k]]
#   truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0,0.8,0.1,-2.197,-2.197)
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
#   
#     
#   for (i in 1:1000){
#     if (is.null(results[[i]])) {
#       next}
#     naive1coef <- rbind(naive1coef, results[[i]]$naive1coef)
#     naive1sd <- rbind(naive1sd,sqrt(diag( results[[i]]$naive1vcov)))
#     naive2coef <- rbind(naive2coef, results[[i]]$naive2coef)
#     naive2sd <- rbind(naive2sd,sqrt(diag( results[[i]]$naive2vcov)))
#     
#     if ((!is.null(results[[i]]$betameasonly)) & (!is.null(results[[i]]$sdmeasonly))) {
#       measonlycoef <- rbind(measonlycoef,as.vector(results[[i]]$betameasonly))
#       measonlysd_i <- results[[i]]$sdmeasonly
#       measonlysd_i <- ifelse(abs(measonlysd_i)<10,measonlysd_i,NA)
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
#       misconlysd_i <- ifelse(abs(misconlysd_i)<10,misconlysd_i,NA)
#       misconlysd <- rbind(misconlysd, misconlysd_i)
#       
#       CILBmisconly <- results[[i]]$betamisconly - 1.96 *(misconlysd_i)
#       CIUBmisconly <- results[[i]]$betamisconly + 1.96 *(misconlysd_i)
#       CImisconly <- rbind(CImisconly,ifelse((truebeta<as.vector(CIUBmisconly)) & (truebeta>as.vector(CILBmisconly)),1,0))
#     }
#     
#     betahat0 <- results[[i]]$betahat
#     sd0 <- results[[i]]$sd
#     sd0 <- ifelse(abs(sd0)<5,sd0,NA)
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
#   biasnaive1 <- colMeans(naive1coef,na.rm=T)-truebeta[1:3]
#   biasnaive2 <- colMeans(naive2coef,na.rm=T)-truebeta[4:6]
#   naive_esd <- apply(cbind(naive1coef,naive2coef), MARGIN = 2 , FUN=sd, na.rm=T)
#   sdnaive1 <- colMeans(naive1sd,na.rm=T)
#   sdnaive2 <- colMeans(naive2sd,na.rm=T)
#   CInaive1 <- colMeans(CI1naive,na.rm=T) 
#   CInaive2 <- colMeans(CI2naive,na.rm=T) 
#   
#   naivebias <- c(biasnaive1,biasnaive2,rep(0,6))
#   naive_esd <- c(naive_esd,rep(0,6))
#   naivesd <- c(sdnaive1,sdnaive2,rep(0,6))
#   naiveCI <- c(CInaive1,CInaive2,rep(0,6))
#   
#   bias_measonly <- colMeans(na.omit(measonlycoef),na.rm = T) - truebeta
#   sd_emp_measonly <- apply(na.omit(measonlycoef),MARGIN = 2, FUN = sd)
#   sd_mod_measonly <- colMeans(na.omit(measonlysd),na.rm = T)
#   CI_measonly <- colMeans(na.omit(CImeasonly),na.rm = T)
#   
#   bias_misconly <- colMeans(na.omit(misconlycoef),na.rm = T) - truebeta
#   sd_emp_misconly <- apply(na.omit(misconlycoef),MARGIN = 2, FUN = sd)
#   sd_mod_misconly <- colMeans(na.omit(misconlysd),na.rm = T)
#   CI_misconly <- colMeans(na.omit(CImisconly),na.rm = T)
#   
#   bias1 <- colMeans(na.omit(betas),na.rm = T) - truebeta
#   sd_emp <- apply(na.omit(betas),MARGIN = 2, FUN = sd)
#   sd_mod <- colMeans(na.omit(sds),na.rm = T)
#   CIrate <- colMeans(na.omit(CIs),na.rm = T)
#   
#   Results0 <- data.frame( omega = omega[k], 
#     naivebias=round(naivebias,3),naive_esd=round(naive_esd,3),naivesd=round(naivesd,3),naiveCI=percent(round(naiveCI,3)),
#     biasmeas=round(bias_measonly,3),esdmeas=round(sd_emp_measonly,3),sdmeas=round(sd_mod_measonly,3),CImeas=percent(round(CI_measonly,3)),
#     biasmisc=round(bias_misconly,3),esdmisc=round(sd_emp_misconly,3),sdmisc=round(sd_mod_misconly,3),CImisc=percent(round(CI_misconly,3)),
#     biasprop=round(bias1,3),propose_esd=round(sd_emp,3),sdpropose=round(sd_mod,3),CI_propose=percent(round(CIrate,3)))
#   
#   Results <- rbind(Results,Results0)
# }
# 
# 
# save(Results,file="WINV_RTable3.RData")
# library(xtable)
# xtable(Results,digits = 3)