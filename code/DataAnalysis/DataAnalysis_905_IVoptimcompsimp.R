
### Project GEEmix

### This is the project with internal validation data

#### 0 Project Information and Configuration ####
ProjectName<- paste0("DataAnalysis_109")
cat(paste0("Note: This is for the data analysis of SNP 109 \n")) 
set.seed(2019)

## 0.1 Set file path ####
# setwd("N:/Work/waterloo/2018/GEEmix")
# WD=paste0(getwd(),"/results/EMsensi/")


## Finished part -- no need to run again
# genotype<-read.table(file="data/geno.txt",header=T)
# save(genotype,file="genotype.RData")
# load("~/Work/waterloo/2018/mixResponse/file/genotype.RData")
# F_PC <- prcomp(genotype[,3:dim(genotype)[2]],center = T,scale.=T)

# load("N:/Work/waterloo/2018/GEEmix/file/F_PC.RData")
# loading <- F_PC$x
# save(loading,file="loading.RData")
# screeplot(F_PC)  # obtain the screeplot


# obtain the screeplot
# Scree <- data.frame(variance = (F_PC$sdev)^2, percentage = cumsum((F_PC$sdev)^2/sum((F_PC$sdev)^2)))
# Scree <- Scree[1:10,]
# rownamesa <- NULL
# for (i in 1:10){
#   rownamesa <- c(rownamesa, paste0("PC",i))
# }
# 
# Scree <- data.frame(rownamesa = rownamesa,Scree)
# Scree$rownamesa <-  factor(Scree$rownamesa,levels = Scree$rownamesa)
# library(ggplot2)
# 
# ggplot(Scree)  + 
#   geom_bar(aes(x=rownamesa, y=variance),stat="identity")+
#   geom_line(aes(x=rownamesa, y=percentage*max(Scree$variance)*10,group = 1),stat="identity")+
#   labs(x = "Principal Components", y = "Variance") +
#   scale_y_continuous(sec.axis = sec_axis(~./max(Scree$variance)/10, name ="Cumulative percentage"))

load("~/Work/waterloo/2018/mixResponse/file/genotrunc.RData")
phenotype<-read.csv(file="~/Work/waterloo/2018/mixResponse/data/pheno.csv")
load("~/Work/waterloo/2018/GEEmix/file/loading.RData")
load("~/Work/waterloo/2018/GEEmix/file/rs27338905.RData")

### Note:!! if use windows, change ~/ back to P:/ or N:/


## 0.2 Prepare package of GWAS toolbox ####
# library(maxLik)
library(MASS)
library(GeneErrorMis)
library(parallel)
library(nleqslv)
# source("code/CGRM.R")

## 0.3 Prepare the datasets ####

## Predict the TAstar
modelTA <- lm(TA~gastroc+EDL+plantaris, data=phenotype)
phenotype$TAstar <- predict(modelTA, newdata=phenotype)

BMD90 <- quantile(phenotype$BMD,na.rm=T, probs = 0.9)
phenotype$BCstar <- ifelse(phenotype$BMD<BMD90,0,1)

## Association Data Frame
data_GWA_Pre <- data.frame(id = phenotype$id,Y1 = phenotype$TA,Y2 = phenotype$BCstar, discard = ifelse(phenotype$discard=="no",1,0),
                           mixup = ifelse(phenotype$mixup=="no",0,1), 
                           batch= as.numeric(phenotype$round %in% c("SW16","SW6","SW17","SW1","SW11","SW3","SW4")),
                           BW = phenotype$bw0, BMD=phenotype$BMD, Y1star=phenotype$TAstar, Y2star=phenotype$abnormalbone)


data_GWAS_unsorted <- data_GWA_Pre[data_GWA_Pre$id %in% genotype$id,]
data_GWAS <- data_GWAS_unsorted[order(data_GWAS_unsorted$id),]

# Response
Y1star <- data_GWAS$Y1star
Y2star <- data_GWAS$Y2star
Y1 <- data_GWAS$Y1
Y2 <- data_GWAS$Y2

## Covariates for main model
Covariates<-cbind(Intercept=rep(1,length(genotype[,2])),X905=rs27338905,data_GWAS[,c("BW","BMD","mixup","batch")],
                  PC1 =loading[,2], PC2 =loading[,3])
Covariates<-as.matrix(Covariates)




## Mark the complete sample
comp<-data.frame(Y1star,Y2star,Covariates[,c("BW","X905","PC1","PC2")])
index.notna<-complete.cases(comp)
index.notna <- ifelse(data_GWAS$discard == 1, index.notna , FALSE)

Y1star<-Y1star[index.notna]
Y2star<-Y2star[index.notna]
Y1<-Y1[index.notna]
Y2<-Y2[index.notna]
Covariates_all<-Covariates[index.notna,]
nsample<-dim(Covariates_all)[1]


CovMis1 <- as.matrix(data.frame(intercept=rep(0,nsample),cov=rep(0,nsample)))
  
CovMis2 <- as.matrix(data.frame(intercept=rep(1,nsample),BW=Covariates_all[,"BW"]))




## 0.4 Function set up ####

## 0.4.1 Functions of Original Version ####
GEE_UI <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, alpha1, alpha0, sigma_e){
  # cat(theta, " \n")
  return(GEE_UfuncIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                      theta[1:4], theta[5:8], sigma = theta[9], xi = theta[10], 
                      gamma1 = 1, gamma=c(0,0),  alpha1, alpha0, sigma_e
  ))
}

GEE_SIGMA <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, gamma, alpha1, alpha0,sigma_e){
  return(GEE_SIGMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                      theta[1:4], theta[5:8], sigma = theta[9], xi = theta[10], 
                      gamma1 = 1, gamma, alpha1, alpha0, sigma_e))
}

GEE_GAMMA <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2){
  GAMMA <- GEE_GAMMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, 
                        beta1=theta[1:4], beta2=theta[5:8], sigma = theta[9], xi = theta[10])
  return(GAMMA)
}

GEE_GAMMA.inv <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2){
  GAMMA <- GEE_GAMMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, 
                        beta1=theta[1:4], beta2=theta[5:8], sigma = theta[9], xi = theta[10])
  GAMMA.inv <- solve(GAMMA,tol=1e-200)
  return(GAMMA.inv)
}

GEE_cov <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, gamma, alpha1, alpha0, sigma_e){
  GAMMA.inv <- GEE_GAMMA.inv(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
  SIGMA <- GEE_SIGMA(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, gamma, alpha1, alpha0, sigma_e)
  covmatrix <- GAMMA.inv %*% SIGMA %*% t(as.matrix(GAMMA.inv))
  return(covmatrix)
}


## 0.4.1 Functions of Original Version ####


GEE_UI_IV <- function(theta,  data.validation, data.mismeasure,
                      gamma1, gamma, alpha1, alpha0, sigma_e, Weight){
  # cat(theta, " \n")
  return(GEE_UfuncInsIVWeight(Y1star=data.mismeasure$Y1star,
                        Y2star=data.mismeasure$Y2star,
                        Y1 = data.validation$Y1,
                        Y2 = data.validation$Y2,
                        DesignMatrix1 = as.matrix(data.mismeasure[,3:6]), 
                        DesignMatrix2 = as.matrix(data.mismeasure[,3:6]), 
                        ValidationMatrix1 = as.matrix(data.validation[,5:8]),
                        ValidationMatrix2 = as.matrix(data.validation[,5:8]),
                        CovMis1 = as.matrix(data.mismeasure[,7:8]),
                        CovMis2 = as.matrix(data.mismeasure[,9]),
                        Weight = Weight,
                        beta1=theta[1:4], beta2=theta[5:8], sigma = theta[9], xi = theta[10], 
                        gamma1=gamma1, gamma=gamma, alpha1=alpha1, alpha0=alpha0, sigma_e=sigma_e
  ))
}


GEE_GAMMA_IV0 <- function(theta, Y1star, Y2star, Y1, Y2, ValidationMatrix1, ValidationMatrix2, CovMis1, CovMis2, 
                          gamma1, gamma, alpha1, alpha0, sigma_e,
                          fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0){
  return(GEE_GAMMAInsIV0(Y1star, Y2star, Y1, Y2,
                         CovMis1, CovMis2, ValidationMatrix1, ValidationMatrix2,
                         beta1=theta[1:4], beta2=theta[5:8], xi=theta[10], sigma=theta[9],
                         gamma1, gamma, alpha1, alpha0, sigma_e, 
                         fixgamma1=fixgamma1, fixgamma=fixgamma, fixsigma_e=fixsigma_e, fixalpha1=fixalpha1, fixalpha0=fixalpha0)
  )
}


GEE_GAMMA_IVI <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, 
                          gamma1, gamma, alpha1, alpha0, sigma_e,
                          fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0){
  return(GEE_GAMMAInsIVI(Y1star, Y2star, DesignMatrix1, DesignMatrix2,
                         CovMis1, CovMis2,
                         beta1=theta[1:4], beta2=theta[5:8], xi=theta[10], sigma=theta[9],
                         gamma1, gamma, alpha1, alpha0, sigma_e,
                         fixgamma1=fixgamma1, fixgamma=fixgamma, fixsigma_e=fixsigma_e, fixalpha1=fixalpha1, fixalpha0=fixalpha0)
  )
}

GEE_SIGMA_IV0 <- function(theta, Y1star, Y2star, Y1, Y2, ValidationMatrix1, ValidationMatrix2, CovMis1, CovMis2, 
                          gamma1, gamma, alpha1, alpha0, sigma_e,
                          fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0){
  return(GEE_GAMMAInsIV0(Y1star, Y2star, Y1, Y2,
                         CovMis1, CovMis2, ValidationMatrix1, ValidationMatrix2,
                         beta1=theta[1:4], beta2=theta[5:8], xi=theta[10], sigma=theta[9],
                         gamma1, gamma, alpha1, alpha0, sigma_e, 
                         fixgamma1=fixgamma1, fixgamma=fixgamma, fixsigma_e=fixsigma_e, fixalpha1=fixalpha1, fixalpha0=fixalpha0)
  )
}

GEE_SIGMA_IV0 <- function(theta, Y1star, Y2star, Y1, Y2, ValidationMatrix1, ValidationMatrix2, CovMis1, CovMis2, 
                          gamma1, gamma, alpha1, alpha0, sigma_e,
                          fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0){
  return(GEE_SIGMAInsIV0(Y1star, Y2star, Y1, Y2, ValidationMatrix1, ValidationMatrix2,
                         CovMis1, CovMis2,
                         beta1=theta[1:4], beta2=theta[5:8], xi=theta[10], sigma=theta[9],
                         gamma1, gamma, alpha1, alpha0, sigma_e,
                         fixgamma1=fixgamma1, fixgamma=fixgamma, fixsigma_e=fixsigma_e, fixalpha1=fixalpha1, fixalpha0=fixalpha0)
  )
}


GEE_SIGMA_IVI <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, 
                          gamma1, gamma, alpha1, alpha0, sigma_e,
                          fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0){
  return(GEE_SIGMAInsIVI(Y1star, Y2star, DesignMatrix1, DesignMatrix2,
                         CovMis1, CovMis2,
                         beta1=theta[1:4], beta2=theta[5:8], xi=theta[10], sigma=theta[9],
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
                      ValidationMatrix1 = as.matrix(data.validation[,5:8]),
                      ValidationMatrix2 = as.matrix(data.validation[,5:8]), 
                      CovMis1 = as.matrix(data.validation[,9:11]), 
                      CovMis2 = as.matrix(data.validation[,12]), 
                      gamma1, gamma, alpha1, alpha0, sigma_e,
                      fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0)
  
  M1 <- GEE_GAMMA_IVI(theta, 
                      Y1star=data.mismeasure$Y1star, 
                      Y2star=data.mismeasure$Y2star, 
                      DesignMatrix1 = as.matrix(data.mismeasure[,3:6]),
                      DesignMatrix2 = as.matrix(data.mismeasure[,3:6]), 
                      CovMis1 = as.matrix(data.mismeasure[,7:9]), 
                      CovMis2 = as.matrix(data.mismeasure[,10]), 
                      gamma1, gamma, alpha1, alpha0, sigma_e,
                      fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0)
  
  GAMMA_IV <- Omega %*% M1 + I_Omega %*% M0
  
  B0 <- GEE_SIGMA_IV0(theta, 
                      Y1star=data.validation$Y1star, 
                      Y2star=data.validation$Y2star, 
                      Y1 = data.validation$Y1, 
                      Y2 = data.validation$Y2, 
                      ValidationMatrix1 = as.matrix(data.validation[,5:8]),
                      ValidationMatrix2 = as.matrix(data.validation[,5:8]), 
                      CovMis1 = as.matrix(data.validation[,9:11]), 
                      CovMis2 = as.matrix(data.validation[,12]), 
                      gamma1, gamma, alpha1, alpha0, sigma_e,
                      fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0)
  
  B1 <- GEE_SIGMA_IVI(theta, 
                      Y1star=data.mismeasure$Y1star, 
                      Y2star=data.mismeasure$Y2star, 
                      DesignMatrix1 = as.matrix(data.mismeasure[,3:6]),
                      DesignMatrix2 = as.matrix(data.mismeasure[,3:6]), 
                      CovMis1 = as.matrix(data.mismeasure[,7:9]), 
                      CovMis2 = as.matrix(data.mismeasure[,10]), 
                      gamma1, gamma, alpha1, alpha0, sigma_e,
                      fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0)
  
  SIGMA_IV <- Omega %*% B1 %*% Omega + I_Omega %*% B0 %*% I_Omega 
  
  GAMMA.inv <- solve(GAMMA_IV,tol=1e-200)
  covmatrix <- GAMMA.inv %*% SIGMA_IV %*% t(as.matrix(GAMMA.inv))
  return(list(M1=M1,M0=M0,B1=B1,B0=B0,GAMMA_IV=GAMMA.inv,SIGMA_IV=SIGMA_IV,covmatrix=covmatrix))
}

Cov01matrix <- function(cov1,cov0,nomegas,nsco){
  Gammaomega1 <- cov0$M0
  Gammaomega1 <- cbind(Gammaomega1,matrix(rep(0,nomegas*(nomegas+nsco)),ncol=nomegas))
  Gammaomega2 <- matrix(rep(0,nomegas*nomegas),ncol=nomegas)
  Gammaomega2 <- cbind(Gammaomega2,cov1$M1[1:nomegas,(nomegas+1):(nomegas+nsco)])
  Gammaomega2 <- cbind(Gammaomega2,cov1$M1[1:nomegas,1:nomegas])
  Gammaomega <- rbind(Gammaomega1,Gammaomega2)
  
  Sigmaomega1 <- cov0$B0
  Sigmaomega1 <- cbind(Sigmaomega1,matrix(rep(0,nomegas*(nomegas+nsco)),ncol=nomegas))
  Sigmaomega2 <- matrix(rep(0,nomegas*(nomegas+nsco)),nrow=nomegas)
  Sigmaomega2 <- cbind(Sigmaomega2,cov1$B1[1:nomegas,1:nomegas])
  Sigmaomega <- rbind(Sigmaomega1,Sigmaomega2)
  
  GAMMAOMG.inv <- solve(Gammaomega,tol=1e-200)
  covmatrix <- GAMMAOMG.inv %*% Sigmaomega %*% t(as.matrix(GAMMAOMG.inv))
  return(list(Var0 = covmatrix[1:(nomegas+nsco),1:(nomegas+nsco)], 
              Var1 = covmatrix[c((nomegas+nsco+1):(2*nomegas+nsco),(nomegas+1):(nomegas+nsco)),
                               c((nomegas+nsco+1):(2*nomegas+nsco),(nomegas+1):(nomegas+nsco))],
              Cov01 = covmatrix[1:(nomegas+nsco),c((nomegas+nsco+1):(2*nomegas+nsco),(nomegas+1):(nomegas+nsco))]))
}

Weightestimate<- function(intial4, omegas, data.validation, data.mismeasure,
                          gamma1, gamma, alpha1, alpha0,sigma_e,
                          fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0,
                          optimal= FALSE){
  nomegas <- length(omegas)
  nsco <- length(c(fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0))-sum(c(fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0))
  NR1 <- nleqslv(intial4, GEE_UI_IV, jacobian=T, control=list(maxit=2000),
                 data.validation = data.validation, data.mismeasure = data.mismeasure, 
                 Weight = c(rep(1,nomegas),rep(0,nsco)),
                 gamma1 = gamma1, gamma = gamma, alpha1= alpha1, alpha0= alpha0, sigma_e = sigma_e)
  
  betahat1 <- ifelse(abs(NR1$x)<100,NR1$x,NA)
  
  ### variance estimation with validation data
  if (!any(is.na(betahat1))) {
    cov1 <- GEE_covIV (betahat1, data.validation, data.mismeasure, Weight = c(rep(1,nomegas),rep(0,nsco)),
                       gamma1 = gamma1, gamma = gamma, alpha1, alpha0, sigma_e,
                       fixgamma1=fixgamma1, fixgamma=fixgamma, fixsigma_e=fixsigma_e, fixalpha1=fixalpha1, fixalpha0=fixalpha0)}
  
  
  NR0 <- nleqslv(intial4, GEE_UI_IV, jacobian=T, control=list(maxit=2000),
                 data.validation = data.validation, data.mismeasure = data.mismeasure, 
                 Weight = c(rep(0,nomegas),rep(0,nsco)),
                 gamma1 = gamma1, gamma=gamma, alpha1= alpha1, alpha0= alpha0, sigma_e = sigma_e)
  
  betahat0 <- ifelse(abs(NR0$x)<100,NR0$x,NA)
  
  ### variance estimation with validation data
  if (!any(is.na(betahat0))) {
    cov0 <- GEE_covIV (betahat0, data.validation, data.mismeasure, Weight = c(rep(0,nomegas),rep(0,nsco)),
                       gamma1=gamma1, gamma = gamma, alpha1=alpha1, alpha0=alpha0, sigma_e=sigma_e,
                       fixgamma1=fixgamma1, fixgamma=fixgamma, fixsigma_e=fixsigma_e, fixalpha1=fixalpha1, fixalpha0=fixalpha0)}
  
  Cov01matrices <- Cov01matrix(cov1,cov0,nomegas,nsco)
  
  if (optimal) {
    weiopttop <- diag(Cov01matrices$Var0) - diag(Cov01matrices$Cov01)
    weioptbot <- diag(Cov01matrices$Var1) + diag(Cov01matrices$Var0) - 2* diag(Cov01matrices$Cov01)
    weiopt <- weiopttop/weioptbot
    weiopt <- ifelse(weiopt<0, 0, weiopt)
    weiopt <- ifelse(weiopt>1, 1, weiopt)
    omegas <- weiopt[1:nomegas]
  }
  
  betahat <- omegas * betahat1 + (1-omegas) * betahat0
  
  Omega <- diag(c(omegas,rep(0,nsco)))
  I_Omega <- diag(1-c(omegas,rep(0,nsco)))
  
  
  Covomega <- Omega %*% Cov01matrices$Var1 %*% Omega + 
    I_Omega %*% Cov01matrices$Var0 %*% I_Omega +
    Omega %*% Cov01matrices$Cov01 %*% I_Omega
  
  return(list(betahat=betahat, Covomega = Covomega))
  
}


#### 1. Data Encryption ####

Covariates1 <- Covariates2 <- as.matrix(cbind(Covariates_all[,c("Intercept","X905","PC1","PC2")])) 

CovMis1 <- as.matrix(data.frame(intercept=rep(1,nsample),cov=rep(1,nsample),BW=Covariates_all[,"BW"]))

CovMis2 <- as.matrix(data.frame(intercept=rep(1,nsample),BW=Covariates_all[,"BW"]))

indexVal <- ifelse(is.na(Covariates_all[,"BMD"])|(Covariates_all[,"mixup"]==1)|(Covariates_all[,"batch"]==1),F,T)

mismeasuredata <- data.frame(Y1star=Y1star[!indexVal], Y2star=Y2star[!indexVal],
                             Covariates1[!indexVal,],CovMis1[!indexVal,],CovMis2[!indexVal,])
validationdata <- data.frame(Y1star=Y1star[indexVal], Y2star=Y2star[indexVal],
                             Y1=Y1[indexVal], Y2=Y2[indexVal],
                             Covariates1[indexVal,],CovMis1[indexVal,],CovMis2[indexVal,])


#### 2. Analysis  of 109 ####


## 2.1 Naive Analysis of 109 ####


beta_Y1_0 <- mean(mismeasuredata$Y1star)
beta_Y2_0 <- log(mean(mismeasuredata$Y2star)/(1-mean(mismeasuredata$Y2star)))



## 2.2 Estimating Procedure ####

# naive.model1 <- lm(Y1star ~ 1)
# naive.model2 <- glm(Y2star ~ 1)

## 2.2.1 Naive model ####
naive.model1 <- lm(Y1star ~ X905 + PC1 + PC2, data = mismeasuredata)
naive.model2 <- glm(Y2star ~ X905 + PC1 + PC2, family = binomial(link = logit), data = mismeasuredata)

intial2 <- c(naive.model1$coefficients,naive.model2$coefficients,0.03,0)
cat("Initial value:", intial2, "\n")


naive.model3 <- lm(Y1 ~ X905 + PC1 + PC2, data = validationdata)
naive.model4 <- glm(Y2 ~ X905 + PC1 + PC2, family = binomial(link = logit), data = validationdata)

intial3 <- c(naive.model3$coefficients,naive.model4$coefficients,1,0)
cat("Initial value:", intial3, "\n")


## 2.2.2 the proposed approach ####

## Measurement Error and Misclassification Parameters
model.measure <- lm(Y1star ~ Y1 + BW,data = validationdata) 
model.class1 <- glm((1-Y2star) ~ 1, data = validationdata[validationdata$Y2==1,],family = binomial(link="logit")) 
model.class0 <- glm(Y2star ~ 1, data = validationdata[validationdata$Y2==0,],family = binomial(link="logit")) 

gamma0 <- model.measure$coefficients[1]
gamma1 <- model.measure$coefficients[2]
gamma3 <- model.measure$coefficients[3] 
sigma_e <- sigma(model.measure)
alpha1 <- model.class1$coefficients
alpha0 <- model.class0$coefficients

## Mainly interested parameters

omegas <- rep(0,length(intial2))
  
NR <- Weightestimate(intial2, omegas = omegas, 
                     data.validation = validationdata, data.mismeasure = mismeasuredata,
                     gamma1=gamma1, gamma=c(gamma0,0,gamma3), alpha1=alpha1, alpha0=alpha0, sigma_e=sigma_e,
                     fixgamma1=0,fixgamma=c(0,1,0),fixsigma_e=0,fixalpha1=0,fixalpha0=0,
                     optimal = T)

betahat <- NR$betahat

### variance estimation with validation data
if (!any(is.na(betahat))) {
  cov <- NR$Covomega
  sd <- sqrt(diag(cov))} else {
    sd <- rep(NA,length(betahat))
  }

betahat <- c(betahat,gamma1,gamma0,gamma3,sigma_e,alpha1,alpha0)

Zvalue <- betahat/sd
pvalue <- 2*(1-pnorm(abs(Zvalue)))

naivebeta <- c(naive.model1$coefficients,naive.model2$coefficients)
sdnaive <- sqrt(c(diag(vcov(naive.model1)),diag(vcov(naive.model1))))

Zvalue_naive <- naivebeta/sdnaive
pvalue_naive <- 2*(1-pnorm(abs(Zvalue_naive)))

Table_numeric <- data.frame(naivebeta=c(naivebeta,rep(0,8)),
                            sdnaive=c(sdnaive,rep(0,8)),
                            Z_naive=c(Zvalue_naive,rep(0,8)),
                            p_naive=c(pvalue_naive,rep(0,8)),
                            propobeta=betahat,
                            proposd = sd,
                            propoZ = Zvalue,
                            propop = pvalue)

library(xtable)
xtable(Table_numeric,digits = 3)

#### 2. Analysis  without considering optimal weights ####

intial2 <- c(betahat[1:8],0.05,0)
nomegas <- length(intial2)

NR1 <- nleqslv(intial2, GEE_UI_IV, jacobian=T, control=list(maxit=2000),
               data.validation = validationdata, data.mismeasure = mismeasuredata, 
               Weight = c(rep(0.5,nomegas),rep(0.5,6)),
               gamma1 = gamma1, gamma = c(gamma0,0,gamma3), alpha1= alpha1, alpha0= alpha0, sigma_e = sigma_e)

betahat2 <- ifelse(abs(NR1$x)<100,NR1$x,NA)

### variance estimation with validation data
if (!any(is.na(betahat))) {
  cov2 <- GEE_covIV (betahat2, data.validation = validationdata, data.mismeasure = mismeasuredata, Weight = c(rep(0.5,nomegas),rep(0.5,6)),
                     gamma1 = gamma1, gamma = c(gamma0,0,gamma3), alpha1 = alpha1, alpha0 = alpha0, sigma_e = sigma_e,
                     fixgamma1=0, fixgamma=c(0,1,0), fixsigma_e=0, fixalpha1=0, fixalpha0=0)
  sd2 <- sqrt(diag(cov2$covmatrix))
  }



betahat2 <- c(betahat2,gamma1,gamma0,gamma3,sigma_e,alpha1,alpha0)

Zvalue2 <- betahat2/sd2
pvalue2 <- 2*(1-pnorm(abs(Zvalue2)))

naivebeta2 <- c(naive.model1$coefficients,naive.model2$coefficients)
sdnaive2 <- sqrt(c(diag(vcov(naive.model1)),diag(vcov(naive.model1))))

Zvalue_naive2 <- naivebeta2/sdnaive2
pvalue_naive2 <- 2*(1-pnorm(abs(Zvalue_naive2)))

Table_numeric2 <- data.frame(naivebeta=c(naivebeta2,rep(0,8)),
                            sdnaive=c(sdnaive2,rep(0,8)),
                            Z_naive=c(Zvalue_naive2,rep(0,8)),
                            p_naive=c(pvalue_naive2,rep(0,8)),
                            propobeta=betahat2,
                            proposd = sd2,
                            propoZ = Zvalue2,
                            propop = pvalue2)



library(xtable)
xtable(Table_numeric2,digits = 3)

## 2.1 Naive Analysis of 331 ####

# Covariates1 <- Covariates2 <- as.matrix(cbind(Covariates_all[,c("Intercept","X331","PC1","PC2")])) 
# 
# beta_Y1_0 <- mean(Y1star)
# beta_Y2_0 <- log(mean(Y2star)/(1-mean(Y2star)))
# 
# 
# intial1<-c(beta_Y1_0,0,0,0,0,beta_Y2_0,0,0,0,0,1,0)
# 
# # naive.model1 <- lm(Y1star ~ 1)
# # naive.model2 <- glm(Y2star ~ 1)
# naive.model1 <- lm(Y1star ~ ., data = data.frame(Covariates1[,-1]))
# naive.model2 <- glm(Y2star ~ ., family = binomial(link = logit), data = data.frame(Covariates2[,-1]))
# 
# intial2 <- c(naive.model1$coefficients,naive.model2$coefficients,1,0)
# cat("Initial value:", intial2, "\n")
# 
# 
# 
# NR <- nleqslv(intial2, GEE_UI, Y1star=Y1star, Y2star=Y2star, DesignMatrix1=Covariates1, DesignMatrix2=Covariates2,
#               CovMis1=CovMis1, CovMis2=CovMis2, alpha1=-2, alpha0=-2, jacobian=T, control=list(maxit=2000), sigma_e=0.37)
# 
# betahat <- ifelse(abs(NR$x)<100,NR$x,NA)
# 
# if (!any(is.na(betahat))) {
#   cov <- GEE_cov(betahat, Y1star, Y2star, DesignMatrix1=Covariates1, DesignMatrix2=Covariates2, 
#                  CovMis1=CovMis1, CovMis2=CovMis2,
#                  gamma=c(0,0), alpha1=-2, alpha0=-2, sigma_e=0.37)
#   sd <- sqrt(diag(cov))
# }
# 
# Zvalue <- betahat/sd
# pvalue <- 2*(1-pnorm(abs(Zvalue)))
# 
# Results109 <- round(data.frame(estimate=betahat, std=sd, pvalue=sd),3)
# save(Results109,file="results/DataAnalysis/Results109.RData")
