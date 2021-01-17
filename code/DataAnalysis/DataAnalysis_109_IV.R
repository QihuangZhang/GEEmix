
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

load("N:/Work/waterloo/2018/GEEmix/file/F_PC.RData")
# loading <- F_PC$x
# save(loading,file="loading.RData")
# screeplot(F_PC)  # obtain the screeplot


# obtain the screeplot
Scree <- data.frame(variance = (F_PC$sdev)^2, percentage = cumsum((F_PC$sdev)^2/sum((F_PC$sdev)^2)))
Scree <- Scree[1:10,]
rownamesa <- NULL
for (i in 1:10){
  rownamesa <- c(rownamesa, paste0("PC",i))
}

Scree <- data.frame(rownamesa = rownamesa,Scree)
Scree$rownamesa <-  factor(Scree$rownamesa,levels = Scree$rownamesa)
library(ggplot2)

# ggplot(Scree)  +
#   geom_bar(aes(x=rownamesa, y=variance),stat="identity")+
#   geom_line(aes(x=rownamesa, y=percentage*max(Scree$variance)*10,group = 1),stat="identity")+
#   labs(x = "Principal Components", y = "Variance") +
#   scale_y_continuous(sec.axis = sec_axis(~./max(Scree$variance)/10, name ="Cumulative percentage"))

# ggplot(Scree)  +
#   geom_bar(aes(x=rownamesa, y=variance),stat="identity")+
#   labs(x = "Principal Components", y = "Variance")



load("~/Work/waterloo/2018/mixResponse/file/genotrunc.RData")
phenotype<-read.csv(file="~/Work/waterloo/2018/mixResponse/data/pheno.csv")
load("file/loading.RData")

### Note:!! if use windows, change ~/ back to P:/ or N:/



## 0.2 Prepare package of GWAS toolbox ####
# library(maxLik)
library(MASS)
library(GeneErrorMis)
library(parallel)
library(nleqslv)
# source("code/CGRM.R")

## 0.3 Prepare the datasets ####
## Association Data Frame
data_GWA_Pre <- data.frame(id = phenotype$id,Y1 = phenotype$tibia,Y2 = phenotype$abnormalbone, discard = ifelse(phenotype$discard=="no",1,0),
                           batch16 = as.numeric(phenotype$round=="SW16"), BW = phenotype$bw0, BMD=phenotype$BMD)


data_GWAS_unsorted <- data_GWA_Pre[data_GWA_Pre$id %in% genotype$id,]
data_GWAS <- data_GWAS_unsorted[order(data_GWAS_unsorted$id),]

# Response
Y1star <- data_GWAS$Y1
Y2star <- data_GWAS$Y2

## Covariates for main model
Covariates<-cbind(Intercept=rep(1,length(genotype[,2])),X109=genotype[,2],X331=genotype[,3],data_GWAS[,c("batch16","BW","BMD")],
                  PC1 =loading[,2], PC2 =loading[,3])
Covariates<-as.matrix(Covariates)




## Mark the complete sample
comp<-data.frame(Y1star,Y2star,Covariates)
index.notna<-complete.cases(comp)
index.notna <- ifelse(data_GWAS$discard == 1, index.notna , FALSE)

Y1<-Y1star[index.notna]
Y2<-Y2star[index.notna]
Covariates_all<-Covariates[index.notna,]
nsample<-dim(Covariates_all)[1]


CovMis1 <- as.matrix(data.frame(intercept=rep(0,nsample),cov=rep(0,nsample)))
  
CovMis2 <- as.matrix(rep(1,nsample))




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
                      gamma1, gamma, alpha1, alpha0, sigma_e){
  # cat(theta, " \n")
  return(GEE_UfuncInsIV(Y1star=data.mismeasure$Y1star,
                        Y2star=data.mismeasure$Y2star,
                        Y1 = data.validation$Y1,
                        Y2 = data.validation$Y2,
                        DesignMatrix1 = as.matrix(data.mismeasure[,3:6]), 
                        DesignMatrix2 = as.matrix(data.mismeasure[,3:6]), 
                        ValidationMatrix1 = as.matrix(data.validation[,5:8]),
                        ValidationMatrix2 = as.matrix(data.validation[,5:8]),
                        CovMis1 = as.matrix(data.mismeasure[,7:8]),
                        CovMis2 = as.matrix(data.mismeasure[,9]),
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

GEE_covIV <- function(theta, data.validation, data.mismeasure,
                      gamma1, gamma, alpha1, alpha0, sigma_e,
                      fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0){
  nvalidation <- dim(data.validation)[1]
  nsample <- dim(data.mismeasure)[1] + nvalidation
  
  M0 <- GEE_GAMMA_IV0(theta, 
                      Y1star=data.validation$Y1star, 
                      Y2star=data.validation$Y2star, 
                      Y1 = data.validation$Y1, 
                      Y2 = data.validation$Y2, 
                      ValidationMatrix1 = as.matrix(data.validation[,5:8]),
                      ValidationMatrix2 = as.matrix(data.validation[,5:8]), 
                      CovMis1 = as.matrix(data.validation[,9:10]), 
                      CovMis2 = as.matrix(data.validation[,11]), 
                      gamma1, gamma, alpha1, alpha0, sigma_e,
                      fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0)
  
  M1 <- GEE_GAMMA_IVI(theta, 
                      Y1star=data.mismeasure$Y1star, 
                      Y2star=data.mismeasure$Y2star, 
                      DesignMatrix1 = as.matrix(data.mismeasure[,3:6]),
                      DesignMatrix2 = as.matrix(data.mismeasure[,3:6]), 
                      CovMis1 = as.matrix(data.mismeasure[,7:8]), 
                      CovMis2 = as.matrix(data.mismeasure[,9]), 
                      gamma1, gamma, alpha1, alpha0, sigma_e,
                      fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0)
  
  GAMMA_IV <- M1 + M0
  
  B0 <- GEE_SIGMA_IV0(theta, 
                      Y1star=data.validation$Y1star, 
                      Y2star=data.validation$Y2star, 
                      Y1 = data.validation$Y1, 
                      Y2 = data.validation$Y2, 
                      ValidationMatrix1 = as.matrix(data.validation[,5:8]),
                      ValidationMatrix2 = as.matrix(data.validation[,5:8]), 
                      CovMis1 = as.matrix(data.validation[,9:10]), 
                      CovMis2 = as.matrix(data.validation[,11]), 
                      gamma1, gamma, alpha1, alpha0, sigma_e,
                      fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0)
  
  B1 <- GEE_SIGMA_IVI(theta, 
                      Y1star=data.mismeasure$Y1star, 
                      Y2star=data.mismeasure$Y2star, 
                      DesignMatrix1 = as.matrix(data.mismeasure[,3:6]),
                      DesignMatrix2 = as.matrix(data.mismeasure[,3:6]), 
                      CovMis1 = as.matrix(data.mismeasure[,7:8]), 
                      CovMis2 = as.matrix(data.mismeasure[,9]), 
                      gamma1, gamma, alpha1, alpha0, sigma_e,
                      fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0)
  
  SIGMA_IV <- B1 + B0
  
  GAMMA.inv <- solve(GAMMA_IV,tol=1e-200)
  covmatrix <- GAMMA.inv %*% SIGMA_IV %*% t(as.matrix(GAMMA.inv))
  return(covmatrix)
}

#### 1. Data Encryption ####

Covariates1 <- Covariates2 <- as.matrix(cbind(Covariates_all[,c("Intercept","X109","PC1","PC2")])) 


noise <- rnorm(nsample,0,0.1)
noise2 <- runif(nsample,0,1)

Y1star <- Y1 + noise
Y2star <- ifelse(noise2>0.05,Y2,1-Y2)

CovMis1 <- as.matrix(data.frame(intercept=rep(0,nsample),cov=rep(0,nsample)))

CovMis2 <- as.matrix(rep(1,nsample))

mismeasuredata <- data.frame(Y1star=Y1star[401:nsample], Y2star=Y2star[401:nsample],
                             Covariates1[401:nsample,],CovMis1[401:nsample,],CovMis2[401:nsample,])
validationdata <- data.frame(Y1star=Y1star[1:400], Y2star=Y2star[1:400],
                             Y1=Y1[1:400], Y2=Y2[1:400],
                             Covariates1[1:400,],CovMis1[1:400,],CovMis2[1:400,])


#### 2. Analysis  of 109 ####


## 2.1 Naive Analysis of 109 ####


beta_Y1_0 <- mean(mismeasuredata$Y1star)
beta_Y2_0 <- log(mean(mismeasuredata$Y2star)/(1-mean(mismeasuredata$Y2star)))



## 2.2 Estimating Procedure ####

# naive.model1 <- lm(Y1star ~ 1)
# naive.model2 <- glm(Y2star ~ 1)

## 2.2.1 Naive model ####
naive.model1 <- lm(Y1star ~ X109 + PC1 + PC2, data = mismeasuredata)
naive.model2 <- glm(Y2star ~ X109 + PC1 + PC2, family = binomial(link = logit), data = mismeasuredata)

intial2 <- c(naive.model1$coefficients,naive.model2$coefficients,0.03,0)
cat("Initial value:", intial2, "\n")


naive.model3 <- lm(Y1 ~ X109 + PC1 + PC2, data = validationdata)
naive.model4 <- glm(Y2 ~ X109 + PC1 + PC2, family = binomial(link = logit), data = validationdata)

intial3 <- c(naive.model3$coefficients,naive.model4$coefficients,1,0)
cat("Initial value:", intial3, "\n")


## 2.2.2 the proposed approach ####

## Measurement Error and Misclassification Parameters
model.measure <- lm(Y1star ~ -1 + offset(Y1),data = validationdata) 
model.class1 <- glm((1-Y2star) ~ 1, data = validationdata[validationdata$Y2==1,],family = binomial(link="logit")) 
model.class0 <- glm(Y2star ~ 1, data = validationdata[validationdata$Y2==0,],family = binomial(link="logit")) 

# gamma1 <- model.measure$coefficients[1] 
sigma_e <- sigma(model.measure)
alpha1 <- model.class1$coefficients
alpha0 <- model.class0$coefficients

## Mainly interested parameters

NR <- nleqslv(intial2, GEE_UI_IV, jacobian=T, control=list(maxit=4000),
              data.validation = validationdata, data.mismeasure = mismeasuredata, 
              gamma1 = 1, gamma=c(0,0), alpha1= alpha1, alpha0= alpha0, sigma_e = sigma_e)

betahat <- NR$x

### variance estimation with validation data
if (!any(is.na(betahat))) {
  cov <- GEE_covIV (betahat, data.validation= validationdata, data.mismeasure = mismeasuredata,
                    gamma1=1, gamma = c(0,0), alpha1, alpha0, sigma_e,
                    fixgamma1=1, fixgamma=c(1,1), fixsigma_e=0, fixalpha1=0, fixalpha0=0)
  sd <- sqrt(diag(cov))} else {
    sd <- rep(NA,length(betahat))
  }

betahat <- c(betahat,sigma_e,alpha1,alpha0)

Zvalue <- betahat/sd
pvalue <- 2*(1-pnorm(abs(Zvalue)))

naivebeta <- c(naive.model1$coefficients,naive.model2$coefficients)
sdnaive <- sqrt(c(diag(vcov(naive.model1)),diag(vcov(naive.model1))))

Zvalue_naive <- naivebeta/sdnaive
pvalue_naive <- 2*(1-pnorm(abs(Zvalue_naive)))

Table_numeric <- data.frame(naivebeta=c(naivebeta,rep(0,5)),
                            sdnaive=c(sdnaive,rep(0,5)),
                            Z_naive=c(Zvalue_naive,rep(0,5)),
                            p_naive=c(pvalue_naive,rep(0,5)),
                            propobeta=betahat,
                            proposd = sd,
                            propoZ = Zvalue,
                            propop = pvalue)

library(xtable)
xtable(Table_numeric,digits = 3)

#### 2. Analysis  of 331 ####


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
