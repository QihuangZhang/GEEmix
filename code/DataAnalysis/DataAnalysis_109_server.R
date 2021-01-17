
### Project GEEmix

#### 0 Project Information and Configuration ####
ProjectName<- paste0("DataAnalysis_109")
cat(paste0("Note: This is for the data analysis of SNP 109 \n")) 

## 0.1 Set file path ####
setwd("~/Work/waterloo/2018/GEEmix")
# WD=paste0(getwd(),"/results/EMsensi/")


## Finished part -- no need to run again
# genotype<-read.table(file="data/geno.txt",header=T)
# save(genotype,file="genotype.RData")
# load("P:/Work/waterloo/2018/mixResponse/file/genotype.RData")
# load("P:/Work/waterloo/2018/mixResponsefile/RR.RData")
# F_PC <- prcomp(genotype[,3:dim(genotype)[2]],center = T,scale.=T)
# load("P:/Work/waterloo/2018/GEEmix/file/F_PC.RData")
# loading <- F_PC$x
# save(loading,file="loading.RData")
# screeplot(F_PC)  # obtain the screeplot

load("~/Work/waterloo/2018/mixResponse/file/genotrunc.RData")
phenotype<-read.csv(file="~/Work/waterloo/2018/mixResponse/data/pheno.csv")
load("~/Work/waterloo/2018/GEEmix/file/loading.RData")




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

Y1star<-Y1star[index.notna]
Y2star<-Y2star[index.notna]
Covariates_all<-Covariates[index.notna,]
nsample<-dim(Covariates_all)[1]


CovMis1 <- as.matrix(data.frame(intercept=rep(0,nsample),cov=rep(0,nsample)))
  
CovMis2 <- as.matrix(rep(1,nsample))




## 0.4 Function set up ####
GEE_UI <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, alpha1, alpha0, sigma_e){
  # cat(theta, " \n")
  return(GEE_UfuncIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                      theta[1:4], theta[5:8], phi = theta[9], xi = theta[10], 
                      gamma1 = 1, gamma=c(0,0),  alpha1, alpha0, sigma_e
                      ))
}

GEE_SIGMA <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, gamma, alpha1, alpha0,sigma_e){
  return(GEE_SIGMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                      theta[1:4], theta[5:8], phi = theta[9], xi = theta[10], 
                      gamma1 = 1, gamma, alpha1, alpha0, sigma_e))
}

GEE_GAMMA <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2){
  GAMMA <- GEE_GAMMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, 
                        beta1=theta[1:4], beta2=theta[5:8], phi = theta[9], xi = theta[10])
  return(GAMMA)
}

GEE_GAMMA.inv <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2){
  GAMMA <- GEE_GAMMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, 
                        beta1=theta[1:4], beta2=theta[5:8], phi = theta[9], xi = theta[10])
  GAMMA.inv <- solve(GAMMA,tol=1e-200)
  return(GAMMA.inv)
}

GEE_cov <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, gamma, alpha1, alpha0, sigma_e){
  GAMMA.inv <- GEE_GAMMA.inv(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
  SIGMA <- GEE_SIGMA(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, gamma, alpha1, alpha0, sigma_e)
  covmatrix <- GAMMA.inv %*% SIGMA %*% t(as.matrix(GAMMA.inv))
  return(covmatrix)
}

expit <- function(x){
    value <- exp(x)/(1+exp(x))
    ifelse(is.na(value),1,value) 
}

logit <- function(x){
  value <- log(x/(1-x))
  return(value)
}

#### 1. Analysis  of 109 ####


## 1.1 Naive Analysis of 109 ####

Covariates1 <- Covariates2 <- as.matrix(cbind(Covariates_all[,c("Intercept","X109","PC1","PC2")])) 

beta_Y1_0 <- mean(Y1star)
beta_Y2_0 <- log(mean(Y2star)/(1-mean(Y2star)))


intial1<-c(beta_Y1_0,0,0,0,0,beta_Y2_0,0,0,0,0,1,0)

# naive.model1 <- lm(Y1star ~ 1)
# naive.model2 <- glm(Y2star ~ 1)

## 1.2 the proposed approach ####
naive.model1 <- lm(Y1star ~ ., data = data.frame(Covariates1[,-1]))
naive.model2 <- glm(Y2star ~ ., family = binomial(link = logit), data = data.frame(Covariates2[,-1]))

intial2 <- c(naive.model1$coefficients,naive.model2$coefficients,1,0)
cat("Initial value:", intial2, "\n")

Large109 <- NULL

for (sigma_e in c(0.72,0.77,0.82)){
  Results109 <- c("(Intercept1)","X","PC1","PC2","(Intercept2)","X","PC1","PC2","sigma_e","xi")
  for (alphavalue in c(-2.442347, -2.376273, -2.313635, -2.254058, -2.197225)){
    NR <- nleqslv(intial2, GEE_UI, Y1star=Y1star, Y2star=Y2star, DesignMatrix1=Covariates1, DesignMatrix2=Covariates2,
                  CovMis1=CovMis1, CovMis2=CovMis2, alpha1=alphavalue, alpha0=alphavalue, jacobian=T, control=list(maxit=2000), sigma_e=sigma_e)
    
    betahat <- ifelse(abs(NR$x)<100,NR$x,NA)
    
    if (!any(is.na(betahat))) {
      cov <- GEE_cov(betahat, Y1star, Y2star, DesignMatrix1=Covariates1, DesignMatrix2=Covariates2, 
                     CovMis1=CovMis1, CovMis2=CovMis2,
                     gamma=c(0,0), alpha1=alphavalue, alpha0=alphavalue, sigma_e=sigma_e)
      sd <- sqrt(diag(cov))
    }
    
    Zvalue <- betahat/sd
    pvalue <- 2*(1-pnorm(abs(Zvalue)))
    
    Results109 <- data.frame(Results109,round(data.frame(alphavalue = alphavalue, estimate=betahat, std=sd, pvalue=pvalue),3))
  }
  Large109 <- rbind(Large109, Results109)
}

library(xtable)
xtable(Large109,digits = 3)


Large109 <- NULL

for (sigma_e in c(0.77)){
  Results109 <- c("(Intercept1)","X","PC1","PC2","(Intercept2)","X","PC1","PC2","sigma_e","xi")
  for (alphavalue in c(-2.338303, -2.325907, -2.313635, -2.301486, -2.289456)){
    NR <- nleqslv(intial2, GEE_UI, Y1star=Y1star, Y2star=Y2star, DesignMatrix1=Covariates1, DesignMatrix2=Covariates2,
                  CovMis1=CovMis1, CovMis2=CovMis2, alpha1=alphavalue, alpha0=alphavalue, jacobian=T, control=list(maxit=2000), sigma_e=sigma_e)
    
    betahat <- ifelse(abs(NR$x)<100,NR$x,NA)
    
    if (!any(is.na(betahat))) {
      cov <- GEE_cov(betahat, Y1star, Y2star, DesignMatrix1=Covariates1, DesignMatrix2=Covariates2, 
                     CovMis1=CovMis1, CovMis2=CovMis2,
                     gamma=c(0,0), alpha1=alphavalue, alpha0=alphavalue, sigma_e=sigma_e)
      sd <- sqrt(diag(cov))
    }
    
    Zvalue <- betahat/sd
    pvalue <- 2*(1-pnorm(abs(Zvalue)))
    
    Results109 <- data.frame(Results109,round(data.frame(alphavalue = alphavalue, estimate=betahat, std=sd, pvalue=pvalue),3))
  }
  Large109 <- rbind(Large109, Results109)
}

library(xtable)
xtable(Large109,digits = 3)
# save(Results109,file="results/DataAnalysis/Results109.RData")

# ResultsDA<-cbind(Results109,Results331)
# library(xtable)
# xtable(ResultsDA)



# #### 2. Analysis  of 331 ####
# 
# 
# ## 2.1 Naive Analysis of 331 ####
# 
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
# # save(Results109,file="results/DataAnalysis/Results109.RData")
