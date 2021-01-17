## This file is the start of the project, targeting on the correctness of the code.


#### 1 Data Generation ####

## parameters
beta1 <- c(0.7,1.5,-1)
beta2 <- c(0.7,-1.5,1)

theta <- c(beta1, beta2, phi, xi)
phi <- 1
xi <- 0.8

### Generate the true data sets

set.seed(2019)
nsample <- 100000
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

Y1star <- Y1 + 0.2 * Y2+  e
Y2star <- ifelse(U2>0.2,Y2,1-Y2)

## Naive model
naive.model1 <- lm(Y1star ~ X + W)
true.model1 <- lm(Y1 ~ X + W)
naive.model2 <- glm(Y2star ~ X + W, family = binomial(link = logit))
true.model2 <- glm(Y2star ~ X + W, family = binomial(link = logit))

#### 2 Implementation of the methods #####

## 2.1 Preparation ###

DesignMatrix1 <-  DesignMatrix2 <- cbind(rep(1,length(X)),X,W)
CovMis1 <- cbind(rep(0,length(X)),rep(1,length(X)))
CovMis2 <- c(rep(1,length(X)))

DesignMatrix1 <- as.matrix(DesignMatrix1)
DesignMatrix2 <- as.matrix(DesignMatrix2)
CovMis1 <- as.matrix(CovMis1)
CovMis2 <- as.matrix (CovMis2)

rho <- 0

results <- GEE_UfuncIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
             beta1, beta2, xi = rho, 
             sigma = 1, gamma1 = 1, gamma=c(0,0.2), alpha1=-1.386294, alpha0=-1.386294,
             sigma_e = 0.5)
# Delta <- results$Delta
# RHS <- mean((Y2-results$Y2starstar)^2)
# mean(Y1*Y2)
# 0.2*Delta
# mean(results$Y1starstar*results$Y2starstar)

GEE_UI <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2){
  cat(theta, " \n")
  # if ((theta[8]>1)|(theta[8]<-1)) {theta[8] <- 0}
  return(GEE_UfuncIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                      theta[1:3], theta[4:6], phi = theta[7], xi = theta[8], 
                       gamma1 = 1, gamma=c(0,0), alpha1=-1.386294, alpha0=-1.386294,
                      sigma_e = 0.5))
}

library(rootSolve)

GEE_UI(intial4, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
GEE_UI(intial2, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
GEE_UI(intial3, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)

intial <- c(0,0,0,0,0,0,1,0)
intial2 <- c(beta1,beta2,1,0)
intial3 <- c(0,0,0,1,0,0,1,0)
intial4 <- c(naive.model1$coefficients,naive.model2$coefficients,1,0)

NR <- multiroot(GEE_UI,intial4,Y1star=Y1star, Y2star=Y2star, DesignMatrix1=DesignMatrix1, DesignMatrix2=DesignMatrix2,
                       CovMis1=CovMis1, CovMis2=CovMis2)


##### 3 The correlated mismeasruement ####
## 3.1 Data Generation ####

gamma <- 0.8

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


## 3.2 Implementation Generation ###

DesignMatrix1 <-  DesignMatrix2 <- cbind(rep(1,length(X)),X,W)
CovMis1 <- cbind(rep(0,length(X)),rep(1,length(X)))
CovMis2 <- c(rep(1,length(X)))

DesignMatrix1 <- as.matrix(DesignMatrix1)
DesignMatrix2 <- as.matrix(DesignMatrix2)
CovMis1 <- as.matrix(CovMis1)
CovMis2 <- as.matrix (CovMis2)

rho <- 0

results <- GEE_UfuncIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                        beta1, beta2, xi = rho, 
                        phi = 1, gamma1 = 1, gamma=c(0,0.8), alpha1=-1.386294, alpha0=-1.386294,
                        sigma_e = sqrt(0.5))

GEE_UI2 <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2){
  cat(theta, " \n")
  # if ((theta[8]>1)|(theta[8]<-1)) {theta[8] <- 0}
  return(GEE_UfuncIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                      theta[1:3], theta[4:6], phi = theta[7], xi = theta[8], 
                      gamma1 = 1, gamma=c(0,0.8), alpha1=-1.386294, alpha0=-1.386294,
                      sigma_e = 0.5))
}

library(rootSolve)

intial <- c(0,0,0,0,0,0,1,0)
intial2 <- c(beta1,beta2,1,0)
intial3 <- c(0,0,0,1,0,0,1,0)
intial4 <- c(naive.model1$coefficients,naive.model2$coefficients,1,0.2)

GEE_UI2(intial, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
GEE_UI2(intial2, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
GEE_UI2(intial3, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
GEE_UI2(intial4, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)



NR <- multiroot(GEE_UI2,intial4,Y1star=Y1star, Y2star=Y2star, DesignMatrix1=DesignMatrix1, DesignMatrix2=DesignMatrix2,
                CovMis1=CovMis1, CovMis2=CovMis2)



### Debugging

## unbiasedness of Y1starstar  ...  Good
# mean(results - mu1)
## unbiasedness of Y2starstar  ...  Good
# mean(results - mu2expit)
## unbiasedness of Y1starstar  ...  Good



### Check if gamma function is corrected specified.

GEE_UI(intial, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)

initial_g1 <- NR$x
initial_g2 <- NR$x + c(0,0,0,0.001,0,0,0,0)

U1 <- GEE_UI(initial_g1, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
U2 <- GEE_UI(initial_g2, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
(U2-U1)/0.001

M2 <- GEE_GAMMA(initial_g1,Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2) * 1000

(NR$jac - M2)/NR$jac
