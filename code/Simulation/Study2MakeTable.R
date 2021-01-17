####
# This file collect the information from all tables and produce a summary table

#### 0 Simulation set up####
library(parallel)
library(scales)
library(xtable)
library(tidyr)
library(ggplot2)

#### 1. Package the data ####

 

for (xx in 1:3) {load(paste0("output/KNOWN_R",xx,".RData"))}
KNOWN <- list(results_1, results_2, results_3)

for (xx in 1:3) {load(paste0("output/INS3_R",xx,".RData"))}
INTERNAL <- list(results_1, results_2, results_3)

for (xx in 1:3) {load(paste0("output/INSEV_R",xx,".RData"))}
EXTERNAL <- list(results_1, results_2, results_3)


#### 2. Main data ####
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

MakeTable <- function(xx){
  Results <- NULL
  Results1 <- NULL
  
  truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0)
  
  if (xx==1) {
    range <- c(0.1,0.5,0.7)
  } else if (xx==2) {
    range <- c(-4.595,-2.197,-1.386)
  } else if (xx==3) {
    range <- c(0,-0.8,0.8)
  }
  
  #### 2.1 The naive table and know table ####
  
  for (k in 1:3) {
    results <- KNOWN[[xx]][[k]]
    naive1coef <- NULL
    naive1sd <- NULL
    naive2coef <- NULL
    naive2sd <- NULL
    CI1naive <- NULL
    CI2naive <- NULL
    
    betas <- NULL
    sds <- NULL
    CIs <- NULL
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
    sd_emp <- apply(na.omit(betas),MARGIN = 2, FUN = function(x){  
      x.noout <-  remove_outliers(x)   
      return( sd(x.noout,na.rm = T ))   })
    sd_mod <-  apply(na.omit(sds),MARGIN = 2, FUN = function(x){  
      x.noout <-  remove_outliers(x) 
      return( mean(x.noout,na.rm = T ))   })
    CIrate <- colMeans(na.omit(CIs),na.rm = T)
    
    Results0 <- data.frame( range = range[k],
                            naivebias=c(round(naivebias,3),rep(NA,4)),
                            naive_esd=c(round(naive_esd,3),rep(NA,4)),
                            naivesd=c(round(naivesd,3),rep(NA,4)),
                            naiveCI=c(  percent(naiveCI, accuracy = .1),rep(NA,4)),
                            biasprop=c(round(bias1,3),rep(NA,4)),
                            propose_esd=c(round(sd_emp,3),rep(NA,4)),
                            sdpropose=c(round(sd_mod,3),rep(NA,4)),
                            CI_propose=c(  percent(CIrate, accuracy = .1),rep(NA,4)))
    
    Results1 <- rbind(Results1,Results0)
  }
  
  Results <- Results1
  
  Results1 <- NULL
  
  for (k in 1:3) {
    results <- INTERNAL[[xx]][[k]]
    if (xx==1) {
      truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0,0.8,range[k],-2.197,-2.197)
    } else if (xx==2) {
      truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0,0.8,0.1,range[k],range[k])
    } else if (xx==3) {
      truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0,range[k],0.1,-2.197,-2.197)
    }
    
    betas <- NULL
    sds <- NULL
    CIs <- NULL
    
    for (i in 1:1000){
      if (is.null(results[[i]])) {
        next}
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
    
    bias1 <- colMeans(na.omit(betas),na.rm = T) - truebeta
    sd_emp <- apply(na.omit(betas),MARGIN = 2, FUN = function(x){  
      x.noout <-  remove_outliers(x) 
      return( sd(x.noout,na.rm = T ))   })
    sd_mod <-  apply(na.omit(sds),MARGIN = 2, FUN = function(x){ 
      x.noout <-  remove_outliers(x)  
      return( mean(x.noout,na.rm = T ))   })
    CIrate <- colMeans(na.omit(CIs),na.rm = T)
    
    Results0 <- data.frame(biasprop=round(bias1,3),propose_esd=round(sd_emp,3),sdpropose=round(sd_mod,3),CI_propose=percent(CIrate, accuracy = .1))
    
    Results1 <- rbind(Results1,Results0)
  }
  
  Results <- cbind(Results,Results1)
  
  Results1 <- NULL
  
  for (k in 1:3) {
    results <- EXTERNAL[[xx]][[k]]
    if (xx==1) {
      truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0,0.8,range[k],-2.197,-2.197)
    } else if (xx==2) {
      truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0,0.8,0.1,range[k],range[k])
    } else if (xx==3) {
      truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0,range[k],0.1,-2.197,-2.197)
    }
    
    betas <- NULL
    sds <- NULL
    CIs <- NULL
    
    for (i in 1:1000){
      if (is.null(results[[i]])) {
        next}
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
    
    bias1 <- colMeans(na.omit(betas),na.rm = T) - truebeta
    sd_emp <- apply(na.omit(betas),MARGIN = 2, FUN = function(x){
      x.noout <-  remove_outliers(x)    
      return( sd(x.noout,na.rm = T ))   })
    sd_mod <-  apply(na.omit(sds),MARGIN = 2, FUN = function(x){
      x.noout <-  remove_outliers(x)   
      return( mean(x.noout,na.rm = T ))   })
    CIrate <- colMeans(na.omit(CIs),na.rm = T)
    
    Results0 <- data.frame(biasprop=round(bias1,3),propose_esd=round(sd_emp,3),sdpropose=round(sd_mod,3),CI_propose=percent(CIrate, accuracy = .1))
    
    Results1 <- rbind(Results1,Results0)
  }
  
  Results <- cbind(Results,Results1)
  
  return(Results)
}


#### 3 Print the Data ####

Table4 <- MakeTable(1)
xtable(Table4,digits = 3)


Table5 <- MakeTable(2)
xtable(Table5,digits = 3)


Table6 <- MakeTable(3)
xtable(Table6,digits = 3)




#### 4 Produce the figure ####
### 4.1 Function for figure ####
MakeTable4fig <- function(xx){
  Results <- NULL
  Results1 <- NULL
  
  truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0)
  
  if (xx==1) {
    range <- c(0.1,0.5,0.7)
  } else if (xx==2) {
    range <- c(-4.595,-2.197,-1.386)
  } else if (xx==3) {
    range <- c(0,-0.8,0.8)
  }
  
  #### 2.1 The naive table and know table ####
  
  for (k in 1:3) {
    results <- KNOWN[[xx]][[k]]
    naive1coef <- NULL
    naive1sd <- NULL
    naive2coef <- NULL
    naive2sd <- NULL
    CI1naive <- NULL
    CI2naive <- NULL
    
    betas <- NULL
    sds <- NULL
    CIs <- NULL
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
    
    biasnaive1 <- abs(colMeans(naive1coef)-truebeta[1:3])
    biasnaive2 <- abs(colMeans(naive2coef)-truebeta[4:6])
    naive_esd <- apply(cbind(naive1coef,naive2coef), MARGIN = 2 , FUN=sd, na.rm=T)
    sdnaive1 <- colMeans(naive1sd)
    sdnaive2 <- colMeans(naive2sd)
    CInaive1 <- colMeans(CI1naive,na.rm=T)
    CInaive2 <- colMeans(CI2naive,na.rm=T)
    
    naivebias <- c(biasnaive1,biasnaive2,0,0)
    naive_esd <- c(naive_esd,0,0)
    naivesd <- c(sdnaive1,sdnaive2,0,0)
    naiveCI <- c(CInaive1,CInaive2,0,0)
    
    bias1 <- abs(colMeans(na.omit(betas),na.rm = T) - truebeta)
    sd_emp <- apply(na.omit(betas),MARGIN = 2, FUN = function(x){ 
      x.noout <-  remove_outliers(x)    
      return( sd(x.noout,na.rm = T ))   })
    sd_mod <-  apply(na.omit(sds),MARGIN = 2, FUN = function(x){ 
      x.noout <-  remove_outliers(x)  
      return( mean(x.noout,na.rm = T ))   })
    CIrate <- colMeans(na.omit(CIs),na.rm = T)
    
    Results0 <- data.frame( range = range[k],
                            naivebias=c(round(naivebias,3),rep(NA,4))/truebeta,
                            naive_esd=c(round(naive_esd,3),rep(NA,4)),
                            naivesd=c(round(naivesd,3),rep(NA,4)),
                            naiveCI=c(  (round(naiveCI,3)),rep(NA,4)),
                            biasprop=c(round(bias1,3),rep(NA,4))/truebeta,
                            propose_esd=c(round(sd_emp,3),rep(NA,4)),
                            sdpropose=c(round(sd_mod,3),rep(NA,4)),
                            CI_propose=c(  percent(CIrate, accuracy = .1),rep(NA,4)))
    
    Results1 <- rbind(Results1,Results0)
  }
  
  Results <- Results1
  
  Results1 <- NULL
  
  for (k in 1:3) {
    results <- INTERNAL[[xx]][[k]]
    if (xx==1) {
      truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0,0.8,range[k],-2.197,-2.197)
    } else if (xx==2) {
      truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0,0.8,0.1,range[k],range[k])
    } else if (xx==3) {
      truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0,range[k],0.1,-2.197,-2.197)
    }
    
    betas <- NULL
    sds <- NULL
    CIs <- NULL
    
    for (i in 1:1000){
      if (is.null(results[[i]])) {
        next}
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
    
    bias1 <- abs((colMeans(na.omit(betas),na.rm = T) - truebeta)/truebeta)
    sd_emp <- apply(na.omit(betas),MARGIN = 2, FUN = function(x){   
      x.noout <-  remove_outliers(x)    
      return( sd(x.noout,na.rm = T ))   })
    sd_mod <-  apply(na.omit(sds),MARGIN = 2, FUN = function(x){ 
      x.noout <-  remove_outliers(x)    
      return( mean(x.noout,na.rm = T ))   })
    CIrate <- colMeans(na.omit(CIs),na.rm = T)
    
    Results0 <- data.frame(biasprop=round(bias1,3),propose_esd=round(sd_emp,3),sdpropose=round(sd_mod,3),CI_propose=percent(CIrate, accuracy = .1))
    
    Results1 <- rbind(Results1,Results0)
  }
  
  Results <- cbind(Results,Results1)
  
  Results1 <- NULL
  
  for (k in 1:3) {
    results <- EXTERNAL[[xx]][[k]]
    if (xx==1) {
      truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0,0.8,range[k],-2.197,-2.197)
    } else if (xx==2) {
      truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0,0.8,0.1,range[k],range[k])
    } else if (xx==3) {
      truebeta <- c(0.7,1.5,-1,0.7,-1.5,1,1,0,range[k],0.1,-2.197,-2.197)
    }
    
    betas <- NULL
    sds <- NULL
    CIs <- NULL
    
    for (i in 1:1000){
      if (is.null(results[[i]])) {
        next}
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
    
    bias1 <- abs((colMeans(na.omit(betas),na.rm = T) - truebeta)/truebeta)
    sd_emp <- apply(na.omit(betas),MARGIN = 2, FUN = function(x){ 
      x.noout <-  remove_outliers(x)    
      return( sd(x.noout,na.rm = T ))   })
    sd_mod <- apply(na.omit(sds),MARGIN = 2, FUN = function(x){ 
      x.noout <-  remove_outliers(x)  
      return( mean(x.noout,na.rm = T ))   })
    CIrate <- colMeans(na.omit(CIs),na.rm = T)
    
    Results0 <- data.frame(biasprop=round(bias1,3),propose_esd=round(sd_emp,3),sdpropose=round(sd_mod,3),CI_propose=percent(CIrate, accuracy = .1))
    
    Results1 <- rbind(Results1,Results0)
  }
  
  Results <- cbind(Results,Results1)
  
  return(Results)
}


### 1.2 Data Preperation ####
Table4 <- MakeTable4fig(1)

Table5 <- MakeTable4fig(2)

Table6 <- MakeTable4fig(3)


### 1.3 Data integration ####


# parnames <- c("beta10","beta11","beta12","beta20","beta21","beta22","","","","","","")
# response0 <- c(rep("Continuous",3),rep("Binary",3),rep("",6))
# colnames(Table5) <- c("range",rep(c("bias","esd", "sdpropose","CI"),4))
# 
# FigData2 <- data.frame(par=rep(parnames,3),response=rep(response0,3),range=Table5[,1],Table5[2:5],Val="Naive")
# FigData2 <- rbind(FigData2,data.frame(par=rep(parnames,3),response=rep(response0,3),range=Table5[,1], Table5[6:9],Val="Known (Case 1)"))
# FigData2 <- rbind(FigData2,data.frame(par=rep(parnames,3),response=rep(response0,3),range=Table5[,1],Table5[10:13],Val="Internal (Case 2)"))
# FigData2 <- FigData2[!FigData2$par=="",c("par","response","range", "bias", "CI", "Val")]
# colnames(FigData2) <- c("Parameter","Response","Misclassification","Bias","Coverage","Method")
# 
# FigData2$Misclassification[FigData2$Misclassification==-4.595] <- "1%"
# FigData2$Misclassification[FigData2$Misclassification==-2.197] <- "10%"
# FigData2$Misclassification[FigData2$Misclassification==-1.386] <- "20%"
# 
# FigData_long2 <- gather(FigData2, Type , Value, c(Bias,Coverage))
# 
# FigData_long2$Type[FigData_long2$Type=="Coverage"] <- "Coverage Rate"
# FigData_long2$Value <- abs(FigData_long2$Value)
# ### 2.2  Plot the Figure ####
# # Grouped
# custom.col <- c("#FFDB6D",  "#C3D7A4", 
#                 "#56B4E9", "#E69F00","#52854C", "#293352")
# pdf(file = "output/SSC.pdf",height = 5.6, width = 9.3)
# figSSC <- ggplot()+geom_bar(data=subset(FigData_long2,Type=="Bias"), aes(fill=Parameter, y=Value, x=Misclassification),position="dodge", stat="identity") + 
#   geom_point(data=subset(FigData_long2,Type=="Coverage Rate"),aes(colour = Parameter, y=Value, x=Misclassification),size=3) + 
#   facet_grid(Type~Method,scales="free_y",switch ="y") +
#   geom_hline(data=subset(FigData_long2,Type=="Coverage Rate"),aes(yintercept=0.95),linetype="dashed", color = "red")+
#   scale_fill_manual(values=custom.col) +
#   scale_colour_manual(values=custom.col) +
#   scale_y_continuous(labels=scales::percent_format(),position = "right") + 
#   theme(text=element_text(size=14, family="sans")) +
#   labs(x = "Misclassification Rate", y="Percentage")
# figSSC
# dev.off()


