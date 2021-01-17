####
# This file collect the information from all tables and produce a summary table

#### 0 Simulation set up####
#### 0.1 Load R packages ####
library(parallel)
library(scales)
library(xtable)
library(tidyr)
library(ggplot2)
library(plyr)

#### 0.2 Global Parameters ####
parnames <- c("beta[10]","beta[11]","beta[12]","beta[20]","beta[21]","beta[22]","","","","","","")
response0 <- c(rep("Continuous",3),rep("Binary",3),rep("",6))


#### 1. Figure 2 - Sample Size ####
### 1.1 Data Preperation ####
load("P:/Work/waterloo/2018/GEEmix/WINV_RTable1.RData")
FigData <- data.frame(par=rep(parnames,3),response=rep(response0,3),Results,Val="Internal")
load("P:/Work/waterloo/2018/GEEmix/SSEV_RTable1.RData")
FigData <- rbind(FigData,data.frame(par=rep(parnames,3),response=rep(response0,3),Results,Val="External"))
FigData <- FigData[!FigData$par=="",c("par","response","nsample", "biasprop", "sdpropose", "Val")]
colnames(FigData) <- c("Parameter","Response","nsample","Bias","S.D.","Validation")

FigData_long <- gather(FigData, Type , Value, c(Bias,S.D.))

### 1.2  Plot the Figure ####
# Grouped
custom.col <- c("#FFDB6D",  "#C3D7A4", 
                "#56B4E9", "#E69F00","#52854C", "#293352")
pdf(file = "output/figure2.pdf",height = 5.5, width = 7)
f1 <- ggplot(FigData_long, aes(fill=Parameter, y=Value, x=Validation)) + 
  facet_grid(Type~nsample, scales = "free_y",switch ="y")+
  theme(text=element_text(size=13, family="mono"))+
  geom_bar(position="dodge", stat="identity") + 
  scale_fill_manual(values=custom.col) 
  

f1

dev.off()


#### 2. Figure 2 - Sample Size Rate ####
### 2.1 Data Preperation ####
load("P:/Work/waterloo/2018/GEEmix/WINV_RTable2.RData")
FigData2 <- data.frame(par=rep(parnames,3),response=rep(response0,3),Results,Val="Internal")
FigData2 <- FigData2[!FigData2$par=="",c("par","response","nrate", "biasprop", "sdpropose", "Val")]
load("P:/Work/waterloo/2018/GEEmix/SSEV_RTable2.RData")
FigData2 <- rbind(FigData2,data.frame(par=rep(parnames,3),response=rep(response0,3),Results,Val="External"))
FigData2 <- FigData2[!FigData2$par=="",c("par","response","nrate", "biasprop", "sdpropose", "Val")]
colnames(FigData2) <- c("Parameter","Response","nrate","Bias","S.D.","Validation")

FigData_long2 <- gather(FigData2, Type , Value, c(Bias,S.D.))
FigData_long2$nrate <- factor(FigData_long2$nrate,levels=c("0.333","0.5","0.667"),labels=c("1:2","1:1","2:1"))

### 2.2  Plot the Figure ####
# Grouped
custom.col <- c("#FFDB6D",  "#C3D7A4", 
                "#56B4E9", "#E69F00","#52854C", "#293352")
pdf(file = "output/figure3.pdf",height = 5.5, width = 7)
f2 <- ggplot(FigData_long2, aes(fill=Parameter, y=Value, x=Validation)) + 
  facet_grid(Type~nrate, scales = "free_y",switch ="y")+
  theme(text=element_text(size=13, family="mono"))+
  geom_bar(position="dodge", stat="identity") + 
  scale_fill_manual(values=custom.col) 


f2

dev.off()





#### 3. Figure 3 - Weight ####
betanames <- c("beta10","beta11","beta12","beta20","beta21","beta22")

S.E. <- c(0.046, 0.022, 0.045, 0.161, 0.133, 0.171, ## Validation data only
          0.045, 0.023, 0.049, 0.362, 0.480, 0.411, ## Non-validation data
          0.036, 0.018, 0.036, 0.215, 0.259, 0.245, ## Equal weight for validation data and non-validation data
          0.030, 0.014, 0.029, 0.140, 0.121, 0.147) ## optimal weighted
scenario <-c(rep("Scenario 1",6),rep("Scenario 2",6),rep("Scenario 3",6),rep("Scenario 4",6))
response<-c(rep(c(rep("Continuous",3),rep("Binary",3)),4))
betaname <- rep(betanames,4)

data_weight <- data.frame(SE=SE, Scenario=scenario, Parameter=betaname, Response = response)

data_weight$Response <- factor(data_weight$Response, 
                               levels = c("Continuous", "Binary"))

### 2.2  Plot the Figure ####
# Grouped
custom.col <- c("#FFDB6D",  "#C3D7A4", 
                "#56B4E9", "#E69F00","#52854C", "#293352")
pdf(file = "output/figure4.pdf",height = 5.5, width = 7)
f4 <- ggplot(data_weight, aes(fill=Parameter, y=S.E., x=Response)) + 
  facet_wrap(~ Scenario, ncol=2)+
  theme(text=element_text(size=13, family="mono"))+
  geom_bar(position="dodge", stat="identity") + 
  scale_fill_manual(values=custom.col) 


f4

dev.off()



#### Figure 2 Redo: Sample size ratio and weight ####

### Figure 4: Bias plot ####

load("WINV_RTable2.RData")
FigData2 <- data.frame(par=rep(parnames,3),response=rep(response0,3),Results,Val="Unweighted")
FigData2 <- FigData2[!FigData2$par=="",c("par","response","nrate", "biasprop", "sdpropose", "Val")]
load("WINV_RTable2_opt.RData")
FigData2f <- data.frame(par=rep(parnames,3),response=rep(response0,3),Results,Val="Weighted")
FigData2f <- FigData2f[!FigData2f$par=="",c("par","response","nrate", "biasprop", "sdpropose", "Val")]
load("WINV_RTable2_opt_GMM.RData")
FigData2e <- data.frame(par=rep(parnames,3),response=rep(response0,3),Results,Val="GMM")
FigData2e <- FigData2e[!FigData2e$par=="",c("par","response","nrate", "biasprop", "sdpropose", "Val")]

FigData2 <- rbind(FigData2,FigData2f,FigData2e)
colnames(FigData2) <- c("Parameter","Response","nrate","Bias","S.D.","Validation")

FigData_long2 <- gather(FigData2, Type , Value, c(Bias,S.D.))
FigData_long2$nrate <- factor(FigData_long2$nrate,levels=c("0.333","0.5","0.667"),labels=c("1:2","1:1","2:1"))
FigData_long2$Value <- abs(FigData_long2$Value)

FigData_long2$Validation <- factor(FigData_long2$Validation, levels=c("Unweighted","Weighted","GMM"))


### 2.2  Plot the Figure 1 ####
# Grouped
custom.col <- c("#DDB247", "#CD7E59", "#457373")
pdf(file = "output/figure16.pdf",height = 3, width = 9)
f4 <- ggplot(FigData_long2[FigData_long2$Type=="Bias",], aes(fill=Validation, y=Value, x=nrate)) + 
  facet_grid(.~Parameter, switch ="y",labeller= labeller(Parameter=label_parsed))+
  theme(text=element_text(size=13, family="mono"))+
  geom_bar(position="dodge", stat="identity") + 
  scale_fill_manual(values=custom.col) +
  labs(x = "Sample Size Ratio", y="Bias", fill="Method") + theme_bw() +
  theme(strip.background =element_rect(fill="#5A6868",color="#5A6868"))+
  theme(strip.text = element_text(colour = 'white'))+
  theme(panel.border = element_rect(colour = '#5A6868'))

f4

dev.off()



### 2.2  Plot the Figure 2 ####
# Grouped
custom.col <- c("#DDB247", "#CD7E59", "#457373")
pdf(file = "output/figure18.pdf",height = 3, width = 9)
f5 <- ggplot(FigData_long2[FigData_long2$Type=="S.D.",], aes(fill=Validation, y=Value, x=nrate)) + 
  facet_grid(.~Parameter, switch ="y",labeller= labeller(Parameter=label_parsed))+
  theme(text=element_text(size=13, family="mono"))+
  geom_bar(position="dodge", stat="identity") + 
  scale_fill_manual(values=custom.col) +
  labs(x = "Sample Size Ratio", y="S.E.", fill="Method") + theme_bw() +
  theme(strip.background =element_rect(fill="#5A6868",color="#5A6868"))+
  theme(strip.text = element_text(colour = 'white')) +
  theme(panel.border = element_rect(colour = '#5A6868'))

f5

dev.off()


