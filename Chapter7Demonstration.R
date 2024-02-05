library(devtools)
#install_github("josephdescallar/compriskMPL")
library(compriskMPL)

simdata_phcsh_mpl <- function(n.datasets = 1000, n.obs = 200, prob_event = 0.2, 
                              rho = 3, lambda1 = 1, lambda2 = 0.5, 
                              betas = list(c(0.5,-2), c(1,-1)), gammaL = 0.5, 
                              gammaR = 2.5){
  simdata <- list()
  for(sim in 1:n.datasets){
    seed <- sim
    set.seed(seed)
    #Generate U~unif(0,1)
    u <- runif(n.obs)
    #Generate x1 and x2
    x1 <- rnorm(n.obs, mean=0, sd=1)
    x2 <- rnorm(n.obs, mean=0, sd=1)
    #x2 <- runif(n.obs, 0, 2)
    x <- data.matrix(data.frame(x1,x2))
    #Calculate t random variable
    t <- (-log(u) / (lambda1*exp(x %*% betas[[1]]) + 
                       lambda2*(exp(x %*% betas[[2]]))))^(1/rho)
    Ul <- runif(n.obs)
    Ur <- runif(n.obs,Ul,1)
    Ue <- runif(n.obs)
    event <- ifelse(Ue < prob_event,1,ifelse(gammaL*Ul <= t & t <= gammaR*Ur,
                                             3,ifelse(t < gammaL*Ul,2,0)))
    #Generate V to determine which risk t belongs
    p1 <- (lambda1*exp(x %*% betas[[1]])) / ((lambda1*exp(x %*% betas[[1]]) + 
                                                lambda2*(exp(x %*% betas[[2]]))))
    p2 <- (lambda2*exp(x %*% betas[[2]])) / ((lambda1*exp(x %*% betas[[1]]) + 
                                                lambda2*(exp(x %*% betas[[2]]))))
    v <- runif(n.obs)
    time1 <- ifelse(event==0,gammaR*Ur,ifelse(event==1,t,ifelse(event==2,
                                                                gammaL*Ul,gammaL*Ul)))
    time2 <- ifelse(event==0,NA,ifelse(event==1,t,ifelse(event==2,NA,
                                                         gammaR*Ur)))
    tmid <- ifelse(event==3,rowMeans(cbind(time1,time2)),ifelse(event==2,
                                                                time1/2, time1))
    risk <- ifelse(event==0,NA,ifelse(v <= p1, 1, 2))
    risk_1 <- ifelse((event==1 | event==2 | event==3) & risk == 1, 1, 0)
    risk_2 <- ifelse((event==1 | event==2 | event==3) & risk == 2, 1, 0)
    risk2 <- ifelse(is.na(risk),0,risk)
    x <- data.frame(x1,x2)
    simdata[[sim]] <- data.frame("time"=time1,time2,x,event,risk)
  }
  simdata
}

# Section 7.4
sim_data <- simdata_phcsh_mpl(n.datasets=1)
sim_data
sim_data_phcsh <- phcsh_mpl(Surv(time, time2, event=event, "interval") ~ x1 + 
                              x2, data=sim_data[[1]], risk=sim_data[[1]]$risk)
summary.phcsh_mpl(sim_data_phcsh)
summary.phcsh_mpl(sim_data_phcsh, grad=TRUE)

#predicitons
predict.phcsh_mpl(sim_data_phcsh, risk=1)
predict.phcsh_mpl(sim_data_phcsh, risk=2)

#plots for risk 1
plot.phcsh_mpl(sim_data_phcsh, risk = 1)
plot.phcsh_mpl(sim_data_phcsh, plots="bh", risk = 1)
plot.phcsh_mpl(sim_data_phcsh, plots="surv", risk = 1)
plot.phcsh_mpl(sim_data_phcsh, plots="cif", risk = 1)

#plots for risk 2
plot.phcsh_mpl(sim_data_phcsh, risk = 2)
plot.phcsh_mpl(sim_data_phcsh, plots="bh", risk = 2)
plot.phcsh_mpl(sim_data_phcsh, plots="surv", risk = 2)
plot.phcsh_mpl(sim_data_phcsh, plots="cif", risk = 2)

# Melanoma cancer data example
wbrt <- read.csv("...//wbrt.csv"). #Data available from Melanoma Institute Australia

#Some data setup
#Data setup
wbrt$timeLastMRI <- wbrt$Last_MRI_Date_Num - wbrt$startDateNum
wbrt$timeLastMRI[wbrt$timeLastMRI < 0] <- 0
wbrt$timeIntra <- wbrt$Fail_MRI_Date_Num - wbrt$startDateNum
wbrt$timeDeath <- wbrt$Date_Death_Num - wbrt$startDateNum
wbrt$timeLast <- wbrt$Last_Date_Num - wbrt$startDateNum
wbrt$CancerDeath <- ifelse(wbrt$Death_Reason=='Intracranial progression of disease' | wbrt$Death_Reason=='Other melanoma specific cause',1,0)
wbrt$risk <- ifelse(is.na(wbrt$timeIntra) & is.na(wbrt$timeDeath), NA,
                    ifelse(is.na(wbrt$timeIntra) & !is.na(wbrt$timeDeath), 2,
                           ifelse(!is.na(wbrt$timeIntra) & is.na(wbrt$timeDeath), 1,
                                  1 )))
wbrt$event <- ifelse(is.na(wbrt$risk),0,
                     ifelse(wbrt$risk==1 & wbrt$timeLastMRI != 0 ,3,
                            ifelse(wbrt$risk==1, 3,1)))
wbrt$time1 <- ifelse(wbrt$event==0, wbrt$timeLast,
                     ifelse(wbrt$event==1, wbrt$timeDeath, wbrt$timeLastMRI))
wbrt$time2 <- ifelse(wbrt$event==0, NA,
                     ifelse(wbrt$event==1, wbrt$timeDeath, wbrt$timeIntra))
wbrt$treatment.n <- ifelse(wbrt$treatment=='Observation',1,2)
wbrt$tmid <- ifelse(wbrt$event==3,rowMeans(cbind(wbrt$time1,wbrt$time2)),ifelse(wbrt$event==2,wbrt$time1/2, wbrt$time1))
wbrt$risk_1 <- ifelse((wbrt$event==1 | wbrt$event==2 | wbrt$event==3) & wbrt$risk == 1, 1, 0)
wbrt$risk_2 <- ifelse((wbrt$event==1 | wbrt$event==2 | wbrt$event==3) & wbrt$risk == 2, 1, 0)
wbrt$treatment.n <- as.factor(wbrt$treatment.n)
wbrt$sex <- as.factor(wbrt$sex)
wbrt$Intracranial_Metastases_No <- as.factor(wbrt$Intracranial_Metastases_No)
wbrt$Extracranial_Dx <- as.factor(wbrt$Extracranial_Dx)
wbrt$systemic <- as.factor(wbrt$systemic)
wbrt$steroid <- as.factor(wbrt$steroid)

#Section 7.5
table("Risk"=wbrt$risk, "Event"=wbrt$event)
table("Treatment Group" = wbrt$treatment)

#Specify the z matrix for the cure fraction
wbrt.z <- as.matrix(cbind(factor(wbrt$treatment.n), 
                          factor(wbrt$sex), factor(wbrt$Intracranial_Metastases_No), 
                          factor(wbrt$Extracranial_Dx), factor(wbrt$system_prior), 
                          factor(wbrt$system_within12),  factor(wbrt$steroid2)))

#provide column names for z cure fraction
colnames(wbrt.z) <- c("treatment.n", "sex", "Intracranial_Metastases_No", 
                      "Extracranial_Dx", "system_prior", "system_within12", "steroid")

#fit wbrt function with a cure fraction
wbrt.analysis <- phcsh_mpl(Surv(time1, time2, event=eventccr, "interval") ~ 
                             treatment + factor(sex) + factor(Intracranial_Metastases_No) + 
                             factor(Extracranial_Dx) + factor(system_prior) + factor(system_within12)  
                           + steroid2, risk = wbrt$riskccr, data = wbrt,
                           z = wbrt.z)

#summary information for the fit model
summary.phcsh_mpl(wbrt.analysis)
summary.phcsh_mpl(wbrt.analysis, grad=TRUE)$risks[[1]]$Theta
summary.phcsh_mpl(wbrt.analysis, grad=TRUE)$risks[[2]]$Theta

# Section 7.5.1
#predictions for reference levels
wbrt.predict.r1 <- predict.phcsh_mpl(wbrt.analysis)

# View first 10 observations from predictions for baseline hazard, survival and CIF
wbrt.predict.r1[1:10,c("t", "pred.h0r", "pred.h0r.lower", "pred.h0r.upper")]
wbrt.predict.r1[1:10,c("t", "pred.S0r", "pred.S0r.lower", "pred.S0r.upper")]
wbrt.predict.r1[1:10,c("t", "pred.F0r", "pred.F0r.lower", "pred.F0r.upper")]

#plots of hazard, survival and CIF for distant intracranial disease
plot.phcsh_mpl(wbrt.analysis, pred=wbrt.predict.r1, plots="bh")
plot.phcsh_mpl(wbrt.analysis, pred=wbrt.predict.r1, plots="surv")
plot.phcsh_mpl(wbrt.analysis, pred=wbrt.predict.r1, plots="cif")

#prediction values for cancer specific death
wbrt.predict.r2 <- predict.phcsh_mpl(wbrt.analysis, risk = 2)

#plots of hazard, survival and CIF for cancer specific death
plot.phcsh_mpl(wbrt.analysis, pred=wbrt.predict.r2, plots="bh", risk = 2)
plot.phcsh_mpl(wbrt.analysis, pred=wbrt.predict.r2, plots="surv", risk = 2)
plot.phcsh_mpl(wbrt.analysis, pred=wbrt.predict.r2, plots="cif", risk = 2)



# setup covariate vector for treatment group, whilst keeping others
# at their default reference levels
int.covs <- c(1,0,0,0,0,0,0)

#predictions for wbrt treatment groups
wbrt.predict.r1.int <- predict.phcsh_mpl(wbrt.analysis, 
                                         covs = int.covs)
wbrt.predict.r2.int <- predict.phcsh_mpl(wbrt.analysis, risk = 2, 
                                         covs = int.covs)

plot.phcsh_mpl(wbrt.analysis, pred=wbrt.predict.r1.int, plots="bh")
plot.phcsh_mpl(wbrt.analysis, pred=wbrt.predict.r1.int, plots="surv")
plot.phcsh_mpl(wbrt.analysis, pred=wbrt.predict.r1.int, plots="cif")

plot.phcsh_mpl(wbrt.analysis, pred=wbrt.predict.r2.int, plots="bh", risk = 2)
plot.phcsh_mpl(wbrt.analysis, pred=wbrt.predict.r2.int, plots="surv", risk = 2)
plot.phcsh_mpl(wbrt.analysis, pred=wbrt.predict.r2.int, plots="cif", risk = 2)


