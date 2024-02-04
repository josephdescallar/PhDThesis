library(compriskMPL)
wbrt <- read.csv("...//wbrt.csv"). #Data available from Melanoma Institute Australia
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

#determine cancer related deaths
table(wbrt$Death_Reason)
wbrt$riskccr <- ifelse(is.na(wbrt$timeIntra) & wbrt$CancerDeath==0,NA, ifelse(is.na(wbrt$timeIntra) & wbrt$CancerDeath==1,2,1))
wbrt$eventccr <- ifelse(is.na(wbrt$riskccr), 0, ifelse(wbrt$riskccr==1, 3, 1))

#Kaplan-Meier curves
wbrt.cox1 <- coxph(Surv(tmid, riskccr_1==1) ~ factor(treatment.n) + factor(sex) + factor(Intracranial_Metastases_No) + factor(Extracranial_Dx) + factor(system_prior) + factor(system_within12) + factor(steroid2), data = wbrt.sub)
wbrt.cox2 <- coxph(Surv(tmid, riskccr_2=1) ~ factor(treatment.n) + factor(sex) + factor(Intracranial_Metastases_No) + factor(Extracranial_Dx) + factor(system_prior) + factor(system_within12) + factor(steroid2), data = wbrt.sub)

wbrt$yearsmid <- wbrt$tmid/365.25

wbrt.km1 <- survfit(Surv(yearsmid, riskccr_1==1) ~ treatment, data=wbrt)
plot(wbrt.km1, col=c("#003366", "#E31B23"), xlab="Years", ylab="Survival Probability")
legend(0, 0.15, legend=c("Control", "Treatment"),
       col=c("#003366", "#E31B23"), lty=1, cex=0.8, lwd=2)

wbrt.km2 <- survfit(Surv(yearsmid, riskccr_2==1) ~ treatment, data=wbrt)
plot(wbrt.km2, col=c("#003366", "#E31B23"), xlab="Years", ylab="Survival Probability")
legend(0, 0.15, legend=c("Control", "Treatment"),
       col=c("#003366", "#E31B23"), lty=1, cex=0.8, lwd=2)

wbrt.z.mat <- as.matrix(cbind(factor(wbrt.sub$treatment.n), factor(wbrt.sub$sex), factor(wbrt.sub$Intracranial_Metastases_No), factor(wbrt.sub$Extracranial_Dx), factor(wbrt.sub$system_prior), factor(wbrt.sub$system_within12),  factor(wbrt.sub$steroid2)))
colnames(wbrt.z.mat) <- c("treatment.n", "sex", "Intracranial_Metastases_No", "Extracranial_Dx", "system_prior", "system_within12", "steroid2")

wbrt.analysis.cure <- phcsh_mpl(Surv(time1, time2, event=eventccr, "interval") ~ factor(treatment.n) + factor(sex) + factor(Intracranial_Metastases_No) + factor(Extracranial_Dx) + factor(system_prior) + factor(system_within12)  + factor(steroid2),
                                risk = wbrt.sub$riskccr, data = wbrt.sub,
                                z = wbrt.z.mat,
                                max.outer = 10,
                                max.iter = 500, lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL,
                                aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = TRUE,
                                gq.points = 50)
summary.phcsh_mpl(wbrt.analysis.cure)

#####plots
#plot baseline hazards
wbrt.cure.1.cif.risk1 <- plot.phcsh_mpl(wbrt.analysis.cure1, risk=1)
wbrt.cure.1.cif.risk2 <- plot.phcsh_mpl(wbrt.analysis.cure1, risk=2)
max.bh <- max(c(wbrt.cure.1.cif.risk1$h0r.upper, wbrt.cure.1.cif.risk2$h0r.upper))
plot(wbrt.cure.1.cif.risk1$t.points[1:650] / 30.4375, wbrt.cure.1.cif.risk1$h0r[1:650], type = "l", col="#003366", ylim=c(0,max.bh), ylab="Baseline Hazard", xlab = "Months", lwd = 2, xlim = c(0,60))
lines(wbrt.cure.1.cif.risk1$t.points[1:650] / 30.4375, wbrt.cure.1.cif.risk1$h0r.lower[1:650], type = "l", lty="dashed", col="#003366", lwd = 2)
lines(wbrt.cure.1.cif.risk1$t.points[1:650] / 30.4375, wbrt.cure.1.cif.risk1$h0r.upper[1:650], type = "l", lty="dashed", col="#003366", lwd = 2)
lines(wbrt.cure.1.cif.risk2$t.points[1:650] / 30.4375, wbrt.cure.1.cif.risk2$h0r[1:650], type = "l", col="#E31B23", lwd = 2)
lines(wbrt.cure.1.cif.risk2$t.points[1:650] / 30.4375, wbrt.cure.1.cif.risk2$h0r.lower[1:650], type = "l", lty="dashed", col="#E31B23", lwd = 2)
lines(wbrt.cure.1.cif.risk2$t.points[1:650] / 30.4375, wbrt.cure.1.cif.risk2$h0r.upper[1:650], type = "l", lty="dashed", col="#E31B23", lwd = 2)
legend("topright", legend=c("Intracranial Failure", "Death"),
       col=c("#003366", "#E31B23"), lty=1, lwd=2, box.lty = 0)

#CIF

plot(wbrt.cure.1.cif.risk1$t.points[1:650] / 30.4375, wbrt.cure.1.cif.risk1$F0r[1:650], type = "l", col="#003366", ylim=c(0,1), ylab="Cumulative Incidence Function", xlab = "Months", lwd = 2, xlim = c(0,60))
lines(wbrt.cure.1.cif.risk1$t.points[1:650] / 30.4375, wbrt.cure.1.cif.risk1$F0r.lower[1:650], type = "l", lty="dashed", col="#003366", lwd = 2)
lines(wbrt.cure.1.cif.risk1$t.points[1:650] / 30.4375, wbrt.cure.1.cif.risk1$F0r.upper[1:650], type = "l", lty="dashed", col="#003366", lwd = 2)
lines(wbrt.cure.1.cif.risk2$t.points[1:650] / 30.4375, wbrt.cure.1.cif.risk2$F0r[1:650], type = "l", col="#E31B23", lwd = 2)
lines(wbrt.cure.1.cif.risk2$t.points[1:650] / 30.4375, wbrt.cure.1.cif.risk2$F0r.lower[1:650], type = "l", lty="dashed", col="#E31B23", lwd = 2)
lines(wbrt.cure.1.cif.risk2$t.points[1:650] / 30.4375, wbrt.cure.1.cif.risk2$F0r.upper[1:650], type = "l", lty="dashed", col="#E31B23", lwd = 2)
legend(0.001, 0.99, legend=c("Intracranial Failure", "Death"),
       col=c("#003366", "#E31B23"), lty=1, lwd=2, box.lty = 0)


#survival plots for treatment 0
par(mfrow=c(1,2))
plot(t.points[1:650] / 30.4375, treat0_surv$plot.S0r[[1]][1:650], ylim=c(0,1), type="l", col="#003366", lwd = 2, main="(a) Distant Intracranial failure", ylab="Probability", xlab="Months")
lines(t.points[1:650] / 30.4375, treat0_surv$plot.S0r.lower[[1]][1:650], col="#003366", lty=2)
lines(t.points[1:650] / 30.4375, treat0_surv$plot.S0r.upper[[1]][1:650], col="#003366", lty=2)

lines(t.points[1:650] / 30.4375, treat1_surv$plot.S0r[[1]][1:650], ylim=c(0,1), lwd=2, col="#E31B23")
lines(t.points[1:650] / 30.4375, treat1_surv$plot.S0r.lower[[1]][1:650], col="#E31B23", lty=2)
lines(t.points[1:650] / 30.4375, treat1_surv$plot.S0r.upper[[1]][1:650], col="#E31B23", lty=2)

legend(0, 0.15, legend=c("Control", "Treatment"),
       col=c("#003366", "#E31B23"), lty=1, cex=0.8, lwd=2)


plot(t.points[1:650] / 30.4375, treat0_surv$plot.S0r[[2]][1:650], ylim=c(0,1), type="l", col="#003366", lwd = 2, main="(b) Death", ylab="Probability", xlab="Months")
lines(t.points[1:650] / 30.4375, treat0_surv$plot.S0r.lower[[2]][1:650], col="#003366", lty=2)
lines(t.points[1:650] / 30.4375, treat0_surv$plot.S0r.upper[[2]][1:650], col="#003366", lty=2)

lines(t.points[1:650] / 30.4375, treat1_surv$plot.S0r[[2]][1:650], ylim=c(0,1), lwd=2, col="#E31B23")
lines(t.points[1:650] / 30.4375, treat1_surv$plot.S0r.lower[[2]][1:650], col="#E31B23", lty=2)
lines(t.points[1:650] / 30.4375, treat1_surv$plot.S0r.upper[[2]][1:650], col="#E31B23", lty=2)

legend(0, 0.15, legend=c("Control", "Treatment"),
       col=c("#003366", "#E31B23"), lty=1, cex=0.8, lwd=2)









