library(devtools)
#install_github("compriskMPL")
library(compriskMPL)
library(parallel)
library(MASS)
library(foreach)
library(doParallel)
library(doSNOW)
library(gbm)
library(Hmisc)

#Function to generate data for simulation
gen.phcsh.data <- function(n.datasets = 1000, n.obs, prob_event, rho, lambda1, 
                          lambda2, betas, gammaL, gammaR){
  simdata <- list()
  for(sim in 1:n.datasets){
    b11 <- betas[[1]][1]
    b12 <- betas[[1]][2]
    b21 <- betas[[2]][1]
    b22 <- betas[[2]][2]
    seed <- sim
    set.seed(seed)
    #Generate U~unif(0,1)
    u <- runif(n.obs)
    #Generate x1 and x2
    x1 <- rnorm(n.obs, mean=0, sd=1)
    x2 <- rnorm(n.obs, mean=0, sd=1)
    x <- data.matrix(data.frame(x1,x2))
    #Calculate t random variable
    t <- (-log(u) / (lambda1*exp(x1*b11 + x2*b12) + 
                       lambda2*(exp(x1*b21 + x2*b22))))^(1/rho)
    Ul <- runif(n.obs)
    Ur <- runif(n.obs,Ul,1)
    Ue <- runif(n.obs)
    event <- ifelse(Ue < prob_event,1,ifelse(gammaL*Ul <= t & t <= gammaR*Ur,3,
                                             ifelse(t < gammaL*Ul,2,0)))
    #Generate V to determine which risk t belongs
    p1 <- (lambda1*exp(x1*b11 + x2*b12)) / ((lambda1*exp(x1*b11 + x2*b12) + 
                                               lambda2*(exp(x1*b21 + x2*b22))))
    p2 <- (lambda2*exp(x1*b21 + x2*b22)) / ((lambda1*exp(x1*b11 + x2*b12) + 
                                               lambda2*(exp(x1*b21 + x2*b22))))
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
    simdata[[sim]] <- data.frame("time"=time1,time2,x,event,risk,tmid,risk2,t,
                                 risk_1,risk_2)
  }
  simdata
}

#Function to display bias, std, cp for cox regression results
sim.results.cox <- function(object, beta, simdata){
  betalist = seBlist = biaslist = mselist = lowerlist = upperlist = list()
  for(sim in 1:simdata){
    betalist[[sim]] = object[[sim]]$coefficients
    biaslist[[sim]] = beta - betalist[[sim]]
    mselist[[sim]] = biaslist[[sim]]^2
    seBlist[[sim]] = sqrt(diag(object[[sim]]$var))
    lowerlist[[sim]] = betalist[[sim]] - 1.96*seBlist[[sim]]
    upperlist[[sim]] = betalist[[sim]] + 1.96*seBlist[[sim]]
  }
  lower = Reduce("rbind", lowerlist)
  upper = Reduce("rbind", upperlist)
  betamat = matrix(beta, nrow = simdata, ncol = length(beta), byrow = TRUE)
  fit = list()
  fit$bias = colMeans(Reduce("rbind",biaslist))
  fit$seB = colMeans(Reduce("rbind",seBlist))
  fit$seBmc = apply(Reduce("rbind",betalist), 2, sd)
  fit$mse = colMeans(Reduce("rbind",mselist))
  fit$cov.prob = colMeans(ifelse(betamat >= lower & betamat <= upper,1,0))
  fit
}

#Function to display bias, std, cp for MPL regression
sim.results.mpl <- function(object, beta, simdata){
  biaslist = valid = seBlist = betalist =  mselist = list()
  lowerlist = upperlist = seBsandlist = lowersandlist = list()
  uppersandlist = list()
  for(sim in 1:simdata){
    valid[[sim]] = object[sim,]$valid
    betalist[[sim]] = unlist(object[sim, ]$"beta")
    biaslist[[sim]] = unlist(beta) - betalist[[sim]]
    seBlist[[sim]] = unlist(object[sim, ]$"seB")
    seBsandlist[[sim]] = unlist(object[sim, ]$"seB.sand")
    mselist[[sim]] = biaslist[[sim]]^2
    lowerlist[[sim]] = betalist[[sim]] - 1.96*seBlist[[sim]]
    upperlist[[sim]] = betalist[[sim]] + 1.96*seBlist[[sim]]
    lowersandlist[[sim]] = betalist[[sim]] - 1.96*seBsandlist[[sim]]
    uppersandlist[[sim]] = betalist[[sim]] + 1.96*seBsandlist[[sim]]
  }
  valid.temp = unlist(valid)
  valid.index = which(valid.temp %in% 1)
  count=1
  biaslist.valid = seBlist.valid = betalist.valid = fit = mselist.valid = list()
  lowerlist.valid = upperlist.valid = seBsandlist.valid = list()
  lowersandlist.valid = uppersandlist.valid = list()
  for(sim in 1:simdata){
    if(valid.temp[sim]==1){
      biaslist.valid[[count]] = biaslist[[sim]]
      seBlist.valid[[count]] = seBlist[[sim]]
      betalist.valid[[count]] = betalist[[sim]]
      mselist.valid[[count]] = mselist[[sim]]
      lowerlist.valid[[count]] = lowerlist[[sim]]
      upperlist.valid[[count]] = upperlist[[sim]]
      seBsandlist.valid[[count]] = seBsandlist[[sim]]
      lowersandlist.valid[[count]] = lowersandlist[[sim]]
      uppersandlist.valid[[count]] = uppersandlist[[sim]]
      count = count + 1
    }
  }
  lower = Reduce("rbind", lowerlist.valid)
  upper = Reduce("rbind", upperlist.valid)
  betamat = matrix(unlist(beta), nrow = count-1, ncol = length(unlist(beta)), 
                   byrow = TRUE)
  lowersand = Reduce("rbind", lowersandlist.valid)
  uppersand = Reduce("rbind", uppersandlist.valid)
  fit$valid = mean(valid.temp)
  fit$bias = colMeans(Reduce("rbind", biaslist.valid))
  fit$seB.asymp = colMeans(Reduce("rbind", seBlist.valid))
  fit$seB.sand = colMeans(Reduce("rbind", seBsandlist.valid))
  fit$seB.mc = apply(Reduce("rbind", betalist.valid), 2, sd)
  fit$mse = colMeans(Reduce("rbind", mselist.valid))
  fit$cov.prob = colMeans(ifelse(betamat >= lower & betamat <= upper,1,0))
  fit$cov.prob.sand = colMeans(ifelse(betamat >= lowersand & betamat <= uppersand,1,0))
  fit
}

#Function to predict baseline hazard from MPL
pred.bh <- function(object, t.points, r, sand = FALSE, rho=0, lambda=0, 
                    var.max=99999, coef){
  psif <- function(x, bknots, iknots){
    mSpline(x, knots = iknots, Boundary = bknots, degree = object$dgr,
            intercept = object$basis.intercept)}
  pred.psi = psif(t.points, object$b.knots, object$i.knots[[r]])
  pred.h0r = as.vector(pred.psi %*% object$"theta"[[r]]) 
  if(sand == FALSE){
    VarCovMat.theta = object$VarCovMat[object$theta.index[[r]],
                                       object$theta.index[[r]]]
  }
  else if(sand == TRUE){ #up to here
    VarCovMat.theta = object$sand[object$theta.index[[r]], object$theta.index[[r]]]
  }
  pred.h0r.var = diag(pred.psi %*% VarCovMat.theta %*% t(pred.psi))
  if(object$pos.def==1){
    pred.var = pred.h0r.var / pred.h0r^2
    pred.var[pred.var > var.max] <- var.max
    pred.h0r.lower = as.vector(pred.h0r) * exp(-1.96*sqrt(pred.var))
    pred.h0r.upper = as.vector(pred.h0r) * exp(1.96*sqrt(pred.var))
  }
  true.h0r <- lambda*rho*t.points^(rho-1)
  if(object$valid==1){
    rlist <- list("pred.h0r"=pred.h0r, "pred.h0r.lower"=pred.h0r.lower,
                  "pred.h0r.upper"=pred.h0r.upper, "true.h0r"=true.h0r)
  }
  else{
    rlist <- list("pred.h0r"=pred.h0r, "true.h0r"=true.h0r)
  }
  rlist
}

#Function to predic CIF from MPL
pred.CIF <- function(object, t.points, sand=FALSE, rho, lambda, var.max, r){
  psif <- function(x, bknots, iknots){
    mSpline(x, knots = iknots, Boundary = bknots, degree = object$dgr,
            intercept = object$basis.intercept)}
  PSIf <- function(x, bknots, iknots)
    iSpline(x, knots = iknots, Boundary = bknots, degree = object$dgr,
            intercept = object$basis.intercept)
  bma = (t.points - object$b.knots[1]) / 2
  bpa = (t.points + object$b.knots[1]) / 2 #up to here 28
  t.gq.change = bma*matrix(object$nodes, nrow=length(bma), 
      ncol = object$gq.points, byrow=TRUE) + bpa
  pred.F0r.psi.gq = pred.F0r.PSI.gq = pred.F0r.h0qt.gq = list()
  pred.F0r.H0qt.gq = list()
  pred.F0r.S0qt.gq = pred.F0r.Integrand.gq = pred.F0r.dhdT.gq = list() 
  pred.F0r.dSdT.gq = pred.F0r.dFdT.gqT = list()
  for(gq in 1:object$gq.points){
    pred.F0r.psi.gq[[gq]] <- psif(t.gq.change[,gq], object$b.knots,
                                  object$i.knots[[r]])
    pred.F0r.PSI.gq[[gq]] <- PSIf(t.gq.change[,gq], object$b.knots,
                                  object$i.knots[[r]])
    pred.F0r.h0qt.gq[[gq]] <- psif(t.gq.change[,gq], object$b.knots,
                                   object$i.knots[[r]]) %*% object$"theta"[[r]]
    pred.F0r.H0qt.gq.r = pred.F0r.S0qt.gq.r = list()
    for(q in 1:object$n.risk){
      pred.F0r.H0qt.gq.r[[q]] = PSIf(t.gq.change[,gq], object$b.knots,
                                    object$i.knots[[q]]) %*% object$"theta"[[q]]
      pred.F0r.S0qt.gq.r[[q]] = exp(-pred.F0r.H0qt.gq.r[[q]])
    }
    pred.F0r.S0qt.gq[[gq]] = Reduce("*", pred.F0r.S0qt.gq.r)
    pred.F0r.Integrand.gq[[gq]] <- object$weights[gq]*pred.F0r.h0qt.gq[[gq]] *
      pred.F0r.S0qt.gq[[gq]]
    pred.F0r.dhdT.gq[[gq]] <- pred.F0r.psi.gq[[gq]]
    pred.F0r.dSdT.gq[[gq]] <- -pred.F0r.PSI.gq[[gq]] *
      as.vector(pred.F0r.S0qt.gq[[gq]])
    pred.F0r.dFdT.gqT[[gq]] <- (as.matrix(pred.F0r.dhdT.gq[[gq]] *
                                as.vector(pred.F0r.S0qt.gq[[gq]]) +
                                as.vector(pred.F0r.h0qt.gq[[gq]]) *
                                as.matrix(pred.F0r.dSdT.gq[[gq]]))) *
      object$weights[gq]
  }
  
  pred.F0r <- as.vector(bma*Reduce("+",pred.F0r.Integrand.gq))
  pred.F0r.dFdT <- bma*Reduce("+",pred.F0r.dFdT.gqT)
  if(sand == FALSE){
    VarCovMat.theta = object$VarCovMat[object$theta.index[[r]],
                                       object$theta.index[[r]]]
  }
  else if(sand == TRUE){
  VarCovMat.theta = object$sand[object$theta.index[[r]], object$theta.index[[r]]]
  }
  pred.F0r.r.var <- diag(pred.F0r.dFdT %*% as.matrix(VarCovMat.theta) %*%
                           t(pred.F0r.dFdT))
  pred.F0r.logOR <- log(pred.F0r / (1-pred.F0r + 1e-12) + 1e-12)
  pred.F0r.logOR.var = ((1/((1-pred.F0r)*pred.F0r))^2)*pred.F0r.r.var
  pred.F0r.log <- log(pred.F0r + 1e-12)
  pred.F0r.log.var <- (1/pred.F0r^2)*pred.F0r.r.var
  if(object$pos.def==1){
    pred.F0r.logOR.lower = pred.F0r.logOR - 1.96*sqrt(pred.F0r.logOR.var)
    pred.F0r.logOR.upper = pred.F0r.logOR + 1.96*sqrt(pred.F0r.logOR.var)
    pred.F0r.lower <- as.vector(exp(pred.F0r.logOR.lower) / 
                      (1 + exp(pred.F0r.logOR.lower)))
    pred.F0r.upper <- as.vector(exp(pred.F0r.logOR.upper) / 
                    (1 + exp(pred.F0r.logOR.upper)))
    
  }
  true.h0r <- lambda*rho*t.points^(rho-1)
  true.H0r <- lambda*(t.points^rho)
  true.S0r <- exp(-true.H0r)
  if(object$valid==1){
    rlist <- list("pred.F0r"=pred.F0r, "pred.F0r.lower"=pred.F0r.lower,
                  "pred.F0r.upper"=pred.F0r.upper)
  }
  else{
    rlist <- list("pred.F0r"=pred.F0r)
  }
  rlist
}

#Study 1, n = 200, left/interval = 47.5%, right cens = 47.5%
n.sim=1000
true.beta = list(c(-1,0.5), c(1,-0.5))
ic.scen1.r47.n200 <- gen.phcsh.data(n.datasets = 1000, n.obs = 200,
prob_event = 0.05, rho = 3, lambda1 = 1, lambda2 = 0.5,
betas = true.beta, gammaL = 0.5, gammaR = 0.89)
table(unlist(lapply(ic.scen1.r47.n200, function(a) a$event))) /
length(unlist(lapply(ic.scen1.r47.n200, function(a) a$event)))
#cox midpoint t with n = 200
ic.scen1.r47.n200.r1 = ic.scen1.r47.n200.r2 = list()
for(sim in 1:n.sim){
  ic.scen1.r47.n200.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, 
  data = ic.scen1.r47.n200[[sim]])
  ic.scen1.r47.n200.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, 
  data = ic.scen1.r47.n200[[sim]])
}
sim.results.cox(ic.scen1.r47.n200.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(ic.scen1.r47.n200.r2, beta = true.beta[[2]], simdata = n.sim)

cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
ic.scen1.r47.n200.time <- system.time({ic.scen1.r47.n200.results <- foreach(i = 1:n.sim, 
.combine = rbind, .options.snow = opts) %dopar% { phcsh_mpl(Surv(time, time2, 
event=event, "interval") ~ x1 + x2, risk = ic.scen1.r47.n200[[i]]$risk, 
data = ic.scen1.r47.n200[[i]], max.outer = 10, max.iter = 3000, 
lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL, aps = TRUE,
knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = FALSE)
}
})
close(pb)
stopCluster(cl)
sim.results.mpl(ic.scen1.r47.n200.results, beta = true.beta, simdata = n.sim)


#Study 1, n = 200, left/interval = 75%, right cens = 20%
true.beta = list(c(-1,0.5), c(1,-0.5))
ic.scen1.r20.n200 <- gen.phcsh.data(n.datasets = 1000, n.obs = 200,
prob_event = 0.05, rho = 3, lambda1 = 1, lambda2 = 0.5,
betas = true.beta, gammaL = 0.5, gammaR = 1.34)
table(unlist(lapply(ic.scen1.r20.n200, function(a) a$event))) /
  length(unlist(lapply(ic.scen1.r20.n200, function(a) a$event)))
#cox midpoint t with n = 200
ic.scen1.r20.n200.r1 = ic.scen1.r20.n200.r2 = list()
for(sim in 1:n.sim){
  ic.scen1.r20.n200.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, data = ic.scen1.r20.n200[[sim]])
  ic.scen1.r20.n200.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, data = ic.scen1.r20.n200[[sim]])
}
sim.results.cox(ic.scen1.r20.n200.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(ic.scen1.r20.n200.r2, beta = true.beta[[2]], simdata = n.sim)

cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
ic.scen1.r20.n200.time <- system.time({ic.scen1.r20.n200.results <- foreach(i = 1:n.sim, 
  .combine = rbind, .options.snow = opts) %dopar% { phcsh_mpl(Surv(time, time2, 
  event=event, "interval") ~ x1 + x2, risk = ic.scen1.r20.n200[[i]]$risk, #up to here
  data = ic.scen1.r20.n200[[i]], max.outer = 10, max.iter = 3000, 
  lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL, aps = TRUE,
  knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = FALSE)
}
})
close(pb)
stopCluster(cl)
sim.results.mpl(ic.scen1.r20.n200.results, beta = true.beta, simdata = n.sim)


#Study 1, n = 1000, left/interval = 47.5%, right cens = 47.5%
n.sim=1000
true.beta = list(c(-1,0.5), c(1,-0.5))
ic.scen1.r47.n1000 <- gen.phcsh.data(n.datasets = 1000, n.obs = 1000,
prob_event = 0.05, rho = 3, lambda1 = 1, lambda2 = 0.5,
betas = true.beta, gammaL = 0.5, gammaR = 0.89)
table(unlist(lapply(ic.scen1.r47.n1000, function(a) a$event))) /
  length(unlist(lapply(ic.scen1.r47.n1000, function(a) a$event)))
#cox midpoint t with n = 1000
ic.scen1.r47.n1000.r1 = ic.scen1.r47.n1000.r2 = list()
for(sim in 1:n.sim){
  ic.scen1.r47.n1000.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, 
  data = ic.scen1.r47.n1000[[sim]])
  ic.scen1.r47.n1000.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, 
  data = ic.scen1.r47.n1000[[sim]])
}
sim.results.cox(ic.scen1.r47.n1000.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(ic.scen1.r47.n1000.r2, beta = true.beta[[2]], simdata = n.sim)

cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
ic.scen1.r47.n1000.time <- system.time({ic.scen1.r47.n1000.results <- foreach(i = 1:n.sim, 
  .combine = rbind, .options.snow = opts) %dopar% {
  phcsh_mpl(Surv(time, time2, event=event, "interval") ~ x1 + x2, 
  risk = ic.scen1.r47.n1000[[i]]$risk, data = ic.scen1.r47.n1000[[i]], 
  max.outer = 10, max.iter = 3000, lambda = c(1e-5,1e-5), n.basis = NULL, 
  iknots.pos = NULL, aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, 
  iter.disp = FALSE)
}
})
close(pb)
stopCluster(cl)
sim.results.mpl(ic.scen1.r47.n1000.results, beta = true.beta, simdata = n.sim)


#Study 1, n = 1000, left/interval = 75%, right cens = 20%
true.beta = list(c(-1,0.5), c(1,-0.5))
ic.scen1.r20.n1000 <- gen.phcsh.data(n.datasets = 1000, n.obs = 1000,
prob_event = 0.05, rho = 3, lambda1 = 1, lambda2 = 0.5,
betas = true.beta, gammaL = 0.5, gammaR = 1.34)
table(unlist(lapply(ic.scen1.r20.n1000, function(a) a$event))) /
  length(unlist(lapply(ic.scen1.r20.n1000, function(a) a$event)))
#cox midpoint t with n = 200
ic.scen1.r20.n1000.r1 = ic.scen1.r20.n1000.r2 = list()
for(sim in 1:n.sim){
  ic.scen1.r20.n1000.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, data = ic.scen1.r20.n1000[[sim]])
  ic.scen1.r20.n1000.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, data = ic.scen1.r20.n1000[[sim]])
}
sim.results.cox(ic.scen1.r20.n1000.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(ic.scen1.r20.n1000.r2, beta = true.beta[[2]], simdata = n.sim)

cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
ic.scen1.r20.n1000.time <- system.time({ic.scen1.r20.n1000.results <- foreach(i = 1:n.sim, 
  .combine = rbind, .options.snow = opts) %dopar% {phcsh_mpl(Surv(time, time2, 
  event=event, "interval") ~ x1 + x2, risk = ic.scen1.r20.n1000[[i]]$risk, 
  data = ic.scen1.r20.n1000[[i]], max.outer = 10, max.iter = 3000,
  lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL, aps = TRUE,
  knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = FALSE)
}
})
close(pb)
stopCluster(cl)
sim.results.mpl(ic.scen1.r20.n1000.results, beta = true.beta, simdata = n.sim)


#Study 2, n = 200, left/interval = 47.5%, right cens = 47.5%
true.beta = list(c(1,0.5), c(0.5,0.5))
ic.scen2.r47.n200 <- gen.phcsh.data(n.datasets = 1000, n.obs = 200,
prob_event = 0.05, rho = 3, lambda1 = 1, lambda2 = 0.5,
betas = true.beta, gammaL = 0.5, gammaR = 1.)
table(unlist(lapply(ic.scen2.r47.n200, function(a) a$event))) /
  length(unlist(lapply(ic.scen2.r47.n200, function(a) a$event)))
#cox midpoint t with n = 200
ic.scen2.r47.n200.r1 = ic.scen2.r47.n200.r2 = list()
for(sim in 1:n.sim){
  ic.scen2.r47.n200.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, data = ic.scen2.r47.n200[[sim]])
  ic.scen2.r47.n200.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, data = ic.scen2.r47.n200[[sim]])
}
sim.results.cox(ic.scen2.r47.n200.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(ic.scen2.r47.n200.r2, beta = true.beta[[2]], simdata = n.sim)
cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
ic.scen2.r47.n200.time <- system.time({ic.scen2.r47.n200.results <- foreach(i = 1:n.sim, 
  .combine = rbind, .options.snow = opts) %dopar% {
  phcsh_mpl(Surv(time, time2, event=event, "interval") ~ x1 + x2,
  risk = ic.scen2.r47.n200[[i]]$risk, data = ic.scen2.r47.n200[[i]], 
  max.outer = 10, max.iter = 3000, lambda = c(1e-5,1e-5), n.basis = NULL, 
  iknots.pos = NULL, aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, 
  iter.disp = FALSE)
}
})
close(pb)
stopCluster(cl)
sim.results.mpl(ic.scen2.r47.n200.results, beta = true.beta, simdata = n.sim)


#Study 2, n = 200, left/interval = 75%, right cens = 20%
true.beta = list(c(1,0.5), c(0.5,0.5))
ic.scen2.r20.n200 <- gen.phcsh.data(n.datasets = 1000, n.obs = 200,
prob_event = 0.05, rho = 3, lambda1 = 1, lambda2 = 0.5,
betas = true.beta, gammaL = 0.5, gammaR = 1.64)
table(unlist(lapply(ic.scen2.r20.n200, function(a) a$event))) /
  length(unlist(lapply(ic.scen2.r20.n200, function(a) a$event)))
#cox midpoint t with n = 200
ic.scen2.r20.n200.r1 = ic.scen2.r20.n200.r2 = list()
for(sim in 1:n.sim){
  ic.scen2.r20.n200.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, 
  data = ic.scen2.r20.n200[[sim]])
  ic.scen2.r20.n200.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, 
  data = ic.scen2.r20.n200[[sim]])
}
sim.results.cox(ic.scen2.r20.n200.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(ic.scen2.r20.n200.r2, beta = true.beta[[2]], simdata = n.sim)

cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
ic.scen2.r20.n200.time <- system.time({ic.scen2.r20.n200.results <- foreach(i = 1:n.sim, 
  .combine = rbind, .options.snow = opts) %dopar% {
  phcsh_mpl(Surv(time, time2, event=event, "interval") ~ x1 + x2,
  risk = ic.scen2.r20.n200[[i]]$risk, data = ic.scen2.r20.n200[[i]], 
  max.outer = 10, max.iter = 3000, lambda = c(1e-5,1e-5), n.basis = NULL, 
  iknots.pos = NULL, aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, 
  iter.disp = FALSE)
}
})
close(pb)
stopCluster(cl)
sim.results.mpl(ic.scen2.r20.n200.results, beta = true.beta, simdata = n.sim)


#Study 2, n = 1000, left/interval = 47.5%, right cens = 47.5%
true.beta = list(c(1,0.5), c(0.5,0.5))
ic.scen2.r47.n1000 <- gen.phcsh.data(n.datasets = 1000, n.obs = 1000,
prob_event = 0.05, rho = 3, lambda1 = 1, lambda2 = 0.5,
betas = true.beta, gammaL = 0.5, gammaR = 1.)
table(unlist(lapply(ic.scen2.r47.n1000, function(a) a$event))) /
  length(unlist(lapply(ic.scen2.r47.n1000, function(a) a$event)))
#cox midpoint t with n = 1000
ic.scen2.r47.n1000.r1 = ic.scen2.r47.n1000.r2 = list()
for(sim in 1:n.sim){
  ic.scen2.r47.n1000.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, 
  data = ic.scen2.r47.n1000[[sim]])
  ic.scen2.r47.n1000.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, 
  data = ic.scen2.r47.n1000[[sim]])
}
sim.results.cox(ic.scen2.r47.n1000.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(ic.scen2.r47.n1000.r2, beta = true.beta[[2]], simdata = n.sim)

cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
ic.scen2.r47.n1000.time <- system.time({ic.scen2.r47.n1000.results <- foreach(i = 1:n.sim, 
  .combine = rbind, .options.snow = opts) %dopar% {
  phcsh_mpl(Surv(time, time2, event=event, "interval") ~ x1 + x2,
  risk = ic.scen2.r47.n1000[[i]]$risk, data = ic.scen2.r47.n1000[[i]],
  max.outer = 10, max.iter = 3000, lambda = c(1e-5,1e-5), n.basis = NULL, 
  iknots.pos = NULL, aps = TRUE, knots.perc.limit = c(0.05,0.95), 
  tmid = TRUE, iter.disp = FALSE)
}
})
close(pb)
stopCluster(cl)
sim.results.mpl(ic.scen2.r47.n1000.results, beta = true.beta, simdata = n.sim)


#Study 2, n = 1000, left/interval = 75%, right cens = 20%
true.beta = list(c(1,0.5), c(0.5,0.5))
ic.scen2.r20.n1000 <- gen.phcsh.data(n.datasets = 1000, n.obs = 1000,
prob_event = 0.05, rho = 3, lambda1 = 1, lambda2 = 0.5,
betas = true.beta, gammaL = 0.5, gammaR = 1.64)
table(unlist(lapply(ic.scen2.r20.n1000, function(a) a$event))) /
  length(unlist(lapply(ic.scen2.r20.n1000, function(a) a$event)))
#cox midpoint t with n = 1000
ic.scen2.r20.n1000.r1 = ic.scen2.r20.n1000.r2 = list()
for(sim in 1:n.sim){
  ic.scen2.r20.n1000.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, 
  data = ic.scen2.r20.n1000[[sim]])
  ic.scen2.r20.n1000.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, #up to here
  data = ic.scen2.r20.n1000[[sim]])
}
sim.results.cox(ic.scen2.r20.n1000.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(ic.scen2.r20.n1000.r2, beta = true.beta[[2]], simdata = n.sim)

cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
ic.scen2.r20.n1000.time <- system.time({ic.scen2.r20.n1000.results <- foreach(i = 1:n.sim, 
  .combine = rbind, .options.snow = opts) %dopar% {
  phcsh_mpl(Surv(time, time2, event=event, "interval") ~ x1 + x2,
  risk = ic.scen2.r20.n1000[[i]]$risk, data = ic.scen2.r20.n1000[[i]], 
  max.outer = 10, max.iter = 3000, lambda = c(1e-5,1e-5), n.basis = NULL, 
  iknots.pos = NULL, aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, 
  iter.disp = FALSE)
}
})
close(pb)
stopCluster(cl)
sim.results.mpl(ic.scen2.r20.n1000.results, beta = true.beta, simdata = n.sim)


# Study 1 baseline hazard
# int/left cens perc = 47.5%, right cens = 47.5%, n = 200
# calculate values for plots
scen1.r47.n200.min <- max(sapply(ic.scen1.r47.n200, function(a)
  min(na.omit(c(a$time, a$time2)))))
scen1.r47.n200.max <- min(sapply(ic.scen1.r47.n200, function(a)
  max(na.omit(c(a$time, a$time2)))))
scen1.plot.t.r47 <- seq(scen1.r47.n200.min, scen1.r47.n200.max, length.out=1000)
scen1.plot.bh.r47 <- scen1.plot.bh.r47.r2 <- list()
for(sim in 1:1000){
  scen1.plot.bh.r47[[sim]] <- pred.bh(ic.scen1.r47.n200.results[sim,],
  scen1.plot.t.r47, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  scen1.plot.bh.r47.r2[[sim]] <- pred.bh(ic.scen1.r47.n200.results[sim,],
   scen1.plot.t.r47, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}
scen1.plot.bh.r47.all <- lapply(scen1.plot.bh.r47, function(a) a$pred.h0r)
scen1.plot.bh.r47.h0r.mean <- colMeans(Reduce("rbind", scen1.plot.bh.r47.all))
scen1.plot.bh.r47.all.lower <- lapply(scen1.plot.bh.r47, 
function(a) a$pred.h0r.lower)
scen1.plot.bh.r47.h0r.mean.lower<- colMeans(Reduce("rbind",
scen1.plot.bh.r47.all.lower))
scen1.plot.bh.r47.all.upper <- lapply(scen1.plot.bh.r47, 
function(a) a$pred.h0r.upper)
scen1.plot.bh.r47.h0r.mean.upper<- colMeans(Reduce("rbind", 
scen1.plot.bh.r47.all.upper))
bh.cp <- colSums(Reduce("rbind", lapply(scen1.plot.bh.r47, 
function(a) ifelse((a$pred.h0r.lower < a$true.h0r)
& (a$true.h0r < a$pred.h0r.upper),1,0)))) / 1000


scen1.plot.bh.r47.r2.all <- lapply(scen1.plot.bh.r47.r2, function(a) a$pred.h0r)
scen1.plot.bh.r47.r2.h0r.mean <- colMeans(Reduce("rbind", 
            scen1.plot.bh.r47.r2.all))
scen1.plot.bh.r47.r2.all.lower <- lapply(scen1.plot.bh.r47.r2, 
function(a) a$pred.h0r.lower)
scen1.plot.bh.r47.r2.h0r.mean.lower<- colMeans(Reduce("rbind", 
scen1.plot.bh.r47.r2.all.lower))
scen1.plot.bh.r47.r2.all.upper <- lapply(scen1.plot.bh.r47.r2, 
function(a) a$pred.h0r.upper)
scen1.plot.bh.r47.r2.h0r.mean.upper<- colMeans(Reduce("rbind", 
scen1.plot.bh.r47.r2.all.upper))
bh.cp.r2 <- lapply(scen1.plot.bh.r47.r2, function(a) 
ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0))
bh.cp.r2 <- colSums(Reduce("rbind", lapply(scen1.plot.bh.r47.r2, 
function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & 
(a$true.h0r < a$pred.h0r.upper),1,0)))) / 1000

par(mfrow=c(2,2))
plot(scen1.plot.t.r47, scen1.plot.bh.r47.h0r.mean, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(A) Risk 1")
lines(scen1.plot.t.r47, scen1.plot.bh.r47.h0r.mean.lower, col='#003366',
      lty='dashed')
lines(scen1.plot.t.r47, scen1.plot.bh.r47.h0r.mean.upper, col='#003366',
      lty='dashed')
lines(scen1.plot.t.r47, scen1.plot.bh.r47[[1]]$true.h0r, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

plot(scen1.plot.t.r47, scen1.plot.bh.r47.r2.h0r.mean, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(B) Risk 2")
lines(scen1.plot.t.r47, scen1.plot.bh.r47.r2.h0r.mean.lower, col='#003366',
      lty='dashed')
lines(scen1.plot.t.r47, scen1.plot.bh.r47.r2.h0r.mean.upper, col='#003366',
      lty='dashed')
lines(scen1.plot.t.r47, scen1.plot.bh.r47.r2[[1]]$true.h0r, 
      col = "black",cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')


plot(scen1.plot.t.r47, bh.cp.r2, ylab="Coverage probability", xlab="t", 
     col='#003366', type='l', ylim=c(0,1), main="(C) Risk 1")
lines(scen1.plot.t.r47, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), 
       lty=c(1,2), bty = 'n')
plot(scen1.plot.t.r47, bh.cp.r2, ylab="Coverage probability", xlab="t", 
     col='#003366', type='l', ylim=c(0,1), main="(D) Risk 2")
lines(scen1.plot.t.r47, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), 
       lty=c(1,2), bty = 'n')


# CIF for study 1,
# int/left cens perc = 47.5%, right cens = 47.5%, n = 200
scen1.plot.cif.r47 <- scen1.plot.r2.r47.r2 <- list()
for(sim in 1:1000){
  scen1.plot.cif.r47[[sim]] <- pred.CIF(ic.scen1.r47.n200.results[sim,],
  scen1.plot.t.r47, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  
  scen1.plot.r2.r47.r2[[sim]] <- pred.CIF(ic.scen1.r47.n200.results[sim,],
  scen1.plot.t.r47, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}
scen1.plot.cif.r47.all <- lapply(scen1.plot.cif.r47, 
function(a) a$pred.F0r)
scen1.plot.cif.r47.F0r.mean <- colMeans(Reduce("rbind",
scen1.plot.cif.r47.all))
scen1.plot.cif.r47.all.lower <- lapply(scen1.plot.cif.r47, 
function(a) a$pred.F0r.lower)
scen1.plot.cif.r47.F0r.mean.lower<- colMeans(Reduce("rbind", 
scen1.plot.cif.r47.all.lower))
scen1.plot.cif.r47.all.upper <- lapply(scen1.plot.cif.r47, 
function(a) a$pred.F0r.upper)
scen1.plot.cif.r47.F0r.mean.upper<- colMeans(Reduce("rbind", 
scen1.plot.cif.r47.all.upper))

scen1.plot.cif.r47.r2.all <- lapply(scen1.plot.r2.r47.r2, 
function(a) a$pred.F0r)
scen1.plot.cif.r47.r2.F0r.mean <- colMeans(Reduce("rbind", 
scen1.plot.cif.r47.r2.all))
scen1.plot.cif.r47.r2.all.lower <- lapply(scen1.plot.r2.r47.r2, 
function(a) a$pred.F0r.lower)
scen1.plot.cif.r47.r2.F0r.mean.lower<- colMeans(Reduce("rbind", 
scen1.plot.cif.r47.r2.all.lower))
scen1.plot.cif.r47.r2.all.upper <- lapply(scen1.plot.r2.r47.r2, 
function(a) a$pred.F0r.upper)
scen1.plot.cif.r47.r2.F0r.mean.upper<- colMeans(Reduce("rbind", 
scen1.plot.cif.r47.r2.all.upper))

integrand.F0r <- function(x, lambda, lambda1, lambda2, rho){
  (lambda*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho))
}
true.F0r.r1 <- sapply(scen1.plot.t.r47, function(a) integrate(integrand.F0r, 0, 
               a, rho=3, lambda=1, lambda1=1, lambda2=0.5)$value)
true.F0r.r2 <- sapply(scen1.plot.t.r47, function(a) integrate(integrand.F0r, 0, 
               a, rho=3, lambda=0.5, lambda1=1, lambda2=0.5)$value)
F0r.plot.mat <- sapply(scen1.plot.cif.r47, function(a) a$pred.F0r)
F0r.plot.mean <- apply(F0r.plot.mat, 1, mean)
F0r.plot.r2.mat <- sapply(scen1.plot.r2.r47.r2, function(a) a$pred.F0r)
F0r.plot.r2.mean <- apply(F0r.plot.r2.mat, 1, mean)

par(mfrow=c(2,1))
plot(scen1.plot.t.r47, scen1.plot.cif.r47.F0r.mean, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), 
     main="(A) Risk 1", lty=1)
lines(scen1.plot.t.r47, scen1.plot.cif.r47.F0r.mean.lower, col='#003366',
      lty=2)
lines(scen1.plot.t.r47, scen1.plot.cif.r47.F0r.mean.upper, col='#003366',
      lty=2)
lines(scen1.plot.t.r47, true.F0r.r1, col = "black", cex=10, lty=1)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=c(1,1), bty = 'n')

plot(scen1.plot.t.r47, scen1.plot.cif.r47.r2.F0r.mean, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), 
     main="(B) Risk 2", lty=1)
lines(scen1.plot.t.r47, scen1.plot.cif.r47.r2.F0r.mean.lower, col='#003366',
      lty=2)
lines(scen1.plot.t.r47, scen1.plot.cif.r47.r2.F0r.mean.upper, col='#003366',
      lty=2)
lines(scen1.plot.t.r47, true.F0r.r2, col = "black", cex=10, lty=1)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=c(1,1), bty = 'n')


#study 1, n = 1000
# int/left cens perc = 47.5%, right cens = 47.5%, n = 1000
#caluclate values for plots
scen1.r47.n1000.min <- max(sapply(ic.scen1.r47.n1000, function(a)
min(na.omit(c(a$time, a$time2)))))
scen1.r47.n1000.max <- min(sapply(ic.scen1.r47.n1000, function(a)
max(na.omit(c(a$time, a$time2)))))
scen1.plot.t.r47.n1000 <- seq(scen1.r47.n1000.min, scen1.r47.n1000.max, length.out=1000)
scen1.plot.bh.r47.n1000 <- scen1.plot.bh.r47.r2.n1000 <- list()
for(sim in 1:1000){
  scen1.plot.bh.r47.n1000[[sim]] <- pred.bh(ic.scen1.r47.n1000.results[sim,],
  scen1.plot.t.r47.n1000, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  scen1.plot.bh.r47.r2.n1000[[sim]] <- pred.bh(ic.scen1.r47.n1000.results[sim,],
  scen1.plot.t.r47.n1000, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}
scen1.plot.bh.r47.all.n1000 <- lapply(scen1.plot.bh.r47.n1000, function(a) a$pred.h0r)
scen1.plot.bh.r47.h0r.mean.n1000 <- colMeans(Reduce("rbind", scen1.plot.bh.r47.all.n1000))
scen1.plot.bh.r47.all.lower.n1000 <- lapply(scen1.plot.bh.r47.n1000, 
                                      function(a) a$pred.h0r.lower)
scen1.plot.bh.r47.h0r.mean.lower.n1000 <- colMeans(Reduce("rbind",
                                          scen1.plot.bh.r47.all.lower.n1000))
scen1.plot.bh.r47.all.upper.n1000 <- lapply(scen1.plot.bh.r47.n1000, 
                                      function(a) a$pred.h0r.upper)
scen1.plot.bh.r47.h0r.mean.upper.n1000 <- colMeans(Reduce("rbind", 
                                        scen1.plot.bh.r47.all.upper.n1000))
bh.cp.n1000 <- colSums(Reduce("rbind", lapply(scen1.plot.bh.r47.n1000, 
         function(a) ifelse((a$pred.h0r.lower < a$true.h0r)
         & (a$true.h0r < a$pred.h0r.upper),1,0)))) / 1000


scen1.plot.bh.r47.r2.all.n1000 <- lapply(scen1.plot.bh.r47.r2.n1000, 
                                  function(a) a$pred.h0r)
scen1.plot.bh.r47.r2.h0r.mean.n1000 <- colMeans(Reduce("rbind", 
                                      scen1.plot.bh.r47.r2.all.n1000))
scen1.plot.bh.r47.r2.all.lower.n1000 <- lapply(scen1.plot.bh.r47.r2.n1000, 
                                         function(a) a$pred.h0r.lower)
scen1.plot.bh.r47.r2.h0r.mean.lower.n1000 <- colMeans(Reduce("rbind", 
                                         scen1.plot.bh.r47.r2.all.lower.n1000))
scen1.plot.bh.r47.r2.all.upper.n1000 <- lapply(scen1.plot.bh.r47.r2.n1000, 
                                         function(a) a$pred.h0r.upper)
scen1.plot.bh.r47.r2.h0r.mean.upper.n1000 <- colMeans(Reduce("rbind", 
                                    scen1.plot.bh.r47.r2.all.upper.n1000))
bh.cp.r2.n1000 <- lapply(scen1.plot.bh.r47.r2.n1000, function(a) 
  ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0))
bh.cp.r2.n1000 <- colSums(Reduce("rbind", lapply(scen1.plot.bh.r47.r2.n1000, 
 function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & 
(a$true.h0r < a$pred.h0r.upper),1,0)))) / 1000

par(mfrow=c(2,2))
plot(scen1.plot.t.r47.n1000, scen1.plot.bh.r47.h0r.mean.n1000, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(A) Risk 1")
lines(scen1.plot.t.r47.n1000, scen1.plot.bh.r47.h0r.mean.lower.n1000, col='#003366',
      lty='dashed')
lines(scen1.plot.t.r47.n1000, scen1.plot.bh.r47.h0r.mean.upper.n1000, col='#003366',
      lty='dashed')
lines(scen1.plot.t.r47.n1000, scen1.plot.bh.r47.n1000[[1]]$true.h0r, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

plot(scen1.plot.t.r47.n1000, scen1.plot.bh.r47.r2.h0r.mean.n1000, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(B) Risk 2")
lines(scen1.plot.t.r47.n1000, scen1.plot.bh.r47.r2.h0r.mean.lower.n1000, col='#003366',
      lty='dashed')
lines(scen1.plot.t.r47.n1000, scen1.plot.bh.r47.r2.h0r.mean.upper.n1000, col='#003366',
      lty='dashed')
lines(scen1.plot.t.r47.n1000, scen1.plot.bh.r47.r2.n1000[[1]]$true.h0r, 
      col = "black",cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')


plot(scen1.plot.t.r47.n1000, bh.cp.r2.n1000, ylab="Coverage probability", xlab="t", 
     col='#003366', type='l', ylim=c(0,1), main="(C) Risk 1")
lines(scen1.plot.t.r47.n1000, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), 
       lty=c(1,2), bty = 'n')
plot(scen1.plot.t.r47.n1000, bh.cp.r2.n1000, ylab="Coverage probability", xlab="t", 
     col='#003366', type='l', ylim=c(0,1), main="(D) Risk 2")
lines(scen1.plot.t.r47.n1000, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), 
       lty=c(1,2), bty = 'n')

#figure 2, CIF for scenario 1,
# int/left cens perc = 47.5%, right cens = 47.5%, n = 1000
scen1.plot.cif.r47.n1000 <- scen1.plot.r2.r47.r2.n1000 <- list()
for(sim in 1:1000){
  scen1.plot.cif.r47.n1000[[sim]] <- pred.CIF(ic.scen1.r47.n1000.results[sim,],
  scen1.plot.t.r47.n1000, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  scen1.plot.r2.r47.r2.n1000[[sim]] <- pred.CIF(ic.scen1.r47.n1000.results[sim,],
  scen1.plot.t.r47.n1000, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}
scen1.plot.cif.r47.all.n1000 <- lapply(scen1.plot.cif.r47.n1000, 
                                 function(a) a$pred.F0r)
scen1.plot.cif.r47.F0r.mean.n1000 <- colMeans(Reduce("rbind",
                                               scen1.plot.cif.r47.all.n1000))
scen1.plot.cif.r47.all.lower.n1000 <- lapply(scen1.plot.cif.r47.n1000, 
                                       function(a) a$pred.F0r.lower)
scen1.plot.cif.r47.F0r.mean.lower.n1000 <- colMeans(Reduce("rbind", 
                                    scen1.plot.cif.r47.all.lower.n1000))
scen1.plot.cif.r47.all.upper.n1000 <- lapply(scen1.plot.cif.r47.n1000, 
                                       function(a) a$pred.F0r.upper)
scen1.plot.cif.r47.F0r.mean.upper.n1000 <- colMeans(Reduce("rbind", 
                                      scen1.plot.cif.r47.all.upper.n1000))

scen1.plot.cif.r47.r2.all.n1000 <- lapply(scen1.plot.r2.r47.r2.n1000, 
                                    function(a) a$pred.F0r)
scen1.plot.cif.r47.r2.F0r.mean.n1000 <- colMeans(Reduce("rbind", 
                                   scen1.plot.cif.r47.r2.all.n1000))
scen1.plot.cif.r47.r2.all.lower.n1000 <- lapply(scen1.plot.r2.r47.r2.n1000, 
                                          function(a) a$pred.F0r.lower)
scen1.plot.cif.r47.r2.F0r.mean.lower.n1000 <- colMeans(Reduce("rbind", 
                                  scen1.plot.cif.r47.r2.all.lower.n1000))
scen1.plot.cif.r47.r2.all.upper.n1000 <- lapply(scen1.plot.r2.r47.r2.n1000, 
                                          function(a) a$pred.F0r.upper)
scen1.plot.cif.r47.r2.F0r.mean.upper.n1000 <- colMeans(Reduce("rbind", 
                          scen1.plot.cif.r47.r2.all.upper.n1000))

integrand.F0r <- function(x, lambda, lambda1, lambda2, rho){
  (lambda*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho))
}
true.F0r.r1.n1000 <- sapply(scen1.plot.t.r47.n1000, function(a) integrate(integrand.F0r, 0, 
        a, rho=3, lambda=1, lambda1=1, lambda2=0.5)$value)
true.F0r.r2.n1000 <- sapply(scen1.plot.t.r47.n1000, function(a) integrate(integrand.F0r, 0, 
        a, rho=3, lambda=0.5, lambda1=1, lambda2=0.5)$value)
F0r.plot.mat.n1000 <- sapply(scen1.plot.cif.r47.n1000, function(a) a$pred.F0r)
F0r.plot.mean.n1000 <- apply(F0r.plot.mat.n1000, 1, mean)
F0r.plot.r2.mat.n1000 <- sapply(scen1.plot.r2.r47.r2.n1000, 
                      function(a) a$pred.F0r)
F0r.plot.r2.mean.n1000 <- apply(F0r.plot.r2.mat.n1000, 1, mean)

par(mfrow=c(2,1))
plot(scen1.plot.t.r47.n1000, scen1.plot.cif.r47.F0r.mean.n1000, type='l', 
    col='#003366',xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), 
     main="(A) Risk 1", lty=1)
lines(scen1.plot.t.r47.n1000, scen1.plot.cif.r47.F0r.mean.lower.n1000, 
      col='#003366', lty=2)
lines(scen1.plot.t.r47.n1000, scen1.plot.cif.r47.F0r.mean.upper.n1000, col='#003366',
      lty=2)
lines(scen1.plot.t.r47.n1000, true.F0r.r1.n1000, col = "black", cex=10, lty=1)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=c(1,1), bty = 'n')

plot(scen1.plot.t.r47.n1000, scen1.plot.cif.r47.r2.F0r.mean.n1000, type='l', 
     col='#003366', xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), 
     main="(B) Risk 2", lty=1)
lines(scen1.plot.t.r47.n1000, scen1.plot.cif.r47.r2.F0r.mean.lower.n1000, 
      col='#003366', lty=2)
lines(scen1.plot.t.r47.n1000, scen1.plot.cif.r47.r2.F0r.mean.upper.n1000, 
      col='#003366', lty=2)
lines(scen1.plot.t.r47.n1000, true.F0r.r2.n1000, col = "black", cex=10, lty=1)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=c(1,1), bty = 'n')






#study 1 baseline hazard,
# int/left cens perc = 75%, right cens = 20%, n = 200
scen1.r20.n200.min <- max(sapply(ic.scen1.r20.n200, function(a)
  min(na.omit(c(a$time, a$time2)))))
scen1.r20.n200.max <- min(sapply(ic.scen1.r20.n200, function(a)
  max(na.omit(c(a$time, a$time2)))))
scen1.plot.t.r20 <- seq(scen1.r20.n200.min, scen1.r20.n200.max, length.out=1000)
scen1.plot.bh.r20 <- scen1.plot.bh.r20.r2 <- list()
for(sim in 1:1000){
  scen1.plot.bh.r20[[sim]] <- pred.bh(ic.scen1.r20.n200.results[sim,],
    scen1.plot.t.r20, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  scen1.plot.bh.r20.r2[[sim]] <- pred.bh(ic.scen1.r20.n200.results[sim,],
   scen1.plot.t.r47, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}
scen1.plot.bh.r20.all <- lapply(scen1.plot.bh.r20, function(a) a$pred.h0r)
scen1.plot.bh.r20.h0r.mean <- colMeans(Reduce("rbind", scen1.plot.bh.r20.all))
scen1.plot.bh.r20.all.lower <- lapply(scen1.plot.bh.r20, 
                                      function(a) a$pred.h0r.lower)
scen1.plot.bh.r20.h0r.mean.lower<- colMeans(Reduce("rbind",
                                                   scen1.plot.bh.r20.all.lower))
scen1.plot.bh.r20.all.upper <- lapply(scen1.plot.bh.r20, 
                                      function(a) a$pred.h0r.upper)
scen1.plot.bh.r20.h0r.mean.upper<- colMeans(Reduce("rbind", 
                                                   scen1.plot.bh.r20.all.upper))
bh.cp <- colSums(Reduce("rbind", lapply(scen1.plot.bh.r20, 
    function(a) ifelse((a$pred.h0r.lower < a$true.h0r)
      & (a$true.h0r < a$pred.h0r.upper),1,0)))) / 1000


scen1.plot.bh.r20.r2.all <- lapply(scen1.plot.bh.r20.r2, function(a) a$pred.h0r)
scen1.plot.bh.r20.r2.h0r.mean <- colMeans(Reduce("rbind", 
    scen1.plot.bh.r20.r2.all))
scen1.plot.bh.r20.r2.all.lower <- lapply(scen1.plot.bh.r20.r2, 
                                         function(a) a$pred.h0r.lower)
scen1.plot.bh.r20.r2.h0r.mean.lower<- colMeans(Reduce("rbind", 
     scen1.plot.bh.r20.r2.all.lower))
scen1.plot.bh.r20.r2.all.upper <- lapply(scen1.plot.bh.r20.r2, 
                                         function(a) a$pred.h0r.upper)
scen1.plot.bh.r20.r2.h0r.mean.upper<- colMeans(Reduce("rbind", 
     scen1.plot.bh.r20.r2.all.upper))
bh.cp.r2 <- lapply(scen1.plot.bh.r20.r2, function(a) 
  ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0))
bh.cp.r2 <- colSums(Reduce("rbind", lapply(scen1.plot.bh.r20.r2, 
  function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & 
  (a$true.h0r < a$pred.h0r.upper),1,0)))) / 1000

par(mfrow=c(2,2))
plot(scen1.plot.t.r20, scen1.plot.bh.r20.h0r.mean, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(A) Risk 1")
lines(scen1.plot.t.r20, scen1.plot.bh.r20.h0r.mean.lower, col='#003366',
      lty='dashed')
lines(scen1.plot.t.r20, scen1.plot.bh.r20.h0r.mean.upper, col='#003366',
      lty='dashed')
lines(scen1.plot.t.r20, scen1.plot.bh.r20[[1]]$true.h0r, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

plot(scen1.plot.t.r20, scen1.plot.bh.r20.r2.h0r.mean, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(B) Risk 2")
lines(scen1.plot.t.r20, scen1.plot.bh.r20.r2.h0r.mean.lower, col='#003366',
      lty='dashed')
lines(scen1.plot.t.r20, scen1.plot.bh.r20.r2.h0r.mean.upper, col='#003366',
      lty='dashed')
lines(scen1.plot.t.r20, scen1.plot.bh.r20.r2[[1]]$true.h0r, 
      col = "black",cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')


plot(scen1.plot.t.r20, bh.cp.r2, ylab="Coverage probability", xlab="t", 
     col='#003366', type='l', ylim=c(0,1), main="(C) Risk 1")
lines(scen1.plot.t.r20, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), 
       lty=c(1,2), bty = 'n')
plot(scen1.plot.t.r20, bh.cp.r2, ylab="Coverage probability", xlab="t", 
     col='#003366', type='l', ylim=c(0,1), main="(D) Risk 2")
lines(scen1.plot.t.r20, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), 
       lty=c(1,2), bty = 'n')


# CIF for study 1,
# int/left cens perc = 75%, right cens = 20%, n = 200
scen1.plot.cif.r20 <- scen1.plot.r2.r20.r2 <- list()
for(sim in 1:1000){
  scen1.plot.cif.r20[[sim]] <- pred.CIF(ic.scen1.r20.n200.results[sim,],
  scen1.plot.t.r20, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  
  scen1.plot.r2.r20.r2[[sim]] <- pred.CIF(ic.scen1.r20.n200.results[sim,],
  scen1.plot.t.r20, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}
scen1.plot.cif.r20.all <- lapply(scen1.plot.cif.r20, 
function(a) a$pred.F0r)
scen1.plot.cif.r20.F0r.mean <- colMeans(Reduce("rbind",
scen1.plot.cif.r20.all))
scen1.plot.cif.r20.all.lower <- lapply(scen1.plot.cif.r20, 
function(a) a$pred.F0r.lower)
scen1.plot.cif.r20.F0r.mean.lower<- colMeans(Reduce("rbind", 
scen1.plot.cif.r20.all.lower))
scen1.plot.cif.r20.all.upper <- lapply(scen1.plot.cif.r20, 
function(a) a$pred.F0r.upper)
scen1.plot.cif.r20.F0r.mean.upper<- colMeans(Reduce("rbind", 
scen1.plot.cif.r20.all.upper))

scen1.plot.cif.r20.r2.all <- lapply(scen1.plot.r2.r20.r2, 
function(a) a$pred.F0r)
scen1.plot.cif.r20.r2.F0r.mean <- colMeans(Reduce("rbind", 
scen1.plot.cif.r20.r2.all))
scen1.plot.cif.r20.r2.all.lower <- lapply(scen1.plot.r2.r20.r2, 
function(a) a$pred.F0r.lower)
scen1.plot.cif.r20.r2.F0r.mean.lower<- colMeans(Reduce("rbind", 
scen1.plot.cif.r20.r2.all.lower))
scen1.plot.cif.r20.r2.all.upper <- lapply(scen1.plot.r2.r20.r2, 
function(a) a$pred.F0r.upper)
scen1.plot.cif.r20.r2.F0r.mean.upper<- colMeans(Reduce("rbind", 
scen1.plot.cif.r20.r2.all.upper))

integrand.F0r <- function(x, lambda, lambda1, lambda2, rho){
  (lambda*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho))
}
true.F0r.r1 <- sapply(scen1.plot.t.r20, function(a) integrate(integrand.F0r, 0, 
               a, rho=3, lambda=1, lambda1=1, lambda2=0.5)$value)
true.F0r.r2 <- sapply(scen1.plot.t.r20, function(a) integrate(integrand.F0r, 0, 
               a, rho=3, lambda=0.5, lambda1=1, lambda2=0.5)$value)
F0r.plot.mat <- sapply(scen1.plot.cif.r20, function(a) a$pred.F0r)
F0r.plot.mean <- apply(F0r.plot.mat, 1, mean)
F0r.plot.r2.mat <- sapply(scen1.plot.r2.r20.r2, function(a) a$pred.F0r)
F0r.plot.r2.mean <- apply(F0r.plot.r2.mat, 1, mean)

par(mfrow=c(2,1))
plot(scen1.plot.t.r20, scen1.plot.cif.r20.F0r.mean, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), 
     main="(A) Risk 1", lty=1)
lines(scen1.plot.t.r20, scen1.plot.cif.r20.F0r.mean.lower, col='#003366',
      lty=2)
lines(scen1.plot.t.r20, scen1.plot.cif.r20.F0r.mean.upper, col='#003366',
      lty=2)
lines(scen1.plot.t.r20, true.F0r.r1, col = "black", cex=10, lty=1)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=c(1,1), bty = 'n')

plot(scen1.plot.t.r20, scen1.plot.cif.r20.r2.F0r.mean, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), 
     main="(B) Risk 2", lty=1)
lines(scen1.plot.t.r20, scen1.plot.cif.r20.r2.F0r.mean.lower, col='#003366',
      lty=2)
lines(scen1.plot.t.r20, scen1.plot.cif.r20.r2.F0r.mean.upper, col='#003366',
      lty=2)
lines(scen1.plot.t.r20, true.F0r.r2, col = "black", cex=10, lty=1)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=c(1,1), bty = 'n')


#study 1 baseline hazard,
# int/left cens perc = 75%, right cens = 20%, n = 1000
scen1.r20.n1000.min <- max(sapply(ic.scen1.r20.n1000, function(a)
  min(na.omit(c(a$time, a$time2)))))
scen1.r20.n1000.max <- min(sapply(ic.scen1.r20.n1000, function(a)
  max(na.omit(c(a$time, a$time2)))))
scen1.plot.t.r20.n1000 <- seq(scen1.r20.n1000.min, scen1.r20.n1000.max, 
                              length.out=1000)
scen1.plot.bh.r20.n1000 <- scen1.plot.bh.r20.r2.n1000 <- list()
for(sim in 1:1000){
  scen1.plot.bh.r20.n1000[[sim]] <- pred.bh(ic.scen1.r20.n1000.results[sim,],
      scen1.plot.t.r20.n1000, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  scen1.plot.bh.r20.r2.n1000[[sim]] <- pred.bh(ic.scen1.r20.n1000.results[sim,],
      scen1.plot.t.r47.n1000, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}
scen1.plot.bh.r20.all.n1000 <- lapply(scen1.plot.bh.r20.n1000, 
                                function(a) a$pred.h0r)
scen1.plot.bh.r20.h0r.mean.n1000 <- colMeans(Reduce("rbind", 
                              scen1.plot.bh.r20.all.n1000))
scen1.plot.bh.r20.all.lower.n1000 <- lapply(scen1.plot.bh.r20.n1000, 
                                      function(a) a$pred.h0r.lower)
scen1.plot.bh.r20.h0r.mean.lower.n1000 <- colMeans(Reduce("rbind",
                                scen1.plot.bh.r20.all.lower.n1000))
scen1.plot.bh.r20.all.upper.n1000 <- lapply(scen1.plot.bh.r20.n1000, 
                                      function(a) a$pred.h0r.upper)
scen1.plot.bh.r20.h0r.mean.upper.n1000 <- colMeans(Reduce("rbind", 
                              scen1.plot.bh.r20.all.upper.n1000))
bh.cp.n1000 <- colSums(Reduce("rbind", lapply(scen1.plot.bh.r20.n1000, 
         function(a) ifelse((a$pred.h0r.lower < a$true.h0r)
         & (a$true.h0r < a$pred.h0r.upper),1,0)))) / 1000


scen1.plot.bh.r20.r2.all.n1000 <- lapply(scen1.plot.bh.r20.r2.n1000, 
                function(a) a$pred.h0r)
scen1.plot.bh.r20.r2.h0r.mean.n1000 <- colMeans(Reduce("rbind", 
                 scen1.plot.bh.r20.r2.all.n1000))
scen1.plot.bh.r20.r2.all.lower.n1000 <- lapply(scen1.plot.bh.r20.r2.n1000, 
                                         function(a) a$pred.h0r.lower)
scen1.plot.bh.r20.r2.h0r.mean.lower.n1000 <- colMeans(Reduce("rbind", 
                   scen1.plot.bh.r20.r2.all.lower.n1000))
scen1.plot.bh.r20.r2.all.upper.n1000 <- lapply(scen1.plot.bh.r20.r2.n1000, 
              function(a) a$pred.h0r.upper)
scen1.plot.bh.r20.r2.h0r.mean.upper.n1000 <- colMeans(Reduce("rbind", 
                     scen1.plot.bh.r20.r2.all.upper.n1000))
bh.cp.r2.n1000 <- lapply(scen1.plot.bh.r20.r2.n1000, function(a) 
  ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0))
bh.cp.r2.n1000 <- colSums(Reduce("rbind", lapply(scen1.plot.bh.r20.r2.n1000, 
    function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & 
      (a$true.h0r < a$pred.h0r.upper),1,0)))) / 1000

par(mfrow=c(2,2))
plot(scen1.plot.t.r20.n1000, scen1.plot.bh.r20.h0r.mean.n1000, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(A) Risk 1")
lines(scen1.plot.t.r20.n1000, scen1.plot.bh.r20.h0r.mean.lower.n1000, col='#003366',
      lty='dashed')
lines(scen1.plot.t.r20.n1000, scen1.plot.bh.r20.h0r.mean.upper.n1000, col='#003366',
      lty='dashed')
lines(scen1.plot.t.r20.n1000, scen1.plot.bh.r20.n1000[[1]]$true.h0r, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

plot(scen1.plot.t.r20.n1000, scen1.plot.bh.r20.r2.h0r.mean.n1000, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(B) Risk 2")
lines(scen1.plot.t.r20.n1000, scen1.plot.bh.r20.r2.h0r.mean.lower.n1000, col='#003366',
      lty='dashed')
lines(scen1.plot.t.r20.n1000, scen1.plot.bh.r20.r2.h0r.mean.upper.n1000, col='#003366',
      lty='dashed')
lines(scen1.plot.t.r20.n1000, scen1.plot.bh.r20.r2.n1000[[1]]$true.h0r, 
      col = "black",cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')


plot(scen1.plot.t.r20.n1000, bh.cp.r2.n1000, ylab="Coverage probability", xlab="t", 
     col='#003366', type='l', ylim=c(0,1), main="(C) Risk 1")
lines(scen1.plot.t.r20.n1000, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), 
       lty=c(1,2), bty = 'n')
plot(scen1.plot.t.r20.n1000, bh.cp.r2.n1000, ylab="Coverage probability", xlab="t", 
     col='#003366', type='l', ylim=c(0,1), main="(D) Risk 2")
lines(scen1.plot.t.r20.n1000, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), 
       lty=c(1,2), bty = 'n')


# CIF for study 1,
# int/left cens perc = 75%, right cens = 20%, n = 1000
scen1.plot.cif.r20.n1000 <- scen1.plot.r2.r20.r2.n1000 <- list()
for(sim in 1:1000){
  scen1.plot.cif.r20.n1000[[sim]] <- pred.CIF(ic.scen1.r20.n1000.results[sim,],
    scen1.plot.t.r20.n1000, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  
  scen1.plot.r2.r20.r2.n1000[[sim]] <- pred.CIF(ic.scen1.r20.n1000.results[sim,],
    scen1.plot.t.r20.n1000, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}
scen1.plot.cif.r20.all.n1000 <- lapply(scen1.plot.cif.r20.n1000, 
                                 function(a) a$pred.F0r)
scen1.plot.cif.r20.F0r.mean.n1000 <- colMeans(Reduce("rbind",
                                               scen1.plot.cif.r20.all.n1000))
scen1.plot.cif.r20.all.lower.n1000 <- lapply(scen1.plot.cif.r20.n1000, 
                                       function(a) a$pred.F0r.lower)
scen1.plot.cif.r20.F0r.mean.lower.n1000 <- colMeans(Reduce("rbind", 
             scen1.plot.cif.r20.all.lower.n1000))
scen1.plot.cif.r20.all.upper.n1000 <- lapply(scen1.plot.cif.r20.n1000, 
                                       function(a) a$pred.F0r.upper)
scen1.plot.cif.r20.F0r.mean.upper.n1000 <- colMeans(Reduce("rbind", 
                               scen1.plot.cif.r20.all.upper.n1000))

scen1.plot.cif.r20.r2.all.n1000 <- lapply(scen1.plot.r2.r20.r2.n1000, 
                                    function(a) a$pred.F0r)
scen1.plot.cif.r20.r2.F0r.mean.n1000 <- colMeans(Reduce("rbind", 
                                    scen1.plot.cif.r20.r2.all.n1000))
scen1.plot.cif.r20.r2.all.lower.n1000 <- lapply(scen1.plot.r2.r20.r2.n1000, 
                                          function(a) a$pred.F0r.lower)
scen1.plot.cif.r20.r2.F0r.mean.lower.n1000 <- colMeans(Reduce("rbind", 
                           scen1.plot.cif.r20.r2.all.lower.n1000))
scen1.plot.cif.r20.r2.all.upper.n1000 <- lapply(scen1.plot.r2.r20.r2.n1000, 
                                          function(a) a$pred.F0r.upper)
scen1.plot.cif.r20.r2.F0r.mean.upper.n1000 <- colMeans(Reduce("rbind", 
                                     scen1.plot.cif.r20.r2.all.upper.n1000))

integrand.F0r <- function(x, lambda, lambda1, lambda2, rho){
  (lambda*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho))
}
true.F0r.r1.n1000 <- sapply(scen1.plot.t.r20.n1000, function(a) integrate(integrand.F0r, 0, 
                         a, rho=3, lambda=1, lambda1=1, lambda2=0.5)$value)
true.F0r.r2.n1000 <- sapply(scen1.plot.t.r20.n1000, function(a) integrate(integrand.F0r, 0, 
                          a, rho=3, lambda=0.5, lambda1=1, lambda2=0.5)$value)
F0r.plot.mat.n1000 <- sapply(scen1.plot.cif.r20.n1000, function(a) a$pred.F0r)
F0r.plot.mean.n1000 <- apply(F0r.plot.mat.n1000, 1, mean)
F0r.plot.r2.mat.n1000 <- sapply(scen1.plot.r2.r20.r2.n1000, function(a) a$pred.F0r)
F0r.plot.r2.mean.n1000 <- apply(F0r.plot.r2.mat.n1000, 1, mean)

par(mfrow=c(2,1))
plot(scen1.plot.t.r20.n1000, scen1.plot.cif.r20.F0r.mean.n1000, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), 
     main="(A) Risk 1", lty=1)
lines(scen1.plot.t.r20.n1000, scen1.plot.cif.r20.F0r.mean.lower.n1000, col='#003366',
      lty=2)
lines(scen1.plot.t.r20.n1000, scen1.plot.cif.r20.F0r.mean.upper.n1000, col='#003366',
      lty=2)
lines(scen1.plot.t.r20.n1000, true.F0r.r1.n1000, col = "black", cex=10, lty=1)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=c(1,1), bty = 'n')

plot(scen1.plot.t.r20.n1000, scen1.plot.cif.r20.r2.F0r.mean.n1000, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), 
     main="(B) Risk 2", lty=1)
lines(scen1.plot.t.r20.n1000, scen1.plot.cif.r20.r2.F0r.mean.lower.n1000, col='#003366',
      lty=2)
lines(scen1.plot.t.r20.n1000, scen1.plot.cif.r20.r2.F0r.mean.upper.n1000, col='#003366',
      lty=2)
lines(scen1.plot.t.r20.n1000, true.F0r.r2.n1000, col = "black", cex=10, lty=1)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=c(1,1), bty = 'n')


#study 2 baseline hazard graph, right censoring = 47.5%, interval = 47.5%,
# n = 200
scen2.r47.n200.min <- max(sapply(ic.scen2.r47.n200, function(a)
  min(na.omit(c(a$time, a$time2)))))
scen2.r47.n200.max <- min(sapply(ic.scen2.r47.n200, function(a)
  max(na.omit(c(a$time, a$time2)))))
scen2.plot.t.r47 <- seq(scen2.r47.n200.min, scen2.r47.n200.max, length.out=1000)
scen2.plot.bh.r47 <- scen2.plot.bh.r47.r2 <- list()
for(sim in 1:1000){
  scen2.plot.bh.r47[[sim]] <- pred.bh(ic.scen2.r47.n200.results[sim,],
                                      scen2.plot.t.r47, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  scen2.plot.bh.r47.r2[[sim]] <- pred.bh(ic.scen2.r47.n200.results[sim,],
                                         scen2.plot.t.r47, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}
scen2.plot.bh.r47.all <- lapply(scen2.plot.bh.r47, function(a) a$pred.h0r)
scen2.plot.bh.r47.h0r.mean <- colMeans(Reduce("rbind", scen2.plot.bh.r47.all))
scen2.plot.bh.r47.all.lower <- lapply(scen2.plot.bh.r47, 
                                      function(a) a$pred.h0r.lower)
scen2.plot.bh.r47.h0r.mean.lower<- colMeans(Reduce("rbind",
                                                   scen2.plot.bh.r47.all.lower))
scen2.plot.bh.r47.all.upper <- lapply(scen2.plot.bh.r47, 
                                      function(a) a$pred.h0r.upper)
scen2.plot.bh.r47.h0r.mean.upper<- colMeans(Reduce("rbind", 
                                                   scen2.plot.bh.r47.all.upper))
bh.cp <- colSums(Reduce("rbind", lapply(scen2.plot.bh.r47, 
                                        function(a) ifelse((a$pred.h0r.lower < a$true.h0r)
                                                           & (a$true.h0r < a$pred.h0r.upper),1,0)))) / 1000


scen2.plot.bh.r47.r2.all <- lapply(scen2.plot.bh.r47.r2, function(a) a$pred.h0r)
scen2.plot.bh.r47.r2.h0r.mean <- colMeans(Reduce("rbind", 
                                                 scen2.plot.bh.r47.r2.all))
scen2.plot.bh.r47.r2.all.lower <- lapply(scen2.plot.bh.r47.r2, 
                                         function(a) a$pred.h0r.lower)
scen2.plot.bh.r47.r2.h0r.mean.lower<- colMeans(Reduce("rbind", 
                                                      scen2.plot.bh.r47.r2.all.lower))
scen2.plot.bh.r47.r2.all.upper <- lapply(scen2.plot.bh.r47.r2, 
                                         function(a) a$pred.h0r.upper)
scen2.plot.bh.r47.r2.h0r.mean.upper<- colMeans(Reduce("rbind", 
                                                      scen2.plot.bh.r47.r2.all.upper))
bh.cp.r2 <- lapply(scen2.plot.bh.r47.r2, function(a) 
  ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0))
bh.cp.r2 <- colSums(Reduce("rbind", lapply(scen2.plot.bh.r47.r2, 
                                           function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & 
                                                                (a$true.h0r < a$pred.h0r.upper),1,0)))) / 1000

par(mfrow=c(2,2))
plot(scen2.plot.t.r47, scen2.plot.bh.r47.h0r.mean, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(A) Risk 1")
lines(scen2.plot.t.r47, scen2.plot.bh.r47.h0r.mean.lower, col='#003366',
      lty='dashed')
lines(scen2.plot.t.r47, scen2.plot.bh.r47.h0r.mean.upper, col='#003366',
      lty='dashed')
lines(scen2.plot.t.r47, scen2.plot.bh.r47[[1]]$true.h0r, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

plot(scen2.plot.t.r47, scen2.plot.bh.r47.r2.h0r.mean, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(B) Risk 2")
lines(scen2.plot.t.r47, scen2.plot.bh.r47.r2.h0r.mean.lower, col='#003366',
      lty='dashed')
lines(scen2.plot.t.r47, scen2.plot.bh.r47.r2.h0r.mean.upper, col='#003366',
      lty='dashed')
lines(scen2.plot.t.r47, scen2.plot.bh.r47.r2[[1]]$true.h0r, 
      col = "black",cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')


plot(scen2.plot.t.r47, bh.cp.r2, ylab="Coverage probability", xlab="t", 
     col='#003366', type='l', ylim=c(0,1), main="(C) Risk 1")
lines(scen2.plot.t.r47, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), 
       lty=c(1,2), bty = 'n')
plot(scen2.plot.t.r47, bh.cp.r2, ylab="Coverage probability", xlab="t", 
     col='#003366', type='l', ylim=c(0,1), main="(D) Risk 2")
lines(scen2.plot.t.r47, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), 
       lty=c(1,2), bty = 'n')


# CIF for study 2
# int/left cens perc = 47.5%, right cens = 47.5%, n = 200
scen2.plot.cif.r47 <- scen2.plot.r2.r47.r2 <- list()
for(sim in 1:1000){
  scen2.plot.cif.r47[[sim]] <- pred.CIF(ic.scen2.r47.n200.results[sim,],
                                        scen2.plot.t.r47, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  
  scen2.plot.r2.r47.r2[[sim]] <- pred.CIF(ic.scen2.r47.n200.results[sim,],
                                          scen2.plot.t.r47, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}
scen2.plot.cif.r47.all <- lapply(scen2.plot.cif.r47, 
                                 function(a) a$pred.F0r)
scen2.plot.cif.r47.F0r.mean <- colMeans(Reduce("rbind",
                                               scen2.plot.cif.r47.all))
scen2.plot.cif.r47.all.lower <- lapply(scen2.plot.cif.r47, 
                                       function(a) a$pred.F0r.lower)
scen2.plot.cif.r47.F0r.mean.lower<- colMeans(Reduce("rbind", 
                                                    scen2.plot.cif.r47.all.lower))
scen2.plot.cif.r47.all.upper <- lapply(scen2.plot.cif.r47, 
                                       function(a) a$pred.F0r.upper)
scen2.plot.cif.r47.F0r.mean.upper<- colMeans(Reduce("rbind", 
                                                    scen2.plot.cif.r47.all.upper))

scen2.plot.cif.r47.r2.all <- lapply(scen2.plot.r2.r47.r2, 
                                    function(a) a$pred.F0r)
scen2.plot.cif.r47.r2.F0r.mean <- colMeans(Reduce("rbind", 
                                                  scen2.plot.cif.r47.r2.all))
scen2.plot.cif.r47.r2.all.lower <- lapply(scen2.plot.r2.r47.r2, 
                                          function(a) a$pred.F0r.lower)
scen2.plot.cif.r47.r2.F0r.mean.lower<- colMeans(Reduce("rbind", 
                                                       scen2.plot.cif.r47.r2.all.lower))
scen2.plot.cif.r47.r2.all.upper <- lapply(scen2.plot.r2.r47.r2, 
                                          function(a) a$pred.F0r.upper)
scen2.plot.cif.r47.r2.F0r.mean.upper<- colMeans(Reduce("rbind", 
                                                       scen2.plot.cif.r47.r2.all.upper))

integrand.F0r <- function(x, lambda, lambda1, lambda2, rho){
  (lambda*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho))
}
true.F0r.r1 <- sapply(scen2.plot.t.r47, function(a) integrate(integrand.F0r, 0, 
                                                              a, rho=3, lambda=1, lambda1=1, lambda2=0.5)$value)
true.F0r.r2 <- sapply(scen2.plot.t.r47, function(a) integrate(integrand.F0r, 0, 
                                                              a, rho=3, lambda=0.5, lambda1=1, lambda2=0.5)$value)
F0r.plot.mat <- sapply(scen2.plot.cif.r47, function(a) a$pred.F0r)
F0r.plot.mean <- apply(F0r.plot.mat, 1, mean)
F0r.plot.r2.mat <- sapply(scen2.plot.r2.r47.r2, function(a) a$pred.F0r)
F0r.plot.r2.mean <- apply(F0r.plot.r2.mat, 1, mean)

par(mfrow=c(2,1))
plot(scen2.plot.t.r47, scen2.plot.cif.r47.F0r.mean, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), 
     main="(A) Risk 1", lty=1)
lines(scen2.plot.t.r47, scen2.plot.cif.r47.F0r.mean.lower, col='#003366',
      lty=2)
lines(scen2.plot.t.r47, scen2.plot.cif.r47.F0r.mean.upper, col='#003366',
      lty=2)
lines(scen2.plot.t.r47, true.F0r.r1, col = "black", cex=10, lty=1)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=c(1,1), bty = 'n')

plot(scen2.plot.t.r47, scen1.plot.cif.r47.r2.F0r.mean, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), 
     main="(B) Risk 2", lty=1)
lines(scen2.plot.t.r47, scen2.plot.cif.r47.r2.F0r.mean.lower, col='#003366',
      lty=2)
lines(scen2.plot.t.r47, scen2.plot.cif.r47.r2.F0r.mean.upper, col='#003366',
      lty=2)
lines(scen2.plot.t.r47, true.F0r.r2, col = "black", cex=10, lty=1)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=c(1,1), bty = 'n')

# study 2 baseline hazard graph, right censoring = 47.5%, interval censoring = 47.5%
# n = 1000
scen2.r47.n1000.min <- max(sapply(ic.scen2.r47.n1000, function(a)
  min(na.omit(c(a$time, a$time2)))))
scen2.r47.n1000.max <- min(sapply(ic.scen2.r47.n1000, function(a)
  max(na.omit(c(a$time, a$time2)))))
scen2.plot.t.r47.n1000 <- seq(scen2.r47.n1000.min, scen2.r47.n1000.max, length.out=1000)
scen2.plot.bh.r47.n1000 <- scen2.plot.bh.r47.r2.n1000 <- list()
for(sim in 1:1000){
  scen2.plot.bh.r47.n1000[[sim]] <- pred.bh(ic.scen2.r47.n1000.results[sim,],
                                      scen2.plot.t.r47.n1000, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  scen2.plot.bh.r47.r2.n1000[[sim]] <- pred.bh(ic.scen2.r47.n1000.results[sim,],
                                         scen2.plot.t.r47.n1000, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}
scen2.plot.bh.r47.all.n1000 <- lapply(scen2.plot.bh.r47.n1000, function(a) a$pred.h0r)
scen2.plot.bh.r47.h0r.mean.n1000 <- colMeans(Reduce("rbind", scen2.plot.bh.r47.all.n1000))
scen2.plot.bh.r47.all.lower.n1000 <- lapply(scen2.plot.bh.r47.n1000, 
                                      function(a) a$pred.h0r.lower)
scen2.plot.bh.r47.h0r.mean.lower.n1000 <- colMeans(Reduce("rbind",
                                                   scen2.plot.bh.r47.all.lower.n1000))
scen2.plot.bh.r47.all.upper.n1000 <- lapply(scen2.plot.bh.r47.n1000, 
                                      function(a) a$pred.h0r.upper)
scen2.plot.bh.r47.h0r.mean.upper.n1000 <- colMeans(Reduce("rbind", 
                                                   scen2.plot.bh.r47.all.upper.n1000))
bh.cp.n1000 <- colSums(Reduce("rbind", lapply(scen2.plot.bh.r47.n1000, 
                                        function(a) ifelse((a$pred.h0r.lower < a$true.h0r)
                                                           & (a$true.h0r < a$pred.h0r.upper),1,0)))) / 1000


scen2.plot.bh.r47.r2.all.n1000 <- lapply(scen2.plot.bh.r47.r2.n1000, function(a) a$pred.h0r)
scen2.plot.bh.r47.r2.h0r.mean.n1000 <- colMeans(Reduce("rbind", 
                                                 scen2.plot.bh.r47.r2.all.n1000))
scen2.plot.bh.r47.r2.all.lower.n1000 <- lapply(scen2.plot.bh.r47.r2.n1000, 
                                         function(a) a$pred.h0r.lower)
scen2.plot.bh.r47.r2.h0r.mean.lower.n1000 <- colMeans(Reduce("rbind", 
                                                      scen2.plot.bh.r47.r2.all.lower.n1000))
scen2.plot.bh.r47.r2.all.upper.n1000 <- lapply(scen2.plot.bh.r47.r2.n1000, 
                                         function(a) a$pred.h0r.upper)
scen2.plot.bh.r47.r2.h0r.mean.upper.n1000 <- colMeans(Reduce("rbind", 
                                                      scen2.plot.bh.r47.r2.all.upper.n1000))
bh.cp.r2.n1000 <- lapply(scen2.plot.bh.r47.r2.n1000, function(a) 
  ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0))
bh.cp.r2.n1000 <- colSums(Reduce("rbind", lapply(scen2.plot.bh.r47.r2.n1000, 
         function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & 
         (a$true.h0r < a$pred.h0r.upper),1,0)))) / 1000

par(mfrow=c(2,2))
plot(scen2.plot.t.r47.n1000, scen2.plot.bh.r47.h0r.mean.n1000, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(A) Risk 1")
lines(scen2.plot.t.r47.n1000, scen2.plot.bh.r47.h0r.mean.lower.n1000, col='#003366',
      lty='dashed')
lines(scen2.plot.t.r47.n1000, scen2.plot.bh.r47.h0r.mean.upper.n1000, col='#003366',
      lty='dashed')
lines(scen2.plot.t.r47.n1000, scen2.plot.bh.r47.n1000[[1]]$true.h0r, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

plot(scen2.plot.t.r47.n1000, scen2.plot.bh.r47.r2.h0r.mean.n1000, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(B) Risk 2")
lines(scen2.plot.t.r47.n1000, scen2.plot.bh.r47.r2.h0r.mean.lower.n1000, col='#003366',
      lty='dashed')
lines(scen2.plot.t.r47.n1000, scen2.plot.bh.r47.r2.h0r.mean.upper.n1000, col='#003366',
      lty='dashed')
lines(scen2.plot.t.r47.n1000, scen2.plot.bh.r47.r2.n1000[[1]]$true.h0r, 
      col = "black",cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')


plot(scen2.plot.t.r47.n1000, bh.cp.r2.n1000, ylab="Coverage probability", xlab="t", 
     col='#003366', type='l', ylim=c(0,1), main="(C) Risk 1")
lines(scen2.plot.t.r47.n1000, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), 
       lty=c(1,2), bty = 'n')
plot(scen2.plot.t.r47.n1000, bh.cp.r2, ylab="Coverage probability", xlab="t", 
     col='#003366', type='l', ylim=c(0,1), main="(D) Risk 2")
lines(scen2.plot.t.r47.n1000, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), 
       lty=c(1,2), bty = 'n')


# CIF for study 2,
# int/left cens perc = 47.5%, right cens = 47.5%, n = 1000
scen2.plot.cif.r47.n1000 <- scen2.plot.r2.r47.r2.n1000 <- list()
for(sim in 1:1000){
  scen2.plot.cif.r47.n1000[[sim]] <- pred.CIF(ic.scen2.r47.n1000.results[sim,],
                                        scen2.plot.t.r47.n1000, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  
  scen2.plot.r2.r47.r2.n1000[[sim]] <- pred.CIF(ic.scen2.r47.n1000.results[sim,],
                                          scen2.plot.t.r47.n1000, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}
scen2.plot.cif.r47.all.n1000 <- lapply(scen2.plot.cif.r47.n1000, 
                                 function(a) a$pred.F0r)
scen2.plot.cif.r47.F0r.mean.n1000 <- colMeans(Reduce("rbind",
                                               scen2.plot.cif.r47.all.n1000))
scen2.plot.cif.r47.all.lower.n1000 <- lapply(scen2.plot.cif.r47.n1000, 
                                       function(a) a$pred.F0r.lower)
scen2.plot.cif.r47.F0r.mean.lower.n1000 <- colMeans(Reduce("rbind", 
                                                    scen2.plot.cif.r47.all.lower.n1000))
scen2.plot.cif.r47.all.upper.n1000 <- lapply(scen2.plot.cif.r47.n1000, 
                                       function(a) a$pred.F0r.upper)
scen2.plot.cif.r47.F0r.mean.upper.n1000 <- colMeans(Reduce("rbind", 
                                                    scen2.plot.cif.r47.all.upper.n1000))

scen2.plot.cif.r47.r2.all.n1000 <- lapply(scen2.plot.r2.r47.r2.n1000, 
                                    function(a) a$pred.F0r)
scen2.plot.cif.r47.r2.F0r.mean.n1000 <- colMeans(Reduce("rbind", 
                                                  scen2.plot.cif.r47.r2.all.n1000))
scen2.plot.cif.r47.r2.all.lower.n1000 <- lapply(scen2.plot.r2.r47.r2.n1000, 
                                          function(a) a$pred.F0r.lower)
scen2.plot.cif.r47.r2.F0r.mean.lower.n1000 <- colMeans(Reduce("rbind", 
                                                       scen2.plot.cif.r47.r2.all.lower.n1000))
scen2.plot.cif.r47.r2.all.upper.n1000 <- lapply(scen2.plot.r2.r47.r2.n1000, 
                                          function(a) a$pred.F0r.upper)
scen2.plot.cif.r47.r2.F0r.mean.upper.n1000 <- colMeans(Reduce("rbind", 
                                                       scen2.plot.cif.r47.r2.all.upper.n1000))

integrand.F0r <- function(x, lambda, lambda1, lambda2, rho){
  (lambda*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho))
}
true.F0r.r1.n1000 <- sapply(scen2.plot.t.r47.n1000, function(a) integrate(integrand.F0r, 0, 
             a, rho=3, lambda=1, lambda1=1, lambda2=0.5)$value)
true.F0r.r2.n1000 <- sapply(scen2.plot.t.r47.n1000, function(a) integrate(integrand.F0r, 0, 
              a, rho=3, lambda=0.5, lambda1=1, lambda2=0.5)$value)
F0r.plot.mat.n1000 <- sapply(scen2.plot.cif.r47.n1000, function(a) a$pred.F0r)
F0r.plot.mean.n1000 <- apply(F0r.plot.mat.n1000, 1, mean)
F0r.plot.r2.mat.n1000 <- sapply(scen2.plot.r2.r47.r2.n1000, function(a) a$pred.F0r)
F0r.plot.r2.mean.n1000 <- apply(F0r.plot.r2.mat.n1000, 1, mean)

par(mfrow=c(2,1))
plot(scen2.plot.t.r47.n1000, scen2.plot.cif.r47.F0r.mean.n1000, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), 
     main="(A) Risk 1", lty=1)
lines(scen2.plot.t.r47.n1000, scen2.plot.cif.r47.F0r.mean.lower.n1000, col='#003366',
      lty=2)
lines(scen2.plot.t.r47.n1000, scen2.plot.cif.r47.F0r.mean.upper.n1000, col='#003366',
      lty=2)
lines(scen2.plot.t.r47.n1000, true.F0r.r1.n1000, col = "black", cex=10, lty=1)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=c(1,1), bty = 'n')

plot(scen2.plot.t.r47.n1000, scen2.plot.cif.r47.r2.F0r.mean.n1000, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), 
     main="(B) Risk 2", lty=1)
lines(scen2.plot.t.r47.n1000, scen2.plot.cif.r47.r2.F0r.mean.lower.n1000, col='#003366',
      lty=2)
lines(scen2.plot.t.r47.n1000, scen2.plot.cif.r47.r2.F0r.mean.upper.n1000, col='#003366',
      lty=2)
lines(scen2.plot.t.r47.n1000, true.F0r.r2.n1000, col = "black", cex=10, lty=1)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=c(1,1), bty = 'n')

# study 2 baseline hazard graph, right censoring = 20%, interval censoring = 75%
# n = 200
scen2.r20.n200.min <- max(sapply(ic.scen2.r20.n200, function(a)
  min(na.omit(c(a$time, a$time2)))))
scen2.r20.n200.max <- min(sapply(ic.scen2.r20.n200, function(a)
  max(na.omit(c(a$time, a$time2)))))
scen2.plot.t.r20 <- seq(scen2.r20.n200.min, scen2.r20.n200.max, length.out=1000)
scen2.plot.bh.r20 <- scen2.plot.bh.r20.r2 <- list()
for(sim in 1:1000){
  scen2.plot.bh.r20[[sim]] <- pred.bh(ic.scen2.r20.n200.results[sim,],
                                      scen2.plot.t.r20, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  scen2.plot.bh.r20.r2[[sim]] <- pred.bh(ic.scen2.r20.n200.results[sim,],
                                         scen2.plot.t.r20, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}
scen2.plot.bh.r20.all <- lapply(scen2.plot.bh.r20, function(a) a$pred.h0r)
scen2.plot.bh.r20.h0r.mean <- colMeans(Reduce("rbind", scen2.plot.bh.r20.all))
scen2.plot.bh.r20.all.lower <- lapply(scen2.plot.bh.r20, 
                                      function(a) a$pred.h0r.lower)
scen2.plot.bh.r20.h0r.mean.lower<- colMeans(Reduce("rbind",
                                                   scen2.plot.bh.r20.all.lower))
scen2.plot.bh.r20.all.upper <- lapply(scen2.plot.bh.r20, 
                                      function(a) a$pred.h0r.upper)
scen2.plot.bh.r20.h0r.mean.upper<- colMeans(Reduce("rbind", 
                                                   scen2.plot.bh.r20.all.upper))
bh.cp <- colSums(Reduce("rbind", lapply(scen2.plot.bh.r20, 
                                        function(a) ifelse((a$pred.h0r.lower < a$true.h0r)
                                                           & (a$true.h0r < a$pred.h0r.upper),1,0)))) / 1000


scen2.plot.bh.r20.r2.all <- lapply(scen2.plot.bh.r20.r2, function(a) a$pred.h0r)
scen2.plot.bh.r20.r2.h0r.mean <- colMeans(Reduce("rbind", 
                                                 scen2.plot.bh.r20.r2.all))
scen2.plot.bh.r20.r2.all.lower <- lapply(scen2.plot.bh.r20.r2, 
                                         function(a) a$pred.h0r.lower)
scen2.plot.bh.r20.r2.h0r.mean.lower<- colMeans(Reduce("rbind", 
                                                      scen2.plot.bh.r20.r2.all.lower))
scen2.plot.bh.r20.r2.all.upper <- lapply(scen2.plot.bh.r20.r2, 
                                         function(a) a$pred.h0r.upper)
scen2.plot.bh.r20.r2.h0r.mean.upper<- colMeans(Reduce("rbind", 
                                                      scen2.plot.bh.r20.r2.all.upper))
bh.cp.r2 <- lapply(scen2.plot.bh.r20.r2, function(a) 
  ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0))
bh.cp.r2 <- colSums(Reduce("rbind", lapply(scen2.plot.bh.r20.r2, 
                                           function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & 
                                                                (a$true.h0r < a$pred.h0r.upper),1,0)))) / 1000

par(mfrow=c(2,2))
plot(scen2.plot.t.r20, scen2.plot.bh.r20.h0r.mean, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(A) Risk 1")
lines(scen2.plot.t.r20, scen2.plot.bh.r20.h0r.mean.lower, col='#003366',
      lty='dashed')
lines(scen2.plot.t.r20, scen2.plot.bh.r20.h0r.mean.upper, col='#003366',
      lty='dashed')
lines(scen2.plot.t.r20, scen2.plot.bh.r20[[1]]$true.h0r, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

plot(scen2.plot.t.r20, scen2.plot.bh.r20.r2.h0r.mean, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(B) Risk 2")
lines(scen2.plot.t.r20, scen2.plot.bh.r20.r2.h0r.mean.lower, col='#003366',
      lty='dashed')
lines(scen2.plot.t.r20, scen2.plot.bh.r20.r2.h0r.mean.upper, col='#003366',
      lty='dashed')
lines(scen2.plot.t.r20, scen2.plot.bh.r20.r2[[1]]$true.h0r, 
      col = "black",cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')


plot(scen2.plot.t.r20, bh.cp.r2, ylab="Coverage probability", xlab="t", 
     col='#003366', type='l', ylim=c(0,1), main="(C) Risk 1")
lines(scen2.plot.t.r20, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), 
       lty=c(1,2), bty = 'n')
plot(scen2.plot.t.r20, bh.cp.r2, ylab="Coverage probability", xlab="t", 
     col='#003366', type='l', ylim=c(0,1), main="(D) Risk 2")
lines(scen2.plot.t.r20, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), 
       lty=c(1,2), bty = 'n')

# CIF for study 2,
# int/left cens perc = 75%, right cens = 20%, n = 200
scen2.plot.cif.r20 <- scen2.plot.r2.r20.r2 <- list()
for(sim in 1:1000){
  scen2.plot.cif.r20[[sim]] <- pred.CIF(ic.scen2.r20.n200.results[sim,],
                                        scen2.plot.t.r20, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  
  scen2.plot.r2.r20.r2[[sim]] <- pred.CIF(ic.scen2.r20.n200.results[sim,],
                                          scen2.plot.t.r20, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}
scen2.plot.cif.r20.all <- lapply(scen2.plot.cif.r20, 
                                 function(a) a$pred.F0r)
scen2.plot.cif.r20.F0r.mean <- colMeans(Reduce("rbind",
                                               scen2.plot.cif.r20.all))
scen2.plot.cif.r20.all.lower <- lapply(scen2.plot.cif.r20, 
                                       function(a) a$pred.F0r.lower)
scen2.plot.cif.r20.F0r.mean.lower<- colMeans(Reduce("rbind", 
                                                    scen2.plot.cif.r20.all.lower))
scen2.plot.cif.r20.all.upper <- lapply(scen2.plot.cif.r20, 
                                       function(a) a$pred.F0r.upper)
scen2.plot.cif.r20.F0r.mean.upper<- colMeans(Reduce("rbind", 
                                                    scen2.plot.cif.r20.all.upper))

scen2.plot.cif.r20.r2.all <- lapply(scen2.plot.r2.r20.r2, 
                                    function(a) a$pred.F0r)
scen2.plot.cif.r20.r2.F0r.mean <- colMeans(Reduce("rbind", 
                                                  scen2.plot.cif.r20.r2.all))
scen2.plot.cif.r20.r2.all.lower <- lapply(scen2.plot.r2.r20.r2, 
                                          function(a) a$pred.F0r.lower)
scen2.plot.cif.r20.r2.F0r.mean.lower<- colMeans(Reduce("rbind", 
                                                       scen2.plot.cif.r20.r2.all.lower))
scen2.plot.cif.r20.r2.all.upper <- lapply(scen2.plot.r2.r20.r2, 
                                          function(a) a$pred.F0r.upper)
scen2.plot.cif.r20.r2.F0r.mean.upper<- colMeans(Reduce("rbind", 
                                                       scen2.plot.cif.r20.r2.all.upper))

integrand.F0r <- function(x, lambda, lambda1, lambda2, rho){
  (lambda*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho))
}
true.F0r.r1 <- sapply(scen2.plot.t.r20, function(a) integrate(integrand.F0r, 0, 
                                                              a, rho=3, lambda=1, lambda1=1, lambda2=0.5)$value)
true.F0r.r2 <- sapply(scen2.plot.t.r20, function(a) integrate(integrand.F0r, 0, 
                                                              a, rho=3, lambda=0.5, lambda1=1, lambda2=0.5)$value)
F0r.plot.mat <- sapply(scen2.plot.cif.r20, function(a) a$pred.F0r)
F0r.plot.mean <- apply(F0r.plot.mat, 1, mean)
F0r.plot.r2.mat <- sapply(scen2.plot.r2.r20.r2, function(a) a$pred.F0r)
F0r.plot.r2.mean <- apply(F0r.plot.r2.mat, 1, mean)

par(mfrow=c(2,1))
plot(scen2.plot.t.r20, scen2.plot.cif.r20.F0r.mean, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), 
     main="(A) Risk 1", lty=1)
lines(scen2.plot.t.r20, scen2.plot.cif.r20.F0r.mean.lower, col='#003366',
      lty=2)
lines(scen2.plot.t.r20, scen2.plot.cif.r20.F0r.mean.upper, col='#003366',
      lty=2)
lines(scen2.plot.t.r20, true.F0r.r1, col = "black", cex=10, lty=1)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=c(1,1), bty = 'n')

plot(scen2.plot.t.r20, scen2.plot.cif.r20.r2.F0r.mean, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), 
     main="(B) Risk 2", lty=1)
lines(scen2.plot.t.r20, scen2.plot.cif.r20.r2.F0r.mean.lower, col='#003366',
      lty=2)
lines(scen2.plot.t.r20, scen2.plot.cif.r20.r2.F0r.mean.upper, col='#003366',
      lty=2)
lines(scen2.plot.t.r20, true.F0r.r2, col = "black", cex=10, lty=1)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=c(1,1), bty = 'n')

#study 2 baseline hazard graph, r20, n1000
scen2.r20.n1000.min <- max(sapply(ic.scen2.r20.n1000, function(a)
  min(na.omit(c(a$time, a$time2)))))
scen2.r20.n1000.max <- min(sapply(ic.scen2.r20.n1000, function(a)
  max(na.omit(c(a$time, a$time2)))))
scen2.plot.t.r20.n1000 <- seq(scen2.r20.n1000.min, scen2.r20.n1000.max, length.out=1000)
scen2.plot.bh.r20.n1000 <- scen2.plot.bh.r20.r2.n1000 <- list()
for(sim in 1:1000){
  scen2.plot.bh.r20.n1000[[sim]] <- pred.bh(ic.scen2.r20.n1000.results[sim,],
                                      scen2.plot.t.r20.n1000, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  scen2.plot.bh.r20.r2.n1000[[sim]] <- pred.bh(ic.scen2.r20.n1000.results[sim,],
                                         scen2.plot.t.r20.n1000, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}
scen2.plot.bh.r20.all.n1000 <- lapply(scen2.plot.bh.r20.n1000, function(a) a$pred.h0r)
scen2.plot.bh.r20.h0r.mean.n1000 <- colMeans(Reduce("rbind", scen2.plot.bh.r20.all.n1000))
scen2.plot.bh.r20.all.lower.n1000 <- lapply(scen2.plot.bh.r20.n1000, 
                                      function(a) a$pred.h0r.lower)
scen2.plot.bh.r20.h0r.mean.lower.n1000 <- colMeans(Reduce("rbind",
                                                   scen2.plot.bh.r20.all.lower.n1000))
scen2.plot.bh.r20.all.upper.n1000 <- lapply(scen2.plot.bh.r20.n1000, 
                                      function(a) a$pred.h0r.upper)
scen2.plot.bh.r20.h0r.mean.upper.n1000 <- colMeans(Reduce("rbind", 
                                                   scen2.plot.bh.r20.all.upper.n1000))
bh.cp.n1000 <- colSums(Reduce("rbind", lapply(scen2.plot.bh.r20.n1000, 
                                        function(a) ifelse((a$pred.h0r.lower < a$true.h0r)
                                                           & (a$true.h0r < a$pred.h0r.upper),1,0)))) / 1000


scen2.plot.bh.r20.r2.all.n1000 <- lapply(scen2.plot.bh.r20.r2.n1000, function(a) a$pred.h0r)
scen2.plot.bh.r20.r2.h0r.mean.n1000 <- colMeans(Reduce("rbind", 
                                                 scen2.plot.bh.r20.r2.all.n1000))
scen2.plot.bh.r20.r2.all.lower.n1000 <- lapply(scen2.plot.bh.r20.r2.n1000, 
                                         function(a) a$pred.h0r.lower)
scen2.plot.bh.r20.r2.h0r.mean.lower.n1000 <- colMeans(Reduce("rbind", 
                                                      scen2.plot.bh.r20.r2.all.lower.n1000))
scen2.plot.bh.r20.r2.all.upper.n1000 <- lapply(scen2.plot.bh.r20.r2.n1000, 
                                         function(a) a$pred.h0r.upper)
scen2.plot.bh.r20.r2.h0r.mean.upper.n1000 <- colMeans(Reduce("rbind", 
                                                      scen2.plot.bh.r20.r2.all.upper.n1000))
bh.cp.r2.n1000 <- lapply(scen2.plot.bh.r20.r2.n1000, function(a) 
  ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0))
bh.cp.r2.n1000 <- colSums(Reduce("rbind", lapply(scen2.plot.bh.r20.r2.n1000, 
                                           function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & 
                                                                (a$true.h0r < a$pred.h0r.upper),1,0)))) / 1000

par(mfrow=c(2,2))
plot(scen2.plot.t.r20.n1000, scen2.plot.bh.r20.h0r.mean.n1000, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(A) Risk 1")
lines(scen2.plot.t.r20.n1000, scen2.plot.bh.r20.h0r.mean.lower.n1000, col='#003366',
      lty='dashed')
lines(scen2.plot.t.r20.n1000, scen2.plot.bh.r20.h0r.mean.upper.n1000, col='#003366',
      lty='dashed')
lines(scen2.plot.t.r20.n1000, scen2.plot.bh.r20.n1000[[1]]$true.h0r, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

plot(scen2.plot.t.r20.n1000, scen2.plot.bh.r20.r2.h0r.mean.n1000, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(B) Risk 2")
lines(scen2.plot.t.r20.n1000, scen2.plot.bh.r20.r2.h0r.mean.lower.n1000, col='#003366',
      lty='dashed')
lines(scen2.plot.t.r20.n1000, scen2.plot.bh.r20.r2.h0r.mean.upper.n1000, col='#003366',
      lty='dashed')
lines(scen2.plot.t.r20.n1000, scen2.plot.bh.r20.r2.n1000[[1]]$true.h0r, 
      col = "black",cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')


plot(scen2.plot.t.r20.n1000, bh.cp.r2.n1000, ylab="Coverage probability", xlab="t", 
     col='#003366', type='l', ylim=c(0,1), main="(C) Risk 1")
lines(scen2.plot.t.r20.n1000, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), 
       lty=c(1,2), bty = 'n')
plot(scen2.plot.t.r20, bh.cp.r2.n1000, ylab="Coverage probability", xlab="t", 
     col='#003366', type='l', ylim=c(0,1), main="(D) Risk 2")
lines(scen2.plot.t.r20.n1000, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), 
       lty=c(1,2), bty = 'n')

# study 2 baseline hazard
# int/left cens perc = 75%, right cens = 20%, n = 1000
scen2.plot.cif.r20.n1000 <- scen2.plot.r2.r20.r2.n1000 <- list()
for(sim in 1:1000){
  scen2.plot.cif.r20.n1000[[sim]] <- pred.CIF(ic.scen2.r20.n1000.results[sim,],
                                        scen2.plot.t.r20.n1000, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  
  scen2.plot.r2.r20.r2.n1000[[sim]] <- pred.CIF(ic.scen2.r20.n1000.results[sim,],
                                          scen2.plot.t.r20.n1000, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}
scen2.plot.cif.r20.all.n1000 <- lapply(scen2.plot.cif.r20.n1000, 
                                 function(a) a$pred.F0r)
scen2.plot.cif.r20.F0r.mean.n1000 <- colMeans(Reduce("rbind",
                                               scen2.plot.cif.r20.all.n1000))
scen2.plot.cif.r20.all.lower.n1000 <- lapply(scen2.plot.cif.r20.n1000, 
                                       function(a) a$pred.F0r.lower)
scen2.plot.cif.r20.F0r.mean.lower.n1000 <- colMeans(Reduce("rbind", 
                                                    scen2.plot.cif.r20.all.lower.n1000))
scen2.plot.cif.r20.all.upper.n1000 <- lapply(scen2.plot.cif.r20.n1000, 
                                       function(a) a$pred.F0r.upper)
scen2.plot.cif.r20.F0r.mean.upper.n1000 <- colMeans(Reduce("rbind", 
                                                    scen2.plot.cif.r20.all.upper.n1000))

scen2.plot.cif.r20.r2.all.n1000 <- lapply(scen2.plot.r2.r20.r2.n1000, 
                                    function(a) a$pred.F0r)
scen2.plot.cif.r20.r2.F0r.mean.n1000 <- colMeans(Reduce("rbind", 
                                                  scen2.plot.cif.r20.r2.all.n1000))
scen2.plot.cif.r20.r2.all.lower.n1000 <- lapply(scen2.plot.r2.r20.r2.n1000, 
                                          function(a) a$pred.F0r.lower)
scen2.plot.cif.r20.r2.F0r.mean.lower.n1000 <- colMeans(Reduce("rbind", 
                                                       scen2.plot.cif.r20.r2.all.lower.n1000))
scen2.plot.cif.r20.r2.all.upper.n1000 <- lapply(scen2.plot.r2.r20.r2.n1000, 
                                          function(a) a$pred.F0r.upper)
scen2.plot.cif.r20.r2.F0r.mean.upper.n1000 <- colMeans(Reduce("rbind", 
                                                       scen2.plot.cif.r20.r2.all.upper.n1000))

integrand.F0r <- function(x, lambda, lambda1, lambda2, rho){
  (lambda*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho))
}
true.F0r.r1.n1000 <- sapply(scen2.plot.t.r20.n1000, function(a) integrate(integrand.F0r, 0, 
                                                              a, rho=3, lambda=1, lambda1=1, lambda2=0.5)$value)
true.F0r.r2.n1000 <- sapply(scen2.plot.t.r20.n1000, function(a) integrate(integrand.F0r, 0, 
                                                              a, rho=3, lambda=0.5, lambda1=1, lambda2=0.5)$value)
F0r.plot.mat.n1000 <- sapply(scen2.plot.cif.r20.n1000, function(a) a$pred.F0r)
F0r.plot.mean.n1000 <- apply(F0r.plot.mat.n1000, 1, mean)
F0r.plot.r2.mat.n1000 <- sapply(scen2.plot.r2.r20.r2.n1000, function(a) a$pred.F0r)
F0r.plot.r2.mean.n1000 <- apply(F0r.plot.r2.mat.n1000, 1, mean)

par(mfrow=c(2,1))
plot(scen2.plot.t.r20.n1000, scen2.plot.cif.r20.F0r.mean.n1000, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), 
     main="(A) Risk 1", lty=1)
lines(scen2.plot.t.r20.n1000, scen2.plot.cif.r20.F0r.mean.lower.n1000, col='#003366',
      lty=2)
lines(scen2.plot.t.r20.n1000, scen2.plot.cif.r20.F0r.mean.upper.n1000, col='#003366',
      lty=2)
lines(scen2.plot.t.r20.n1000, true.F0r.r1.n1000, col = "black", cex=10, lty=1)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=c(1,1), bty = 'n')

plot(scen2.plot.t.r20.n1000, scen2.plot.cif.r20.r2.F0r.mean.n1000, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), 
     main="(B) Risk 2", lty=1)
lines(scen2.plot.t.r20.n1000, scen2.plot.cif.r20.r2.F0r.mean.lower.n1000, col='#003366',
      lty=2)
lines(scen2.plot.t.r20.n1000, scen2.plot.cif.r20.r2.F0r.mean.upper.n1000, col='#003366',
      lty=2)
lines(scen2.plot.t.r20.n1000, true.F0r.r2.n1000, col = "black", cex=10, lty=1)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=c(1,1), bty = 'n')






#     Cure fraction simulations      #######################

#Function to generate cure fraction data
phcshcf.gen.data <- function(no.of.sims=1000, cureMin = 0.5, cureMax = 4.5,
                             pi.r.z.coef = c(3, 1, -0.5), n.obs = 1000,
                             gammaL = 0.5, gammaR = 2.5, beta = list(c(0.5,-2), c(1,-1)),
                             prob_event){
  simdata <- list()
  for(s in 1:no.of.sims){
    lambda = c(1, 0.5)
    rho = 3
    seed <- s
    set.seed(seed)
    #Generate probability of cured
    z1 <- rnorm(n.obs)
    z2 <- rnorm(n.obs)
    pi.U <- runif(n.obs)
    pi.z = exp(pi.r.z.coef[1] + pi.r.z.coef[2]*z1 + pi.r.z.coef[3]*z2) /
      (1 + exp(pi.r.z.coef[1] + pi.r.z.coef[2]*z1 + pi.r.z.coef[3]*z2))
    noncure <- pi.U < pi.z
    #Make Z and Z covariates the same
    #x1 <- z1
    #x2 <- z2
    #make Z and X difference
    x1 <- rnorm(n.obs)
    x2 <- rnorm(n.obs)
    z <- data.matrix(data.frame(z1,z2,noncure))
    x <- data.matrix(data.frame(x1,x2,noncure))
    #split z and x covariates into cured and non-cured
    z.noncured <- z[noncure==1,]
    z.cured <- z[noncure==0,]
    n.cured <- nrow(z.cured)
    n.noncured <- nrow(z.noncured)
    x.noncured <- x[noncure==1,]
    x.cured <- x[noncure==0,]
    #for cured, generate censored observation times
    Ul.cured <- runif(n.cured)
    Ur.cured <- runif(n.cured,Ul.cured,1)
    # Uc <- gammaR*Ur.cured
    #Uc <- rexp(n.cured, 1/2)
    Uc <- runif(n.cured, cureMin, cureMax)
    Uc.time2 <- NA
    risk <- NA
    event <- 0
    cured <- cbind(z.cured, "time"=Uc, "time2"=Uc.time2, risk, event, 
                   x.cured[, -ncol(x.cured)])
    #for noncured, generate observation times
    u.noncured <- runif(n.noncured)
    t.noncured <- (-log(u.noncured) / 
                     (lambda[1]*exp(x.noncured[, -ncol(x.noncured)] %*% beta[[1]])
                      + lambda[2]*exp(x.noncured[, -ncol(x.noncured)] %*% beta[[2]])))^(1/rho)
    Ul.noncured <- runif(n.noncured)
    Ur.noncured <- runif(n.noncured,Ul.noncured,1)
    Ue.noncured <- runif(n.noncured)
    event <- ifelse(Ue.noncured < prob_event,1,
                    ifelse(gammaL*Ul.noncured <= t.noncured &
                             t.noncured <= gammaR*Ur.noncured,3,ifelse(t.noncured < gammaL*Ul.noncured,
                                                                       2,0)))
    p1 <- (lambda[1]*exp(x.noncured[, -ncol(x.noncured)] %*% beta[[1]])) / 
      ((lambda[1]*exp(x.noncured[, -ncol(x.noncured)] %*% beta[[1]]) +
          lambda[2]*exp(x.noncured[, -ncol(x.noncured)] %*% beta[[2]])))
    p2 <- (lambda[2]*exp(x.noncured[, -ncol(x.noncured)] %*% beta[[2]])) / 
      ((lambda[1]*exp(x.noncured[, -ncol(x.noncured)] %*% beta[[1]]) +
          lambda[2]*exp(x.noncured[, -ncol(x.noncured)] %*% beta[[2]])))
    v <- runif(n.noncured)
    time1 <- ifelse(event==0,gammaR*Ur.noncured,ifelse(event==1,t.noncured,
                                                       ifelse(event==2,gammaL*Ul.noncured,gammaL*Ul.noncured)))
    time2 <- ifelse(event==0,NA,ifelse(event==1,t.noncured,ifelse(event==2,NA,
                                                                  gammaR*Ur.noncured)))
    risk.noncured <- ifelse(event==0,NA,ifelse(v <= p1, 1, 2))
    noncured <- data.frame(cbind(z.noncured, "time"=time1, time2, 
                                 "risk"=risk.noncured, event, x.noncured[, -ncol(x.noncured)]))
    #combine cured and noncured data
    simdata[[s]] <- data.frame(rbind(cured, noncured))
    simdata[[s]]$tmid <- ifelse(simdata[[s]]$event==0, simdata[[s]]$time,
                                ifelse(simdata[[s]]$event==1, simdata[[s]]$time,
                                       ifelse(simdata[[s]]$event==2, simdata[[s]]$time/2,
                                              (simdata[[s]]$time +simdata[[s]]$time2) / 2)))
    simdata[[s]]$risk_1 <- ifelse(is.na(simdata[[s]]$risk), 0,
                                  ifelse(simdata[[s]]$risk==1, 1, 0))
    simdata[[s]]$risk_2 <- ifelse(is.na(simdata[[s]]$risk), 0,
                                  ifelse(simdata[[s]]$risk==2, 1, 0))
  }
  simdata
}

#summarize cure fraction simulation results
cf.mpl.results <- function(object, true.gamma, true.beta){
  #calculate bias
  all.gamma.list <- all.beta1.list <- all.beta2.list <- list()
  for(i in 1:n.sim){
    if(object[i,]$valid==1){
      all.gamma.list[[i]] <- object[i,]$gamma
      all.beta1.list[[i]] <- object[i,]$beta[[1]]
      all.beta2.list[[i]] <- object[i,]$beta[[2]]
    }
  }
  all.gamma <- Reduce("rbind", all.gamma.list)
  all.beta1 <- Reduce("rbind", all.beta1.list)
  all.beta2 <- Reduce("rbind", all.beta2.list)
  gamma.bias <- colMeans(all.gamma) - true.gamma
  beta1.bias <- colMeans(all.beta1) - true.beta[[1]]
  beta2.bias <- colMeans(all.beta2) - true.beta[[2]]
  all.gamma.std.list <- all.beta1.std.list <- all.beta2.std.list <- list()
  for(i in 1:n.sim){
    if(object[i,]$valid==1){
      all.gamma.std.list[[i]] <- object[i,]$seG
      all.beta1.std.list[[i]] <- object[i,]$seB[[1]]
      all.beta2.std.list[[i]] <- object[i,]$seB[[2]]
    }
  }
  all.gamma.std <- Reduce("rbind", all.gamma.std.list)
  all.beta1.std <- Reduce("rbind", all.beta1.std.list)
  all.beta2.std <- Reduce("rbind", all.beta2.std.list)
  gamma.std.mc <- apply(all.gamma, 2, sd)
  gamma.std.asymp <- colMeans(all.gamma.std)
  beta1.std.mc <- apply(all.beta1, 2, sd)
  beta1.std.asymp <- colMeans(all.beta1.std)
  beta2.std.mc <- apply(all.beta2, 2, sd)
  beta2.std.asymp <- colMeans(all.beta2.std)
  #coverage probability
  all.gamma.lower.list <- all.gamma.upper.list <- all.gamma.cp.list <- list()
  all.beta1.lower.list <- all.beta1.upper.list <- all.beta1.cp.list <- list()
  all.beta2.lower.list <- all.beta2.upper.list <- all.beta2.cp.list <- list()
  for(i in 1:n.sim){
    if(object[i,]$valid==1){
      all.gamma.lower.list[[i]] <- object[i,]$gamma - 1.96*object[i,]$seG
      all.gamma.upper.list[[i]] <- object[i,]$gamma + 1.96*object[i,]$seG
      all.gamma.cp.list[[i]] <- ifelse(all.gamma.lower.list[[i]] < true.gamma & true.gamma < all.gamma.upper.list[[i]], 1, 0)
      
      all.beta1.lower.list[[i]] <- object[i,]$beta[[1]] - 1.96*object[i,]$seB[[1]]
      all.beta1.upper.list[[i]] <- object[i,]$beta[[1]] + 1.96*object[i,]$seB[[1]]
      all.beta1.cp.list[[i]] <- ifelse(all.beta1.lower.list[[i]] < true.beta[[1]] & true.beta[[1]] < all.beta1.upper.list[[i]], 1, 0)
      
      all.beta2.lower.list[[i]] <- object[i,]$beta[[2]] - 1.96*object[i,]$seB[[2]]
      all.beta2.upper.list[[i]] <- object[i,]$beta[[2]] + 1.96*object[i,]$seB[[2]]
      all.beta2.cp.list[[i]] <- ifelse(all.beta2.lower.list[[i]] < true.beta[[2]] & true.beta[[2]] < all.beta2.upper.list[[i]], 1, 0)
    }
  }
  all.gamma.lower <- Reduce("rbind", all.gamma.lower.list)
  all.gamma.upper <- Reduce("rbind", all.gamma.upper.list)
  all.gamma.cp <- Reduce("rbind", all.gamma.cp.list)
  all.beta1.lower <- Reduce("rbind", all.beta1.lower.list)
  all.beta1.upper <- Reduce("rbind", all.beta1.upper.list)
  all.beta1.cp <- Reduce("rbind", all.beta1.cp.list)
  all.beta2.lower <- Reduce("rbind", all.beta2.lower.list)
  all.beta2.upper <- Reduce("rbind", all.beta2.upper.list)
  all.beta2.cp <- Reduce("rbind", all.beta2.cp.list)
  gamma.cp <- colMeans(all.gamma.cp, na.rm=TRUE)
  beta1.cp <- colMeans(all.beta1.cp, na.rm=TRUE)
  beta2.cp <- colMeans(all.beta2.cp, na.rm=TRUE)
  gamma.list <- list("bias"=gamma.bias, "std.asymp"=gamma.std.asymp, 
                     "std.mc"=gamma.std.mc, "gamma.cp"=gamma.cp)
  beta1.list <- list("bias"=beta1.bias, "std.asymp"=beta1.std.asymp, 
                     "std.mc"=beta1.std.mc, "cp"=beta1.cp)
  beta2.list <- list("bias"=beta2.bias, "std.asymp"=beta2.std.asymp, 
                     "std.mc"=beta2.std.mc, "cp"=beta2.cp)
  rlist <- list("gamma"=gamma.list, "beta1"=beta1.list, "beta2"=beta2.list)
  return(rlist)
}

#cure fraction proportion = 20%, n = 300
n.sim = 1000
#sim11, n = 300, cm = 12
true.gamma <- c(1.75, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.sim11.data.n300 <- phcshcf.gen.data(no.of.sims=n.sim, cureMin = 0.5, cureMax = 4.378516,
                                       pi.r.z.coef = true.gamma, n.obs = 300, prob_event = 0.2,
                                       beta = true.beta)
total.noncure <- lapply(cf.sim11.data.n300, function(a) a$noncure)
table(unlist(total.noncure)) / length(unlist(total.noncure))

table(unlist(lapply(cf.sim11.data.n300, function(a) a$event))) /
  length(unlist(lapply(cf.sim11.data.n300, function(a) a$event)))

#Use parallel processing
cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cf.sim11n300.results.time <- system.time({cf.sim11n300.results <- foreach(i = 1:n.sim, .combine = rbind, .options.snow = opts) %dopar% {
  zmat <- as.matrix(cbind(cf.sim11.data.n300[[i]]$z1, cf.sim11.data.n300[[i]]$z2))
  colnames(zmat) <- c("z1", "z2")
  phcsh_mpl(Surv(time, time2, event=event, "interval") ~ x1 + x2,
            risk = cf.sim11.data.n300[[i]]$risk, data = cf.sim11.data.n300[[i]],
            z = zmat, max.outer = 10,
            max.iter = 500, lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL,
            aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = FALSE,
            gq.points = 50)
}
})
close(pb)
stopCluster(cl)

true.gamma <- c(1.75, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.mpl.results(cf.sim11n300.results, true.gamma, true.beta)

#cure fraction proportion = 20%, baseline hazard
cf.sim11.data.n300.min <- max(sapply(cf.sim11.data.n300, function(a)
  min(na.omit(c(a$time, a$time2)))))
cf.sim11.data.n300.max <- min(sapply(cf.sim11.data.n300, function(a)
  max(na.omit(c(a$time, a$time2)))))
cf1.scen1.plot.t <- seq(cf.sim11.data.n300.min, cf.sim11.data.n300.max, length.out=1000)
cf1.plot.bh.r1 <- cf1.plot.bh.r2 <- list()
for(sim in 1:1000){
  cf1.plot.bh.r1[[sim]] <- pred.bh(cf.sim11n300.results[sim,],
                                   cf1.scen1.plot.t, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  cf1.plot.bh.r2[[sim]] <- pred.bh(cf.sim11n300.results[sim,],
                                   cf1.scen1.plot.t, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}

temp.valid = vector()
for(i in 1:1000){
  temp.valid[i] = cf.sim11n300.results[i,]$valid
}

cf1.plot.bh.r1.all <- lapply(cf1.plot.bh.r1, function(a) a$pred.h0r)
cf1.plot.bh.r1.h0r.mean <- colMeans(Reduce("rbind", cf1.plot.bh.r1.all))
cf1.plot.bh.r1.all.lower <- lapply(cf1.plot.bh.r1, function(a) a$pred.h0r.lower)
cf1.plot.bh.r1.h0r.mean.lower<- colMeans(Reduce("rbind", cf1.plot.bh.r1.all.lower))
cf1.plot.bh.r1.all.upper <- lapply(cf1.plot.bh.r1, function(a) a$pred.h0r.upper)
cf1.plot.bh.r1.h0r.mean.upper<- colMeans(Reduce("rbind", cf1.plot.bh.r1.all.upper))
#cf1.bh.cp.r1.all <- lapply(cf1.plot.bh.r1, function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0))
cf1.bh.cp.r1 <- colSums(Reduce("rbind", lapply(cf1.plot.bh.r1, function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0)))) / sum(temp.valid)


cf1.plot.bh.r2.all <- lapply(cf1.plot.bh.r2, function(a) a$pred.h0r)
cf1.plot.bh.r2.h0r.mean <- colMeans(Reduce("rbind", cf1.plot.bh.r2.all))
cf1.plot.bh.r2.all.lower <- lapply(cf1.plot.bh.r2, function(a) a$pred.h0r.lower)
cf1.plot.bh.r2.h0r.mean.lower<- colMeans(Reduce("rbind", cf1.plot.bh.r2.all.lower))
cf1.plot.bh.r2.all.upper <- lapply(cf1.plot.bh.r2, function(a) a$pred.h0r.upper)
cf1.plot.bh.r2.h0r.mean.upper<- colMeans(Reduce("rbind", cf1.plot.bh.r2.all.upper))
cf1.bh.cp.r2 <- colSums(Reduce("rbind", lapply(cf1.plot.bh.r2, function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0)))) / sum(temp.valid)

par(mfrow=c(2,2))
plot(cf1.scen1.plot.t, cf1.plot.bh.r1.h0r.mean, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(A) Risk 1", ylim=c(0,250))
lines(cf1.scen1.plot.t, cf1.plot.bh.r1.h0r.mean.lower, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t, cf1.plot.bh.r1.h0r.mean.upper, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t, cf1.plot.bh.r1[[1]]$true.h0r, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

plot(cf1.scen1.plot.t, cf1.plot.bh.r2.h0r.mean, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(B) Risk 2", ylim=c(0,250))
lines(cf1.scen1.plot.t, cf1.plot.bh.r2.h0r.mean.lower, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t, cf1.plot.bh.r2.h0r.mean.upper, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t, cf1.plot.bh.r2[[1]]$true.h0r, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')


plot(cf1.scen1.plot.t, cf1.bh.cp.r1, ylab="Coverage probability", xlab="t", col='#003366',
     type='l', ylim=c(0,1), main="(C) Risk 1")
lines(cf1.scen1.plot.t, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), lty=c(1,2),
       bty = 'n')

plot(cf1.scen1.plot.t, cf1.bh.cp.r2, ylab="Coverage probability", xlab="t", col='#003366',
     type='l', ylim=c(0,1), main="(D) Risk 2")
lines(cf1.scen1.plot.t, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), lty=c(1,2),
       bty = 'n')

#cure fraction proportion = 20%, CIF graph
cf1.plot.cif.r1 <- cf1.plot.cif.r2 <- list()
for(sim in 1:1000){
  if(cf.sim11n300.results[sim,]$valid==1){
    cf1.plot.cif.r1[[sim]] <- pred.CIF(cf.sim11n300.results[sim,],
                                       cf1.scen1.plot.t, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
    
    cf1.plot.cif.r2[[sim]] <- pred.CIF(cf.sim11n300.results[sim,],
                                       cf1.scen1.plot.t, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
  }
}
cf1.plot.cif.r1.all <- lapply(cf1.plot.cif.r1, function(a) a$pred.F0r)
cf1.plot.cif.r1.F0r.mean <- colMeans(Reduce("rbind", cf1.plot.cif.r1.all))
cf1.plot.cif.r1.all.lower <- lapply(cf1.plot.cif.r1, function(a) a$pred.F0r.lower)
cf1.plot.cif.r1.F0r.mean.lower<- colMeans(Reduce("rbind", cf1.plot.cif.r1.all.lower))
cf1.plot.cif.r1.all.upper <- lapply(cf1.plot.cif.r1, function(a) a$pred.F0r.upper)
cf1.plot.cif.r1.F0r.mean.upper<- colMeans(Reduce("rbind", cf1.plot.cif.r1.all.upper))

cf1.plot.cif.r2.all <- lapply(cf1.plot.cif.r2, function(a) a$pred.F0r)
cf1.plot.cif.r2.F0r.mean <- colMeans(Reduce("rbind", cf1.plot.cif.r2.all))
cf1.plot.cif.r2.all.lower <- lapply(cf1.plot.cif.r2, function(a) a$pred.F0r.lower)
cf1.plot.cif.r2.F0r.mean.lower<- colMeans(Reduce("rbind", cf1.plot.cif.r2.all.lower))
cf1.plot.cif.r2.all.upper <- lapply(cf1.plot.cif.r2, function(a) a$pred.F0r.upper)
cf1.plot.cif.r2.F0r.mean.upper<- colMeans(Reduce("rbind", cf1.plot.cif.r2.all.upper))

integrand.F01 <- function(x, lambda1, lambda2, rho){
  (lambda1*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho))
}

integrand.F02 <- function(x, lambda1, lambda2, rho){
  (lambda2*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho))
}
true.F0r.r1 <- sapply(cf1.scen1.plot.t, function(a) integrate(integrand.F01, 0, a, rho=3,
                                                              lambda1=1, lambda2 = 0.5)$value)
true.F0r.r2 <- sapply(cf1.scen1.plot.t, function(a) integrate(integrand.F02, 0, a, rho=3,
                                                              lambda1=1, lambda2=0.5)$value)
F0r.plot.mat <- sapply(cf1.plot.cif.r1, function(a) a$pred.F0r)
#calculate length of each list
lengthr1 <- lengthr2 <- list()
for(sim in 1:1000){
  lengthr1[[sim]] <- length(cf1.plot.cif.r1[[sim]]$pred.F0r)
  lengthr2[[sim]] <- length(cf1.plot.cif.r2[[sim]]$pred.F0r)
}
F0r.plot.mat.list <- list()
count=1
for(sim in 1:1000){
  if(length(cf1.plot.cif.r1[[sim]]$pred.F0r)==1000){
    F0r.plot.mat.list[[count]] <- cf1.plot.cif.r1[[sim]]$pred.F0r
  }
  count=count+1
}
F0r.plot.mat <- Reduce("cbind", F0r.plot.mat.list)
F0r.plot.mean <- apply(F0r.plot.mat, 1, mean)
F0r.plot.l95 <- apply(F0r.plot.mat, 1, quantile, 0.025)
F0r.plot.u95 <- apply(F0r.plot.mat, 1, quantile, 0.975)
F0r.plot.r2.mat.list <- list()
count=1
for(sim in 1:1000){
  if(length(cf1.plot.cif.r2[[sim]]$pred.F0r)==1000){
    F0r.plot.r2.mat.list[[count]] <- cf1.plot.cif.r2[[sim]]$pred.F0r
  }
  count=count+1
}
F0r.plot.r2.mat <- Reduce("cbind", F0r.plot.r2.mat.list)
F0r.plot.r2.mean <- apply(F0r.plot.r2.mat, 1, mean)
F0r.plot.r2.l95 <- apply(F0r.plot.r2.mat, 1, quantile, 0.025)
F0r.plot.r2.u95 <- apply(F0r.plot.r2.mat, 1, quantile, 0.975)

par(mfrow=c(2,1))
plot(cf1.scen1.plot.t, cf1.plot.cif.r1.F0r.mean, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), main="(A) Risk 1", xlim=c(0, 2))
lines(cf1.scen1.plot.t, cf1.plot.cif.r1.F0r.mean.lower, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t, cf1.plot.cif.r1.F0r.mean.upper, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t, true.F0r.r1, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

plot(cf1.scen1.plot.t, cf1.plot.cif.r2.F0r.mean, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), main="(B) Risk 2", xlim=c(0, 2))
lines(cf1.scen1.plot.t, cf1.plot.cif.r2.F0r.mean.lower, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t, cf1.plot.cif.r2.F0r.mean.upper, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t, true.F0r.r2, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

#create variables for Cox Regression application
for(i in 1:n.sim){
  cf.sim11.data.n300[[i]]$tmid <- ifelse(cf.sim11.data.n300[[i]]$event==0, cf.sim11.data.n300[[i]]$time,
                                         ifelse(cf.sim11.data.n300[[i]]$event==1, cf.sim11.data.n300[[i]]$time,
                                                ifelse(cf.sim11.data.n300[[i]]$event==2, cf.sim11.data.n300[[i]]$time/2,
                                                       (cf.sim11.data.n300[[i]]$time +cf.sim11.data.n300[[i]]$time2) / 2)))
  cf.sim11.data.n300[[i]]$risk_1 <- ifelse(is.na(cf.sim11.data.n300[[i]]$risk), 0,
                                           ifelse(cf.sim11.data.n300[[i]]$risk==1, 1, 0))
  cf.sim11.data.n300[[i]]$risk_2 <- ifelse(is.na(cf.sim11.data.n300[[i]]$risk), 0,
                                           ifelse(cf.sim11.data.n300[[i]]$risk==2, 1, 0))
}

cf.sim1.r1.n300.cox = cf.sim1.r2.n300.cox = list()
for(sim in 1:n.sim){
  cf.sim1.r1.n300.cox[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, data = cf.sim11.data.n300[[sim]])
  cf.sim1.r2.n300.cox[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, data = cf.sim11.data.n300[[sim]])
}
sim.results.cox(cf.sim1.r1.n300.cox, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(cf.sim1.r2.n300.cox, beta = true.beta[[2]], simdata = n.sim)

#cure fraction = 20%, n = 1000
true.gamma <- c(1.75, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.sim11.data.n1000 <- phcshcf.gen.data(no.of.sims=n.sim, cureMin = 0.5, cureMax = 4.939527,
                                        pi.r.z.coef = true.gamma, n.obs = 1000, prob_event = 0.2,
                                        beta = true.beta)
total.noncure <- lapply(cf.sim11.data.n1000, function(a) a$noncure)
table(unlist(total.noncure)) / length(unlist(total.noncure))

table(unlist(lapply(cf.sim11.data.n1000, function(a) a$event))) /
  length(unlist(lapply(cf.sim11.data.n1000, function(a) a$event)))

#Use parallel processing
cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cf.sim11n1000.results.time <- system.time({cf.sim11n1000.results <- foreach(i = 1:n.sim, .combine = rbind, .options.snow = opts) %dopar% {
  zmat <- as.matrix(cbind(cf.sim11.data.n1000[[i]]$z1, cf.sim11.data.n1000[[i]]$z2))
  colnames(zmat) <- c("z1", "z2")
  phcsh_mpl(Surv(time, time2, event=event, "interval") ~ x1 + x2,
            risk = cf.sim11.data.n1000[[i]]$risk, data = cf.sim11.data.n1000[[i]],
            z = zmat, max.outer = 10,
            max.iter = 500, lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL,
            aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = FALSE,
            gq.points = 50)
}
})
close(pb)
stopCluster(cl)

true.gamma <- c(1.75, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.mpl.results(cf.sim11n1000.results, true.gamma, true.beta)

# cure fraction = 20% baseline hazard plot
cf.sim11.data.n1000.min <- max(sapply(cf.sim11.data.n1000, function(a)
  min(na.omit(c(a$time, a$time2)))))
cf.sim11.data.n1000.max <- min(sapply(cf.sim11.data.n1000, function(a)
  max(na.omit(c(a$time, a$time2)))))
cf1.scen1.plot.t.n1000 <- seq(cf.sim11.data.n1000.min, cf.sim11.data.n1000.max, length.out=1000)
cf1.plot.bh.r1.n1000 <- cf1.plot.bh.r2.n1000 <- list()
for(sim in 1:1000){
  cf1.plot.bh.r1.n1000[[sim]] <- pred.bh(cf.sim11n1000.results[sim,],
                                   cf1.scen1.plot.t.n1000, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  cf1.plot.bh.r2.n1000[[sim]] <- pred.bh(cf.sim11n1000.results[sim,],
                                   cf1.scen1.plot.t.n1000, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}

temp.valid = vector()
for(i in 1:1000){
  temp.valid[i] = cf.sim11n1000.results[i,]$valid
}

cf1.plot.bh.r1.all.n1000 <- lapply(cf1.plot.bh.r1.n1000, function(a) a$pred.h0r)
cf1.plot.bh.r1.h0r.mean.n1000 <- colMeans(Reduce("rbind", cf1.plot.bh.r1.all.n1000))
cf1.plot.bh.r1.all.lower.n1000 <- lapply(cf1.plot.bh.r1.n1000, function(a) a$pred.h0r.lower)
cf1.plot.bh.r1.h0r.mean.lower.n1000 <- colMeans(Reduce("rbind", cf1.plot.bh.r1.all.lower.n1000))
cf1.plot.bh.r1.all.upper.n1000 <- lapply(cf1.plot.bh.r1.n1000, function(a) a$pred.h0r.upper)
cf1.plot.bh.r1.h0r.mean.upper.n1000 <- colMeans(Reduce("rbind", cf1.plot.bh.r1.all.upper.n1000))
cf1.bh.cp.r1.n1000 <- colSums(Reduce("rbind", lapply(cf1.plot.bh.r1.n1000, function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0)))) / sum(temp.valid)


cf1.plot.bh.r2.all.n1000 <- lapply(cf1.plot.bh.r2.n1000, function(a) a$pred.h0r)
cf1.plot.bh.r2.h0r.mean.n1000 <- colMeans(Reduce("rbind", cf1.plot.bh.r2.all.n1000))
cf1.plot.bh.r2.all.lower.n1000 <- lapply(cf1.plot.bh.r2.n1000, function(a) a$pred.h0r.lower)
cf1.plot.bh.r2.h0r.mean.lower.n1000 <- colMeans(Reduce("rbind", cf1.plot.bh.r2.all.lower.n1000))
cf1.plot.bh.r2.all.upper.n1000 <- lapply(cf1.plot.bh.r2.n1000, function(a) a$pred.h0r.upper)
cf1.plot.bh.r2.h0r.mean.upper.n1000 <- colMeans(Reduce("rbind", cf1.plot.bh.r2.all.upper.n1000 ))
cf1.bh.cp.r2.n1000 <- colSums(Reduce("rbind", lapply(cf1.plot.bh.r2.n1000, function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0)))) / sum(temp.valid)

par(mfrow=c(2,2))
plot(cf1.scen1.plot.t.n1000, cf1.plot.bh.r1.h0r.mean.n1000, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(A) Risk 1", ylim=c(0,250))
lines(cf1.scen1.plot.t.n1000, cf1.plot.bh.r1.h0r.mean.lower.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.n1000, cf1.plot.bh.r1.h0r.mean.upper.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.n1000, cf1.plot.bh.r1.n1000[[1]]$true.h0r, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

plot(cf1.scen1.plot.t.n1000, cf1.plot.bh.r2.h0r.mean.n1000, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(B) Risk 2", ylim=c(0,250))
lines(cf1.scen1.plot.t.n1000, cf1.plot.bh.r2.h0r.mean.lower.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.n1000, cf1.plot.bh.r2.h0r.mean.upper.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.n1000, cf1.plot.bh.r2.n1000[[1]]$true.h0r, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')


plot(cf1.scen1.plot.t.n1000, cf1.bh.cp.r1.n1000, ylab="Coverage probability", xlab="t", col='#003366',
     type='l', ylim=c(0,1), main="(C) Risk 1")
lines(cf1.scen1.plot.t.n1000, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), lty=c(1,2),
       bty = 'n')

plot(cf1.scen1.plot.t.n1000, cf1.bh.cp.r2.n1000, ylab="Coverage probability", xlab="t", col='#003366',
     type='l', ylim=c(0,1), main="(D) Risk 2")
lines(cf1.scen1.plot.t.n1000, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), lty=c(1,2),
       bty = 'n')

#cure fraction 20% CIF plot
cf1.plot.cif.r1.n1000 <- cf1.plot.cif.r2.n1000 <- list()
for(sim in 1:1000){
  if(cf.sim11n1000.results[sim,]$valid==1){
    cf1.plot.cif.r1.n1000[[sim]] <- pred.CIF(cf.sim11n1000.results[sim,],
                                       cf1.scen1.plot.t.n1000, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
    
    cf1.plot.cif.r2.n1000[[sim]] <- pred.CIF(cf.sim11n1000.results[sim,],
                                       cf1.scen1.plot.t.n1000, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
  }
}
cf1.plot.cif.r1.all.n1000 <- lapply(cf1.plot.cif.r1.n1000, function(a) a$pred.F0r)
cf1.plot.cif.r1.F0r.mean.n1000 <- colMeans(Reduce("rbind", cf1.plot.cif.r1.all.n1000))
cf1.plot.cif.r1.all.lower.n1000 <- lapply(cf1.plot.cif.r1.n1000, function(a) a$pred.F0r.lower)
cf1.plot.cif.r1.F0r.mean.lower.n1000 <- colMeans(Reduce("rbind", cf1.plot.cif.r1.all.lower.n1000))
cf1.plot.cif.r1.all.upper.n1000 <- lapply(cf1.plot.cif.r1.n1000, function(a) a$pred.F0r.upper)
cf1.plot.cif.r1.F0r.mean.upper.n1000 <- colMeans(Reduce("rbind", cf1.plot.cif.r1.all.upper.n1000))

cf1.plot.cif.r2.all.n1000 <- lapply(cf1.plot.cif.r2.n1000, function(a) a$pred.F0r)
cf1.plot.cif.r2.F0r.mean.n1000 <- colMeans(Reduce("rbind", cf1.plot.cif.r2.all.n1000))
cf1.plot.cif.r2.all.lower.n1000 <- lapply(cf1.plot.cif.r2.n1000, function(a) a$pred.F0r.lower)
cf1.plot.cif.r2.F0r.mean.lower.n1000 <- colMeans(Reduce("rbind", cf1.plot.cif.r2.all.lower.n1000))
cf1.plot.cif.r2.all.upper.n1000 <- lapply(cf1.plot.cif.r2.n1000, function(a) a$pred.F0r.upper)
cf1.plot.cif.r2.F0r.mean.upper.n1000 <- colMeans(Reduce("rbind", cf1.plot.cif.r2.all.upper.n1000))

integrand.F01 <- function(x, lambda1, lambda2, rho){
  (lambda1*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho))
}

integrand.F02 <- function(x, lambda1, lambda2, rho){
  (lambda2*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho))
}
true.F0r.r1.n1000 <- sapply(cf1.scen1.plot.t.n1000, function(a) integrate(integrand.F01, 0, a, rho=3,
                                                              lambda1=1, lambda2 = 0.5)$value)
true.F0r.r2.n1000 <- sapply(cf1.scen1.plot.t.n1000, function(a) integrate(integrand.F02, 0, a, rho=3,
                                                              lambda1=1, lambda2=0.5)$value)
F0r.plot.mat.n1000 <- sapply(cf1.plot.cif.r1.n1000, function(a) a$pred.F0r)

lengthr1.n1000 <- lengthr2.n1000 <- list()
for(sim in 1:1000){
  lengthr1.n1000[[sim]] <- length(cf1.plot.cif.r1.n1000[[sim]]$pred.F0r)
  lengthr2.n1000[[sim]] <- length(cf1.plot.cif.r2.n1000[[sim]]$pred.F0r)
}
F0r.plot.mat.list.n1000 <- list()
count=1
for(sim in 1:1000){
  if(length(cf1.plot.cif.r1.n1000[[sim]]$pred.F0r)==1000){
    F0r.plot.mat.list.n1000[[count]] <- cf1.plot.cif.r1.n1000[[sim]]$pred.F0r
  }
  count=count+1
}
F0r.plot.mat.n1000 <- Reduce("cbind", F0r.plot.mat.list.n1000)
F0r.plot.mean.n1000 <- apply(F0r.plot.mat.n1000, 1, mean)
F0r.plot.l95.n1000 <- apply(F0r.plot.mat.n1000, 1, quantile, 0.025)
F0r.plot.u95.n1000 <- apply(F0r.plot.mat.n1000, 1, quantile, 0.975)
F0r.plot.r2.mat.list.n1000 <- list()
count=1
for(sim in 1:1000){
  if(length(cf1.plot.cif.r2.n1000[[sim]]$pred.F0r)==1000){
    F0r.plot.r2.mat.list.n1000[[count]] <- cf1.plot.cif.r2.n1000[[sim]]$pred.F0r
  }
  count=count+1
}
F0r.plot.r2.mat.n1000 <- Reduce("cbind", F0r.plot.r2.mat.list.n1000)
F0r.plot.r2.mean.n1000 <- apply(F0r.plot.r2.mat.n1000, 1, mean)
F0r.plot.r2.l95.n1000 <- apply(F0r.plot.r2.mat.n1000, 1, quantile, 0.025)
F0r.plot.r2.u95.n1000 <- apply(F0r.plot.r2.mat.n1000, 1, quantile, 0.975)

par(mfrow=c(2,1))
plot(cf1.scen1.plot.t.n1000, cf1.plot.cif.r1.F0r.mean.n1000, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), main="(A) Risk 1", xlim=c(0, 2))
lines(cf1.scen1.plot.t.n1000, cf1.plot.cif.r1.F0r.mean.lower.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.n1000, cf1.plot.cif.r1.F0r.mean.upper.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.n1000, true.F0r.r1.n1000, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

plot(cf1.scen1.plot.t.n1000, cf1.plot.cif.r2.F0r.mean.n1000, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), main="(B) Risk 2", xlim=c(0, 2))
lines(cf1.scen1.plot.t.n1000, cf1.plot.cif.r2.F0r.mean.lower.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.n1000, cf1.plot.cif.r2.F0r.mean.upper.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.n1000, true.F0r.r2.n1000, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

cf.sim1.r1.n1000.cox = cf.sim1.r2.n1000.cox = list()
for(sim in 1:n.sim){
  cf.sim1.r1.n1000.cox[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, data = cf.sim11.data.n1000[[sim]])
  cf.sim1.r2.n1000.cox[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, data = cf.sim11.data.n1000[[sim]])
}
sim.results.cox(cf.sim1.r1.n1000.cox, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(cf.sim1.r2.n1000.cox, beta = true.beta[[2]], simdata = n.sim)

#cure fraction = 40%, n = 300
true.gamma <- c(0.5, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.sim.data.n300.c40 <- phcshcf.gen.data(no.of.sims=n.sim, cureMin = 0.5, cureMax = 4.504286,
                                         pi.r.z.coef = true.gamma, n.obs = 300, prob_event = 0.2,
                                         beta = true.beta)
total.noncure <- lapply(cf.sim.data.n300.c40, function(a) a$noncure)
table(unlist(total.noncure)) / length(unlist(total.noncure))

table(unlist(lapply(cf.sim.data.n300.c40, function(a) a$event))) /
  length(unlist(lapply(cf.sim.data.n300.c40, function(a) a$event)))

allt <- unlist(lapply(cf.sim.data.n300.c40, function(a) na.omit(c(a[which(a$event != 0), ]$time,
                                                                  a[which(a$event != 0), ]$time2))))
max(allt)

#Use parallel processing
cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cf.sim.results.n300.c40.time <- system.time({cf.sim.results.n300.c40 <- foreach(i = 1:n.sim, .combine = rbind, .options.snow = opts) %dopar% {
  zmat <- as.matrix(cbind(cf.sim.data.n300.c40[[i]]$z1, cf.sim.data.n300.c40[[i]]$z2))
  colnames(zmat) <- c("z1", "z2")
  phcsh_mpl(Surv(time, time2, event=event, "interval") ~ x1 + x2,
            risk = cf.sim.data.n300.c40[[i]]$risk, data = cf.sim.data.n300.c40[[i]],
            z = zmat, max.outer = 10,
            max.iter = 500, lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL,
            aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = FALSE,
            gq.points = 50)
  
}
})
close(pb)
stopCluster(cl)

true.gamma <- c(0.5, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.mpl.results(cf.sim.results.n300.c40, true.gamma, true.beta)

#cure fraction 40%, baseline hazard plot
cf.sim.data.n300.c40.min <- max(sapply(cf.sim.data.n300.c40, function(a)
  min(na.omit(c(a$time, a$time2)))))
cf.sim.data.n300.c40.max <- min(sapply(cf.sim.data.n300.c40, function(a)
  max(na.omit(c(a$time, a$time2)))))
cf1.scen1.plot.t.c40 <- seq(cf.sim.data.n300.c40.min, cf.sim.data.n300.c40.max, length.out=1000)
cf1.plot.bh.r1.c40 <- cf1.plot.bh.r2.c40 <- list()
for(sim in 1:1000){
  cf1.plot.bh.r1.c40[[sim]] <- pred.bh(cf.sim.results.n300.c40[sim,],
                                       cf1.scen1.plot.t.c40, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  cf1.plot.bh.r2.c40[[sim]] <- pred.bh(cf.sim.results.n300.c40[sim,],
                                       cf1.scen1.plot.t.c40, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}

temp.valid = vector()
for(i in 1:1000){
  temp.valid[i] = cf.sim.results.n300.c40[i,]$valid
}

cf1.plot.bh.r1.all.c40 <- lapply(cf1.plot.bh.r1.c40, function(a) a$pred.h0r)
cf1.plot.bh.r1.h0r.mean.c40 <- colMeans(Reduce("rbind", cf1.plot.bh.r1.all.c40))
cf1.plot.bh.r1.all.lower.c40 <- lapply(cf1.plot.bh.r1.c40, function(a) a$pred.h0r.lower)
cf1.plot.bh.r1.h0r.mean.lower.c40 <- colMeans(Reduce("rbind", cf1.plot.bh.r1.all.lower.c40))
cf1.plot.bh.r1.all.upper.c40 <- lapply(cf1.plot.bh.r1.c40, function(a) a$pred.h0r.upper)
cf1.plot.bh.r1.h0r.mean.upper.c40 <- colMeans(Reduce("rbind", cf1.plot.bh.r1.all.upper.c40))
#cf1.bh.cp.r1.all <- lapply(cf1.plot.bh.r1, function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0))
cf1.bh.cp.r1.c40 <- colSums(Reduce("rbind", lapply(cf1.plot.bh.r1.c40, function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0)))) / sum(temp.valid)


cf1.plot.bh.r2.all.c40 <- lapply(cf1.plot.bh.r2.c40, function(a) a$pred.h0r)
cf1.plot.bh.r2.h0r.mean.c40 <- colMeans(Reduce("rbind", cf1.plot.bh.r2.all.c40))
cf1.plot.bh.r2.all.lower.c40 <- lapply(cf1.plot.bh.r2.c40, function(a) a$pred.h0r.lower)
cf1.plot.bh.r2.h0r.mean.lower.c40 <- colMeans(Reduce("rbind", cf1.plot.bh.r2.all.lower.c40))
cf1.plot.bh.r2.all.upper.c40 <- lapply(cf1.plot.bh.r2.c40, function(a) a$pred.h0r.upper)
cf1.plot.bh.r2.h0r.mean.upper.c40 <- colMeans(Reduce("rbind", cf1.plot.bh.r2.all.upper.c40))
#cf1.bh.cp.r2 <- lapply(cf1.plot.bh.r1, function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0))
cf1.bh.cp.r2.c40 <- colSums(Reduce("rbind", lapply(cf1.plot.bh.r2.c40, function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0)))) / sum(temp.valid)

par(mfrow=c(2,2))
plot(cf1.scen1.plot.t.c40, cf1.plot.bh.r1.h0r.mean.c40, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(A) Risk 1", ylim=c(0,250))
lines(cf1.scen1.plot.t.c40, cf1.plot.bh.r1.h0r.mean.lower.c40, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c40, cf1.plot.bh.r1.h0r.mean.upper.c40, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c40, cf1.plot.bh.r1.c40[[1]]$true.h0r, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

plot(cf1.scen1.plot.t.c40, cf1.plot.bh.r2.h0r.mean.c40, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(B) Risk 2", ylim=c(0,250))
lines(cf1.scen1.plot.t.c40, cf1.plot.bh.r2.h0r.mean.lower.c40, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c40, cf1.plot.bh.r2.h0r.mean.upper.c40, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c40, cf1.plot.bh.r2.c40[[1]]$true.h0r, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')


plot(cf1.scen1.plot.t.c40, cf1.bh.cp.r1.c40, ylab="Coverage probability", xlab="t", col='#003366',
     type='l', ylim=c(0,1), main="(C) Risk 1")
lines(cf1.scen1.plot.t.c40, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), lty=c(1,2),
       bty = 'n')

plot(cf1.scen1.plot.t.c40, cf1.bh.cp.r2.c40, ylab="Coverage probability", xlab="t", col='#003366',
     type='l', ylim=c(0,1), main="(D) Risk 2")
lines(cf1.scen1.plot.t.c40, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), lty=c(1,2),
       bty = 'n')

#cure fraction 40% CIF graph
cf1.plot.cif.r1.c40 <- cf1.plot.cif.r2.c40 <- list()
for(sim in 1:1000){
  if(cf.sim.results.n300.c40[sim,]$valid==1){
    cf1.plot.cif.r1.c40[[sim]] <- pred.CIF(cf.sim.results.n300.c40[sim,],
                                       cf1.scen1.plot.t.c40, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
    
    cf1.plot.cif.r2.c40[[sim]] <- pred.CIF(cf.sim.results.n300.c40[sim,],
                                       cf1.scen1.plot.t.c40, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
  }
}
cf1.plot.cif.r1.all.c40 <- lapply(cf1.plot.cif.r1.c40, function(a) a$pred.F0r)
cf1.plot.cif.r1.F0r.mean.c40 <- colMeans(Reduce("rbind", cf1.plot.cif.r1.all.c40))
cf1.plot.cif.r1.all.lower.c40 <- lapply(cf1.plot.cif.r1.c40, function(a) a$pred.F0r.lower)
cf1.plot.cif.r1.F0r.mean.lower.c40 <- colMeans(Reduce("rbind", cf1.plot.cif.r1.all.lower.c40))
cf1.plot.cif.r1.all.upper.c40 <- lapply(cf1.plot.cif.r1.c40, function(a) a$pred.F0r.upper)
cf1.plot.cif.r1.F0r.mean.upper.c40 <- colMeans(Reduce("rbind", cf1.plot.cif.r1.all.upper.c40))

cf1.plot.cif.r2.all.c40 <- lapply(cf1.plot.cif.r2.c40, function(a) a$pred.F0r)
cf1.plot.cif.r2.F0r.mean.c40 <- colMeans(Reduce("rbind", cf1.plot.cif.r2.all.c40))
cf1.plot.cif.r2.all.lower.c40 <- lapply(cf1.plot.cif.r2.c40, function(a) a$pred.F0r.lower)
cf1.plot.cif.r2.F0r.mean.lower.c40 <- colMeans(Reduce("rbind", cf1.plot.cif.r2.all.lower.c40))
cf1.plot.cif.r2.all.upper.c40 <- lapply(cf1.plot.cif.r2.c40, function(a) a$pred.F0r.upper)
cf1.plot.cif.r2.F0r.mean.upper.c40 <- colMeans(Reduce("rbind", cf1.plot.cif.r2.all.upper.c40))

integrand.F01 <- function(x, lambda1, lambda2, rho){
  (lambda1*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho))
}

integrand.F02 <- function(x, lambda1, lambda2, rho){
  (lambda2*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho))
}
true.F0r.r1.c40 <- sapply(cf1.scen1.plot.t.c40, function(a) integrate(integrand.F01, 0, a, rho=3,
                                                              lambda1=1, lambda2 = 0.5)$value)
true.F0r.r2.c40 <- sapply(cf1.scen1.plot.t.c40, function(a) integrate(integrand.F02, 0, a, rho=3,
                                                              lambda1=1, lambda2=0.5)$value)
F0r.plot.mat.c40 <- sapply(cf1.plot.cif.r1.c40, function(a) a$pred.F0r)

lengthr1.c40 <- lengthr2.c40 <- list()
for(sim in 1:1000){
  lengthr1.c40[[sim]] <- length(cf1.plot.cif.r1.c40[[sim]]$pred.F0r)
  lengthr2.c40[[sim]] <- length(cf1.plot.cif.r2.c40[[sim]]$pred.F0r)
}
F0r.plot.mat.list.c40 <- list()
count=1
for(sim in 1:1000){
  if(length(cf1.plot.cif.r1.c40[[sim]]$pred.F0r)==1000){
    F0r.plot.mat.list.c40[[count]] <- cf1.plot.cif.r1.c40[[sim]]$pred.F0r
  }
  count=count+1
}
F0r.plot.mat.c40 <- Reduce("cbind", F0r.plot.mat.list.c40)
F0r.plot.mean.c40 <- apply(F0r.plot.mat.c40, 1, mean)
F0r.plot.l95.c40 <- apply(F0r.plot.mat.c40, 1, quantile, 0.025)
F0r.plot.u95.c40 <- apply(F0r.plot.mat.c40, 1, quantile, 0.975)
F0r.plot.r2.mat.list.c40 <- list()
count=1
for(sim in 1:1000){
  if(length(cf1.plot.cif.r2.c40[[sim]]$pred.F0r)==1000){
    F0r.plot.r2.mat.list.c40[[count]] <- cf1.plot.cif.r2.c40[[sim]]$pred.F0r
  }
  count=count+1
}
F0r.plot.r2.mat.c40 <- Reduce("cbind", F0r.plot.r2.mat.list.c40)
F0r.plot.r2.mean.c40 <- apply(F0r.plot.r2.mat.c40, 1, mean)
F0r.plot.r2.l95.c40 <- apply(F0r.plot.r2.mat.c40, 1, quantile, 0.025)
F0r.plot.r2.u95.c40 <- apply(F0r.plot.r2.mat.c40, 1, quantile, 0.975)

par(mfrow=c(2,1))
plot(cf1.scen1.plot.t.c40, cf1.plot.cif.r1.F0r.mean.c40, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), main="(A) Risk 1", xlim=c(0, 2))
lines(cf1.scen1.plot.t.c40, cf1.plot.cif.r1.F0r.mean.lower.c40, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c40, cf1.plot.cif.r1.F0r.mean.upper.c40, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c40, true.F0r.r1.c40, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

plot(cf1.scen1.plot.t.c40, cf1.plot.cif.r2.F0r.mean.c40, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), main="(B) Risk 2", xlim=c(0, 2))
lines(cf1.scen1.plot.t.c40, cf1.plot.cif.r2.F0r.mean.lower.c40, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c40, cf1.plot.cif.r2.F0r.mean.upper.c40, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c40, true.F0r.r2.c40, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

cf.sim1.c40.n300.cox.r1 = cf.sim1.c40.n300.cox.r2 = list()
for(sim in 1:n.sim){
  cf.sim1.c40.n300.cox.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, data = cf.sim.data.n300.c40[[sim]])
  cf.sim1.c40.n300.cox.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, data = cf.sim.data.n300.c40[[sim]])
}
sim.results.cox(cf.sim1.c40.n300.cox.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(cf.sim1.c40.n300.cox.r2, beta = true.beta[[2]], simdata = n.sim)

#cure fraction 40%, n = 1000
true.gamma <- c(0.5, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.sim.data.n1000.c40 <- phcshcf.gen.data(no.of.sims=n.sim, cureMin = 0.5, cureMax = 5.507274,
                                          pi.r.z.coef = true.gamma, n.obs = 1000, prob_event = 0.2,
                                          beta = true.beta)
total.noncure <- lapply(cf.sim.data.n1000.c40, function(a) a$noncure)
table(unlist(total.noncure)) / length(unlist(total.noncure))

table(unlist(lapply(cf.sim.data.n1000.c40, function(a) a$event))) /
  length(unlist(lapply(cf.sim.data.n1000.c40, function(a) a$event)))

#Use parallel processing
cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cf.sim.results.n1000.c40.time <- system.time({cf.sim.results.n1000.c40 <- foreach(i = 1:n.sim, .combine = rbind, .options.snow = opts) %dopar% {
  zmat <- as.matrix(cbind(cf.sim.data.n1000.c40[[i]]$z1, cf.sim.data.n1000.c40[[i]]$z2))
  colnames(zmat) <- c("z1", "z2")
  phcsh_mpl(Surv(time, time2, event=event, "interval") ~ x1 + x2,
            risk = cf.sim.data.n1000.c40[[i]]$risk, data = cf.sim.data.n1000.c40[[i]],
            z = zmat, max.outer = 10,
            max.iter = 500, lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL,
            aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = FALSE,
            gq.points = 50)
}
})
close(pb)
stopCluster(cl)

true.gamma <- c(0.5, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.mpl.results(cf.sim.results.n1000.c40, true.gamma, true.beta)

#cure fraction 40%, n = 1000, baseline hazard plot
cf.sim.data.n1000.c40.min <- max(sapply(cf.sim.data.n1000.c40, function(a)
  min(na.omit(c(a$time, a$time2)))))
cf.sim.data.n1000.c40.max <- min(sapply(cf.sim.data.n1000.c40, function(a)
  max(na.omit(c(a$time, a$time2)))))
cf1.scen1.plot.t.c40.n1000 <- seq(cf.sim.data.n1000.c40.min, cf.sim.data.n1000.c40.max, length.out=1000)
cf1.plot.bh.r1.c40.n1000 <- cf1.plot.bh.r2.c40.n1000 <- list()
for(sim in 1:1000){
  cf1.plot.bh.r1.c40.n1000[[sim]] <- pred.bh(cf.sim.results.n1000.c40[sim,],
                                       cf1.scen1.plot.t.c40.n1000, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  cf1.plot.bh.r2.c40.n1000[[sim]] <- pred.bh(cf.sim.results.n1000.c40[sim,],
                                       cf1.scen1.plot.t.c40.n1000, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}

temp.valid = vector()
for(i in 1:1000){
  temp.valid[i] = cf.sim.results.n1000.c40[i,]$valid
}

cf1.plot.bh.r1.all.c40.n1000 <- lapply(cf1.plot.bh.r1.c40.n1000, function(a) a$pred.h0r)
cf1.plot.bh.r1.h0r.mean.c40.n1000 <- colMeans(Reduce("rbind", cf1.plot.bh.r1.all.c40.n1000))
cf1.plot.bh.r1.all.lower.c40.n1000 <- lapply(cf1.plot.bh.r1.c40.n1000, function(a) a$pred.h0r.lower)
cf1.plot.bh.r1.h0r.mean.lower.c40.n1000 <- colMeans(Reduce("rbind", cf1.plot.bh.r1.all.lower.c40.n1000))
cf1.plot.bh.r1.all.upper.c40.n1000 <- lapply(cf1.plot.bh.r1.c40.n1000, function(a) a$pred.h0r.upper)
cf1.plot.bh.r1.h0r.mean.upper.c40.n1000 <- colMeans(Reduce("rbind", cf1.plot.bh.r1.all.upper.c40.n1000))
cf1.bh.cp.r1.c40.n1000 <- colSums(Reduce("rbind", lapply(cf1.plot.bh.r1.c40.n1000, function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0)))) / sum(temp.valid)


cf1.plot.bh.r2.all.c40.n1000 <- lapply(cf1.plot.bh.r2.c40.n1000, function(a) a$pred.h0r)
cf1.plot.bh.r2.h0r.mean.c40.n1000 <- colMeans(Reduce("rbind", cf1.plot.bh.r2.all.c40.n1000))
cf1.plot.bh.r2.all.lower.c40.n1000 <- lapply(cf1.plot.bh.r2.c40.n1000, function(a) a$pred.h0r.lower)
cf1.plot.bh.r2.h0r.mean.lower.c40.n1000 <- colMeans(Reduce("rbind", cf1.plot.bh.r2.all.lower.c40.n1000))
cf1.plot.bh.r2.all.upper.c40.n1000 <- lapply(cf1.plot.bh.r2.c40.n1000, function(a) a$pred.h0r.upper)
cf1.plot.bh.r2.h0r.mean.upper.c40.n1000 <- colMeans(Reduce("rbind", cf1.plot.bh.r2.all.upper.c40.n1000))
cf1.bh.cp.r2.c40.n1000 <- colSums(Reduce("rbind", lapply(cf1.plot.bh.r2.c40.n1000, function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0)))) / sum(temp.valid)

par(mfrow=c(2,2))
plot(cf1.scen1.plot.t.c40.n1000, cf1.plot.bh.r1.h0r.mean.c40.n1000, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(A) Risk 1", ylim=c(0,250))
lines(cf1.scen1.plot.t.c40.n1000, cf1.plot.bh.r1.h0r.mean.lower.c40.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c40.n1000, cf1.plot.bh.r1.h0r.mean.upper.c40.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c40.n1000, cf1.plot.bh.r1.c40.n1000[[1]]$true.h0r, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

plot(cf1.scen1.plot.t.c40.n1000, cf1.plot.bh.r2.h0r.mean.c40.n1000, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(B) Risk 2", ylim=c(0,250))
lines(cf1.scen1.plot.t.c40.n1000, cf1.plot.bh.r2.h0r.mean.lower.c40.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c40.n1000, cf1.plot.bh.r2.h0r.mean.upper.c40.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c40.n1000, cf1.plot.bh.r2.c40.n1000[[1]]$true.h0r, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')


plot(cf1.scen1.plot.t.c40.n1000, cf1.bh.cp.r1.c40.n1000, ylab="Coverage probability", xlab="t", col='#003366',
     type='l', ylim=c(0,1), main="(C) Risk 1")
lines(cf1.scen1.plot.t.c40.n1000, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), lty=c(1,2),
       bty = 'n')

plot(cf1.scen1.plot.t.c40.n1000, cf1.bh.cp.r2.c40.n1000, ylab="Coverage probability", xlab="t", col='#003366',
     type='l', ylim=c(0,1), main="(D) Risk 2")
lines(cf1.scen1.plot.t.c40.n1000, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), lty=c(1,2),
       bty = 'n')

#cure fraction 40% CIF graph
cf1.plot.cif.r1.c40.n1000 <- cf1.plot.cif.r2.c40.n1000 <- list()
for(sim in 1:1000){
  if(cf.sim.results.n1000.c40[sim,]$valid==1){
    cf1.plot.cif.r1.c40.n1000[[sim]] <- pred.CIF(cf.sim.results.n1000.c40[sim,],
                                           cf1.scen1.plot.t.c40.n1000, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
    
    cf1.plot.cif.r2.c40.n1000[[sim]] <- pred.CIF(cf.sim.results.n1000.c40[sim,],
                                           cf1.scen1.plot.t.c40.n1000, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
  }
}
cf1.plot.cif.r1.all.c40.n1000 <- lapply(cf1.plot.cif.r1.c40.n1000, function(a) a$pred.F0r)
cf1.plot.cif.r1.F0r.mean.c40.n1000 <- colMeans(Reduce("rbind", cf1.plot.cif.r1.all.c40.n1000))
cf1.plot.cif.r1.all.lower.c40.n1000 <- lapply(cf1.plot.cif.r1.c40.n1000, function(a) a$pred.F0r.lower)
cf1.plot.cif.r1.F0r.mean.lower.c40.n1000 <- colMeans(Reduce("rbind", cf1.plot.cif.r1.all.lower.c40.n1000))
cf1.plot.cif.r1.all.upper.c40.n1000 <- lapply(cf1.plot.cif.r1.c40.n1000, function(a) a$pred.F0r.upper)
cf1.plot.cif.r1.F0r.mean.upper.c40.n1000 <- colMeans(Reduce("rbind", cf1.plot.cif.r1.all.upper.c40.n1000))

cf1.plot.cif.r2.all.c40.n1000 <- lapply(cf1.plot.cif.r2.c40.n1000, function(a) a$pred.F0r)
cf1.plot.cif.r2.F0r.mean.c40.n1000 <- colMeans(Reduce("rbind", cf1.plot.cif.r2.all.c40.n1000))
cf1.plot.cif.r2.all.lower.c40.n1000 <- lapply(cf1.plot.cif.r2.c40.n1000, function(a) a$pred.F0r.lower)
cf1.plot.cif.r2.F0r.mean.lower.c40.n1000 <- colMeans(Reduce("rbind", cf1.plot.cif.r2.all.lower.c40.n1000))
cf1.plot.cif.r2.all.upper.c40.n1000 <- lapply(cf1.plot.cif.r2.c40.n1000, function(a) a$pred.F0r.upper)
cf1.plot.cif.r2.F0r.mean.upper.c40.n1000 <- colMeans(Reduce("rbind", cf1.plot.cif.r2.all.upper.c40.n1000))

integrand.F01 <- function(x, lambda1, lambda2, rho){
  (lambda1*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho))
}

integrand.F02 <- function(x, lambda1, lambda2, rho){
  (lambda2*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho))
}
true.F0r.r1.c40.n1000 <- sapply(cf1.scen1.plot.t.c40.n1000, function(a) integrate(integrand.F01, 0, a, rho=3,
                                                                  lambda1=1, lambda2 = 0.5)$value)
true.F0r.r2.c40.n1000 <- sapply(cf1.scen1.plot.t.c40.n1000, function(a) integrate(integrand.F02, 0, a, rho=3,
                                                                  lambda1=1, lambda2=0.5)$value)
F0r.plot.mat.c40.n1000 <- sapply(cf1.plot.cif.r1.c40.n1000, function(a) a$pred.F0r)
#calculate length of each list
lengthr1.c40.n1000 <- lengthr2.c40.n1000 <- list()
for(sim in 1:1000){
  lengthr1.c40.n1000[[sim]] <- length(cf1.plot.cif.r1.c40.n1000[[sim]]$pred.F0r)
  lengthr2.c40.n1000[[sim]] <- length(cf1.plot.cif.r2.c40.n1000[[sim]]$pred.F0r)
}
F0r.plot.mat.list.c40.n1000 <- list()
count=1
for(sim in 1:1000){
  if(length(cf1.plot.cif.r1.c40.n1000[[sim]]$pred.F0r)==1000){
    F0r.plot.mat.list.c40.n1000[[count]] <- cf1.plot.cif.r1.c40.n1000[[sim]]$pred.F0r
  }
  count=count+1
}
F0r.plot.mat.c40.n1000 <- Reduce("cbind", F0r.plot.mat.list.c40.n1000)
F0r.plot.mean.c40.n1000 <- apply(F0r.plot.mat.c40.n1000, 1, mean)
F0r.plot.l95.c40.n1000 <- apply(F0r.plot.mat.c40.n1000, 1, quantile, 0.025)
F0r.plot.u95.c40.n1000 <- apply(F0r.plot.mat.c40.n1000, 1, quantile, 0.975)
F0r.plot.r2.mat.list.c40.n1000 <- list()
count=1
for(sim in 1:1000){
  if(length(cf1.plot.cif.r2.c40.n1000[[sim]]$pred.F0r)==1000){
    F0r.plot.r2.mat.list.c40.n1000[[count]] <- cf1.plot.cif.r2.c40.n1000[[sim]]$pred.F0r
  }
  count=count+1
}
F0r.plot.r2.mat.c40.n1000 <- Reduce("cbind", F0r.plot.r2.mat.list.c40.n1000)
F0r.plot.r2.mean.c40.n1000 <- apply(F0r.plot.r2.mat.c40.n1000, 1, mean)
F0r.plot.r2.l95.c40.n1000 <- apply(F0r.plot.r2.mat.c40.n1000, 1, quantile, 0.025)
F0r.plot.r2.u95.c40.n1000 <- apply(F0r.plot.r2.mat.c40.n1000, 1, quantile, 0.975)

par(mfrow=c(2,1))
plot(cf1.scen1.plot.t.c40.n1000, cf1.plot.cif.r1.F0r.mean.c40.n1000, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), main="(A) Risk 1", xlim=c(0, 2))
lines(cf1.scen1.plot.t.c40.n1000, cf1.plot.cif.r1.F0r.mean.lower.c40.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c40.n1000, cf1.plot.cif.r1.F0r.mean.upper.c40.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c40.n1000, true.F0r.r1.c40.n1000, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

plot(cf1.scen1.plot.t.c40.n1000, cf1.plot.cif.r2.F0r.mean.c40.n1000, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), main="(B) Risk 2", xlim=c(0, 2))
lines(cf1.scen1.plot.t.c40.n1000, cf1.plot.cif.r2.F0r.mean.lower.c40.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c40.n1000, cf1.plot.cif.r2.F0r.mean.upper.c40.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c40.n1000, true.F0r.r2.c40.n1000, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

cf.sim1.c40.n300.cox.r1 = cf.sim1.c40.n300.cox.r2 = list()
for(sim in 1:n.sim){
  cf.sim1.c40.n300.cox.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, data = cf.sim.data.n300.c40[[sim]])
  cf.sim1.c40.n300.cox.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, data = cf.sim.data.n300.c40[[sim]])
}
sim.results.cox(cf.sim1.c40.n300.cox.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(cf.sim1.c40.n300.cox.r2, beta = true.beta[[2]], simdata = n.sim)

cf.sim1.c40.n1000.cox.r1 = cf.sim1.c40.n1000.cox.r2 = list()
for(sim in 1:n.sim){
  cf.sim1.c40.n1000.cox.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, data = cf.sim.data.n1000.c40[[sim]])
  cf.sim1.c40.n1000.cox.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, data = cf.sim.data.n1000.c40[[sim]])
}
sim.results.cox(cf.sim1.c40.n1000.cox.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(cf.sim1.c40.n1000.cox.r2, beta = true.beta[[2]], simdata = n.sim)

#cure fraction 60% n = 300
true.gamma <- c(-0.5, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.sim.data.n300.c60 <- phcshcf.gen.data(no.of.sims=n.sim, cureMin = 0.5, cureMax = 4.672316,
                                         pi.r.z.coef = true.gamma, n.obs = 300, prob_event = 0.2,
                                         beta = true.beta)
total.noncure <- lapply(cf.sim.data.n300.c60, function(a) a$noncure)
table(unlist(total.noncure)) / length(unlist(total.noncure))

table(unlist(lapply(cf.sim.data.n300.c60, function(a) a$event))) /
  length(unlist(lapply(cf.sim.data.n300.c60, function(a) a$event)))

#Use parallel processing
cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cf.sim.results.n300.c60.time <- system.time({cf.sim.results.n300.c60 <- foreach(i = 1:n.sim, .combine = rbind, .options.snow = opts) %dopar% {
  zmat <- as.matrix(cbind(cf.sim.data.n300.c60[[i]]$z1, cf.sim.data.n300.c60[[i]]$z2))
  colnames(zmat) <- c("z1", "z2")
  phcsh_mpl(Surv(time, time2, event=event, "interval") ~ x1 + x2,
            risk = cf.sim.data.n300.c60[[i]]$risk, data = cf.sim.data.n300.c60[[i]],
            z = zmat, max.outer = 10,
            max.iter = 500, lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL,
            aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = FALSE,
            gq.points = 50)
}
})
close(pb)
stopCluster(cl)

true.gamma <- c(-0.5, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.mpl.results(cf.sim.results.n300.c60, true.gamma, true.beta)

#cure fraction 60%, n = 300, baseline hazard plot
cf.sim.data.n300.c60.min <- max(sapply(cf.sim.data.n300.c60, function(a)
  min(na.omit(c(a$time, a$time2)))))
cf.sim.data.n300.c60.max <- min(sapply(cf.sim.data.n300.c60, function(a)
  max(na.omit(c(a$time, a$time2)))))
cf1.scen1.plot.t.c60 <- seq(cf.sim.data.n300.c60.min, cf.sim.data.n300.c60.max, length.out=1000)
cf1.plot.bh.r1.c60 <- cf1.plot.bh.r2.c60 <- list()
for(sim in 1:1000){
  cf1.plot.bh.r1.c60[[sim]] <- pred.bh(cf.sim.results.n300.c60[sim,],
                                       cf1.scen1.plot.t.c60, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  cf1.plot.bh.r2.c60[[sim]] <- pred.bh(cf.sim.results.n300.c60[sim,],
                                       cf1.scen1.plot.t.c60, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}

temp.valid = vector()
for(i in 1:1000){
  temp.valid[i] = cf.sim.results.n300.c60[i,]$valid
}

cf1.plot.bh.r1.all.c60 <- lapply(cf1.plot.bh.r1.c60, function(a) a$pred.h0r)
cf1.plot.bh.r1.h0r.mean.c60 <- colMeans(Reduce("rbind", cf1.plot.bh.r1.all.c60))
cf1.plot.bh.r1.all.lower.c60 <- lapply(cf1.plot.bh.r1.c60, function(a) a$pred.h0r.lower)
cf1.plot.bh.r1.h0r.mean.lower.c60 <- colMeans(Reduce("rbind", cf1.plot.bh.r1.all.lower.c60))
cf1.plot.bh.r1.all.upper.c60 <- lapply(cf1.plot.bh.r1.c60, function(a) a$pred.h0r.upper)
cf1.plot.bh.r1.h0r.mean.upper.c60 <- colMeans(Reduce("rbind", cf1.plot.bh.r1.all.upper.c60))
cf1.bh.cp.r1.c60 <- colSums(Reduce("rbind", lapply(cf1.plot.bh.r1.c60, function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0)))) / sum(temp.valid)


cf1.plot.bh.r2.all.c60 <- lapply(cf1.plot.bh.r2.c60, function(a) a$pred.h0r)
cf1.plot.bh.r2.h0r.mean.c60 <- colMeans(Reduce("rbind", cf1.plot.bh.r2.all.c60))
cf1.plot.bh.r2.all.lower.c60 <- lapply(cf1.plot.bh.r2.c60, function(a) a$pred.h0r.lower)
cf1.plot.bh.r2.h0r.mean.lower.c60 <- colMeans(Reduce("rbind", cf1.plot.bh.r2.all.lower.c60))
cf1.plot.bh.r2.all.upper.c60 <- lapply(cf1.plot.bh.r2.c60, function(a) a$pred.h0r.upper)
cf1.plot.bh.r2.h0r.mean.upper.c60 <- colMeans(Reduce("rbind", cf1.plot.bh.r2.all.upper.c60))
cf1.bh.cp.r2.c60 <- colSums(Reduce("rbind", lapply(cf1.plot.bh.r2.c60, function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0)))) / sum(temp.valid)

par(mfrow=c(2,2))
plot(cf1.scen1.plot.t.c60, cf1.plot.bh.r1.h0r.mean.c60, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(A) Risk 1", ylim=c(0,250))
lines(cf1.scen1.plot.t.c60, cf1.plot.bh.r1.h0r.mean.lower.c60, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c60, cf1.plot.bh.r1.h0r.mean.upper.c60, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c60, cf1.plot.bh.r1.c60[[1]]$true.h0r, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

plot(cf1.scen1.plot.t.c60, cf1.plot.bh.r2.h0r.mean.c60, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(B) Risk 2", ylim=c(0,250))
lines(cf1.scen1.plot.t.c60, cf1.plot.bh.r2.h0r.mean.lower.c60, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c60, cf1.plot.bh.r2.h0r.mean.upper.c60, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c60, cf1.plot.bh.r2.c60[[1]]$true.h0r, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')


plot(cf1.scen1.plot.t.c60, cf1.bh.cp.r1.c60, ylab="Coverage probability", xlab="t", col='#003366',
     type='l', ylim=c(0,1), main="(C) Risk 1")
lines(cf1.scen1.plot.t.c60, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), lty=c(1,2),
       bty = 'n')

plot(cf1.scen1.plot.t.c60, cf1.bh.cp.r2.c60, ylab="Coverage probability", xlab="t", col='#003366',
     type='l', ylim=c(0,1), main="(D) Risk 2")
lines(cf1.scen1.plot.t.c60, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), lty=c(1,2),
       bty = 'n')

#cure fraction 60%, n = 300, CIF graph
cf1.plot.cif.r1.c60 <- cf1.plot.cif.r2.c60 <- list()
for(sim in 1:1000){
  if(cf.sim.results.n300.c60[sim,]$valid==1){
    cf1.plot.cif.r1.c60[[sim]] <- pred.CIF(cf.sim.results.n300.c60[sim,],
                                           cf1.scen1.plot.t.c60, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
    
    cf1.plot.cif.r2.c60[[sim]] <- pred.CIF(cf.sim.results.n300.c60[sim,],
                                           cf1.scen1.plot.t.c60, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
  }
}
cf1.plot.cif.r1.all.c60 <- lapply(cf1.plot.cif.r1.c60, function(a) a$pred.F0r)
cf1.plot.cif.r1.F0r.mean.c60 <- colMeans(Reduce("rbind", cf1.plot.cif.r1.all.c60))
cf1.plot.cif.r1.all.lower.c60 <- lapply(cf1.plot.cif.r1.c60, function(a) a$pred.F0r.lower)
cf1.plot.cif.r1.F0r.mean.lower.c60 <- colMeans(Reduce("rbind", cf1.plot.cif.r1.all.lower.c60))
cf1.plot.cif.r1.all.upper.c60 <- lapply(cf1.plot.cif.r1.c60, function(a) a$pred.F0r.upper)
cf1.plot.cif.r1.F0r.mean.upper.c60 <- colMeans(Reduce("rbind", cf1.plot.cif.r1.all.upper.c60))

cf1.plot.cif.r2.all.c60 <- lapply(cf1.plot.cif.r2.c60, function(a) a$pred.F0r)
cf1.plot.cif.r2.F0r.mean.c60 <- colMeans(Reduce("rbind", cf1.plot.cif.r2.all.c60))
cf1.plot.cif.r2.all.lower.c60 <- lapply(cf1.plot.cif.r2.c60, function(a) a$pred.F0r.lower)
cf1.plot.cif.r2.F0r.mean.lower.c60 <- colMeans(Reduce("rbind", cf1.plot.cif.r2.all.lower.c60))
cf1.plot.cif.r2.all.upper.c60 <- lapply(cf1.plot.cif.r2.c60, function(a) a$pred.F0r.upper)
cf1.plot.cif.r2.F0r.mean.upper.c60 <- colMeans(Reduce("rbind", cf1.plot.cif.r2.all.upper.c60))

integrand.F01 <- function(x, lambda1, lambda2, rho){
  (lambda1*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho))
}

integrand.F02 <- function(x, lambda1, lambda2, rho){
  (lambda2*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho))
}
true.F0r.r1.c60 <- sapply(cf1.scen1.plot.t.c60, function(a) integrate(integrand.F01, 0, a, rho=3,
                                                                      lambda1=1, lambda2 = 0.5)$value)
true.F0r.r2.c60 <- sapply(cf1.scen1.plot.t.c60, function(a) integrate(integrand.F02, 0, a, rho=3,
                                                                      lambda1=1, lambda2=0.5)$value)
F0r.plot.mat.c60 <- sapply(cf1.plot.cif.r1.c60, function(a) a$pred.F0r)
#calculate length of each list
lengthr1.c60 <- lengthr2.c60 <- list()
for(sim in 1:1000){
  lengthr1.c60[[sim]] <- length(cf1.plot.cif.r1.c60[[sim]]$pred.F0r)
  lengthr2.c60[[sim]] <- length(cf1.plot.cif.r2.c60[[sim]]$pred.F0r)
}
F0r.plot.mat.list.c60 <- list()
count=1
for(sim in 1:1000){
  if(length(cf1.plot.cif.r1.c60[[sim]]$pred.F0r)==1000){
    F0r.plot.mat.list.c60[[count]] <- cf1.plot.cif.r1.c60[[sim]]$pred.F0r
  }
  count=count+1
}
F0r.plot.mat.c60 <- Reduce("cbind", F0r.plot.mat.list.c60)
F0r.plot.mean.c60 <- apply(F0r.plot.mat.c60, 1, mean)
F0r.plot.l95.c60 <- apply(F0r.plot.mat.c60, 1, quantile, 0.025)
F0r.plot.u95.c60 <- apply(F0r.plot.mat.c60, 1, quantile, 0.975)
F0r.plot.r2.mat.list.c60 <- list()
count=1
for(sim in 1:1000){
  if(length(cf1.plot.cif.r2.c60[[sim]]$pred.F0r)==1000){
    F0r.plot.r2.mat.list.c60[[count]] <- cf1.plot.cif.r2.c60[[sim]]$pred.F0r
  }
  count=count+1
}
F0r.plot.r2.mat.c60 <- Reduce("cbind", F0r.plot.r2.mat.list.c60)
F0r.plot.r2.mean.c60 <- apply(F0r.plot.r2.mat.c60, 1, mean)
F0r.plot.r2.l95.c60 <- apply(F0r.plot.r2.mat.c60, 1, quantile, 0.025)
F0r.plot.r2.u95.c60 <- apply(F0r.plot.r2.mat.c60, 1, quantile, 0.975)

par(mfrow=c(2,1))
plot(cf1.scen1.plot.t.c60, cf1.plot.cif.r1.F0r.mean.c60, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), main="(A) Risk 1", xlim=c(0, 2))
lines(cf1.scen1.plot.t.c60, cf1.plot.cif.r1.F0r.mean.lower.c60, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c60, cf1.plot.cif.r1.F0r.mean.upper.c60, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c60, true.F0r.r1.c60, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

plot(cf1.scen1.plot.t.c60, cf1.plot.cif.r2.F0r.mean.c60, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), main="(B) Risk 2", xlim=c(0, 2))
lines(cf1.scen1.plot.t.c60, cf1.plot.cif.r2.F0r.mean.lower.c60, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c60, cf1.plot.cif.r2.F0r.mean.upper.c60, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c60, true.F0r.r2.c60, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

cf.sim1.c60.n300.cox.r1 = cf.sim1.c60.n300.cox.r2 = list()
for(sim in 1:n.sim){
  cf.sim1.c60.n300.cox.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, data = cf.sim.data.n300.c60[[sim]])
  cf.sim1.c60.n300.cox.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, data = cf.sim.data.n300.c60[[sim]])
}

sim.results.cox(cf.sim1.c60.n300.cox.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(cf.sim1.c60.n300.cox.r2, beta = true.beta[[2]], simdata = n.sim)

#cure fraction 60%, n = 1000
true.gamma <- c(-0.5, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.sim.data.n1000.c60 <- phcshcf.gen.data(no.of.sims=n.sim, cureMin = 0.5, cureMax = 4.553649,
                                          pi.r.z.coef = true.gamma, n.obs = 1000, prob_event = 0.2,
                                          beta = true.beta)
total.noncure <- lapply(cf.sim.data.n1000.c60, function(a) a$noncure)
table(unlist(total.noncure)) / length(unlist(total.noncure))

table(unlist(lapply(cf.sim.data.n1000.c60, function(a) a$event))) /
  length(unlist(lapply(cf.sim.data.n1000.c60, function(a) a$event)))


#Use parallel processing
cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cf.sim.results.n1000.c60.time <- system.time({cf.sim.results.n1000.c60 <- foreach(i = 1:n.sim, .combine = rbind, .options.snow = opts) %dopar% {
  zmat <- as.matrix(cbind(cf.sim.data.n1000.c60[[i]]$z1, cf.sim.data.n1000.c60[[i]]$z2))
  colnames(zmat) <- c("z1", "z2")
  phcsh_mpl(Surv(time, time2, event=event, "interval") ~ x1 + x2,
            risk = cf.sim.data.n1000.c60[[i]]$risk, data = cf.sim.data.n1000.c60[[i]],
            z = zmat, max.outer = 10,
            max.iter = 500, lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL,
            aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = FALSE,
            gq.points = 50)
}
})
close(pb)
stopCluster(cl)

true.gamma <- c(-0.5, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.mpl.results(cf.sim.results.n1000.c60, true.gamma, true.beta)

#figure 1
#calculate values for plots
cf.sim.data.n1000.c60.min <- max(sapply(cf.sim.data.n1000.c60, function(a)
  min(na.omit(c(a$time, a$time2)))))
cf.sim.data.n1000.c60.max <- min(sapply(cf.sim.data.n1000.c60, function(a)
  max(na.omit(c(a$time, a$time2)))))
cf1.scen1.plot.t.c60.n1000 <- seq(cf.sim.data.n1000.c60.min, cf.sim.data.n1000.c60.max, length.out=1000)
cf1.plot.bh.r1.c60.n1000 <- cf1.plot.bh.r2.c60.n1000 <- list()
for(sim in 1:1000){
  cf1.plot.bh.r1.c60.n1000[[sim]] <- pred.bh(cf.sim.results.n1000.c60[sim,],
                                       cf1.scen1.plot.t.c60.n1000, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  cf1.plot.bh.r2.c60.n1000[[sim]] <- pred.bh(cf.sim.results.n1000.c60[sim,],
                                       cf1.scen1.plot.t.c60.n1000, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}

temp.valid = vector()
for(i in 1:1000){
  temp.valid[i] = cf.sim.results.n1000.c60[i,]$valid
}

cf1.plot.bh.r1.all.c60.n1000 <- lapply(cf1.plot.bh.r1.c60.n1000, function(a) a$pred.h0r)
cf1.plot.bh.r1.h0r.mean.c60.n1000 <- colMeans(Reduce("rbind", cf1.plot.bh.r1.all.c60.n1000))
cf1.plot.bh.r1.all.lower.c60.n1000 <- lapply(cf1.plot.bh.r1.c60.n1000, function(a) a$pred.h0r.lower)
cf1.plot.bh.r1.h0r.mean.lower.c60.n1000 <- colMeans(Reduce("rbind", cf1.plot.bh.r1.all.lower.c60.n1000))
cf1.plot.bh.r1.all.upper.c60.n1000 <- lapply(cf1.plot.bh.r1.c60.n1000, function(a) a$pred.h0r.upper)
cf1.plot.bh.r1.h0r.mean.upper.c60.n1000 <- colMeans(Reduce("rbind", cf1.plot.bh.r1.all.upper.c60.n1000))
cf1.bh.cp.r1.c60.n1000 <- colSums(Reduce("rbind", lapply(cf1.plot.bh.r1.c60.n1000, function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0)))) / sum(temp.valid)


cf1.plot.bh.r2.all.c60.n1000 <- lapply(cf1.plot.bh.r2.c60.n1000, function(a) a$pred.h0r)
cf1.plot.bh.r2.h0r.mean.c60.n1000 <- colMeans(Reduce("rbind", cf1.plot.bh.r2.all.c60.n1000))
cf1.plot.bh.r2.all.lower.c60.n1000 <- lapply(cf1.plot.bh.r2.c60.n1000, function(a) a$pred.h0r.lower)
cf1.plot.bh.r2.h0r.mean.lower.c60.n1000 <- colMeans(Reduce("rbind", cf1.plot.bh.r2.all.lower.c60.n1000))
cf1.plot.bh.r2.all.upper.c60.n1000 <- lapply(cf1.plot.bh.r2.c60.n1000, function(a) a$pred.h0r.upper)
cf1.plot.bh.r2.h0r.mean.upper.c60.n1000 <- colMeans(Reduce("rbind", cf1.plot.bh.r2.all.upper.c60.n1000))
cf1.bh.cp.r2.c60.n1000 <- colSums(Reduce("rbind", lapply(cf1.plot.bh.r2.c60.n1000, function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0)))) / sum(temp.valid)

par(mfrow=c(2,2))
plot(cf1.scen1.plot.t.c60.n1000, cf1.plot.bh.r1.h0r.mean.c60.n1000, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(A) Risk 1", ylim=c(0,250))
lines(cf1.scen1.plot.t.c60.n1000, cf1.plot.bh.r1.h0r.mean.lower.c60.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c60.n1000, cf1.plot.bh.r1.h0r.mean.upper.c60.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c60.n1000, cf1.plot.bh.r1.c60.n1000[[1]]$true.h0r, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

plot(cf1.scen1.plot.t.c60.n1000, cf1.plot.bh.r2.h0r.mean.c60.n1000, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(B) Risk 2", ylim=c(0,250))
lines(cf1.scen1.plot.t.c60.n1000, cf1.plot.bh.r2.h0r.mean.lower.c60.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c60.n1000, cf1.plot.bh.r2.h0r.mean.upper.c60.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c60.n1000, cf1.plot.bh.r2.c60.n1000[[1]]$true.h0r, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')


plot(cf1.scen1.plot.t.c60.n1000, cf1.bh.cp.r1.c60.n1000, ylab="Coverage probability", xlab="t", col='#003366',
     type='l', ylim=c(0,1), main="(C) Risk 1")
lines(cf1.scen1.plot.t.c60.n1000, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), lty=c(1,2),
       bty = 'n')

plot(cf1.scen1.plot.t.c60.n1000, cf1.bh.cp.r2.c60.n1000, ylab="Coverage probability", xlab="t", col='#003366',
     type='l', ylim=c(0,1), main="(D) Risk 2")
lines(cf1.scen1.plot.t.c60.n1000, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), lty=c(1,2),
       bty = 'n')

#cure fraction 60%, n = 1000, CIF graph
cf1.plot.cif.r1.c60.n1000 <- cf1.plot.cif.r2.c60.n1000 <- list()
for(sim in 1:1000){
  if(cf.sim.results.n1000.c60[sim,]$valid==1){
    cf1.plot.cif.r1.c60.n1000[[sim]] <- pred.CIF(cf.sim.results.n1000.c60[sim,],
                                           cf1.scen1.plot.t.c60.n1000, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
    
    cf1.plot.cif.r2.c60.n1000[[sim]] <- pred.CIF(cf.sim.results.n1000.c60[sim,],
                                           cf1.scen1.plot.t.c60.n1000, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
  }
}
cf1.plot.cif.r1.all.c60.n1000 <- lapply(cf1.plot.cif.r1.c60.n1000, function(a) a$pred.F0r)
cf1.plot.cif.r1.F0r.mean.c60.n1000 <- colMeans(Reduce("rbind", cf1.plot.cif.r1.all.c60.n1000))
cf1.plot.cif.r1.all.lower.c60.n1000 <- lapply(cf1.plot.cif.r1.c60.n1000, function(a) a$pred.F0r.lower)
cf1.plot.cif.r1.F0r.mean.lower.c60.n1000 <- colMeans(Reduce("rbind", cf1.plot.cif.r1.all.lower.c60.n1000))
cf1.plot.cif.r1.all.upper.c60.n1000 <- lapply(cf1.plot.cif.r1.c60.n1000, function(a) a$pred.F0r.upper)
cf1.plot.cif.r1.F0r.mean.upper.c60.n1000 <- colMeans(Reduce("rbind", cf1.plot.cif.r1.all.upper.c60.n1000))

cf1.plot.cif.r2.all.c60.n1000 <- lapply(cf1.plot.cif.r2.c60.n1000, function(a) a$pred.F0r)
cf1.plot.cif.r2.F0r.mean.c60.n1000 <- colMeans(Reduce("rbind", cf1.plot.cif.r2.all.c60.n1000))
cf1.plot.cif.r2.all.lower.c60.n1000 <- lapply(cf1.plot.cif.r2.c60.n1000, function(a) a$pred.F0r.lower)
cf1.plot.cif.r2.F0r.mean.lower.c60.n1000 <- colMeans(Reduce("rbind", cf1.plot.cif.r2.all.lower.c60.n1000))
cf1.plot.cif.r2.all.upper.c60.n1000 <- lapply(cf1.plot.cif.r2.c60.n1000, function(a) a$pred.F0r.upper)
cf1.plot.cif.r2.F0r.mean.upper.c60.n1000 <- colMeans(Reduce("rbind", cf1.plot.cif.r2.all.upper.c60.n1000))

integrand.F01 <- function(x, lambda1, lambda2, rho){
  (lambda1*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho))
}

integrand.F02 <- function(x, lambda1, lambda2, rho){
  (lambda2*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho))
}
true.F0r.r1.c60.n1000 <- sapply(cf1.scen1.plot.t.c60.n1000, function(a) integrate(integrand.F01, 0, a, rho=3,
                                                                      lambda1=1, lambda2 = 0.5)$value)
true.F0r.r2.c60.n1000 <- sapply(cf1.scen1.plot.t.c60.n1000, function(a) integrate(integrand.F02, 0, a, rho=3,
                                                                      lambda1=1, lambda2=0.5)$value)
F0r.plot.mat.c60.n1000 <- sapply(cf1.plot.cif.r1.c60.n1000, function(a) a$pred.F0r)
#calculate length of each list
lengthr1.c60.n1000 <- lengthr2.c60.n1000 <- list()
for(sim in 1:1000){
  lengthr1.c60.n1000[[sim]] <- length(cf1.plot.cif.r1.c60.n1000[[sim]]$pred.F0r)
  lengthr2.c60.n1000[[sim]] <- length(cf1.plot.cif.r2.c60.n1000[[sim]]$pred.F0r)
}
F0r.plot.mat.list.c60.n1000 <- list()
count=1
for(sim in 1:1000){
  if(length(cf1.plot.cif.r1.c60.n1000[[sim]]$pred.F0r)==1000){
    F0r.plot.mat.list.c60.n1000[[count]] <- cf1.plot.cif.r1.c60.n1000[[sim]]$pred.F0r
  }
  count=count+1
}
F0r.plot.mat.c60.n1000 <- Reduce("cbind", F0r.plot.mat.list.c60.n1000)
F0r.plot.mean.c60.n1000 <- apply(F0r.plot.mat.c60.n1000, 1, mean)
F0r.plot.l95.c60.n1000 <- apply(F0r.plot.mat.c60.n1000, 1, quantile, 0.025)
F0r.plot.u95.c60.n1000 <- apply(F0r.plot.mat.c60.n1000, 1, quantile, 0.975)
F0r.plot.r2.mat.list.c60.n1000 <- list()
count=1
for(sim in 1:1000){
  if(length(cf1.plot.cif.r2.c60.n1000[[sim]]$pred.F0r)==1000){
    F0r.plot.r2.mat.list.c60.n1000[[count]] <- cf1.plot.cif.r2.c60.n1000[[sim]]$pred.F0r
  }
  count=count+1
}
F0r.plot.r2.mat.c60.n1000 <- Reduce("cbind", F0r.plot.r2.mat.list.c60.n1000)
F0r.plot.r2.mean.c60.n1000 <- apply(F0r.plot.r2.mat.c60.n1000, 1, mean)
F0r.plot.r2.l95.c60.n1000 <- apply(F0r.plot.r2.mat.c60.n1000, 1, quantile, 0.025)
F0r.plot.r2.u95.c60.n1000 <- apply(F0r.plot.r2.mat.c60.n1000, 1, quantile, 0.975)

par(mfrow=c(2,1))
plot(cf1.scen1.plot.t.c60.n1000, cf1.plot.cif.r1.F0r.mean.c60.n1000, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), main="(A) Risk 1", xlim=c(0, 2))
lines(cf1.scen1.plot.t.c60.n1000, cf1.plot.cif.r1.F0r.mean.lower.c60.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c60.n1000, cf1.plot.cif.r1.F0r.mean.upper.c60.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c60.n1000, true.F0r.r1.c60.n1000, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

plot(cf1.scen1.plot.t.c60.n1000, cf1.plot.cif.r2.F0r.mean.c60.n1000, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), main="(B) Risk 2", xlim=c(0, 2))
lines(cf1.scen1.plot.t.c60.n1000, cf1.plot.cif.r2.F0r.mean.lower.c60.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c60.n1000, cf1.plot.cif.r2.F0r.mean.upper.c60.n1000, col='#003366',
      lty='dashed')
lines(cf1.scen1.plot.t.c60.n1000, true.F0r.r2.c60.n1000, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

cf.sim1.c60.n300.cox.r1 = cf.sim1.c60.n300.cox.r2 = list()
for(sim in 1:n.sim){
  cf.sim1.c60.n300.cox.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, data = cf.sim.data.n300.c60[[sim]])
  cf.sim1.c60.n300.cox.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, data = cf.sim.data.n300.c60[[sim]])
}

cf.sim1.c60.n1000.cox.r1 = cf.sim1.c60.n1000.cox.r2 = list()
for(sim in 1:n.sim){
  cf.sim1.c60.n1000.cox.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, data = cf.sim.data.n1000.c60[[sim]])
  cf.sim1.c60.n1000.cox.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, data = cf.sim.data.n1000.c60[[sim]])
}
sim.results.cox(cf.sim1.c60.n1000.cox.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(cf.sim1.c60.n1000.cox.r2, beta = true.beta[[2]], simdata = n.sim)