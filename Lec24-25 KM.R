# this file contains code for parametric and nonparametric estimation of survival function
# for acute leukemia clinical trial
# By Gen Li, 1/4/2018

library(flexsurv)
library(survival)
library(MASS)
data(gehan)
Surv(gehan$time,gehan$cens,type='right') # 0: censored, 1: observed


# parametric survival function (for 6-mp group)
########################################################
param1 <- flexsurvreg(Surv(time, cens) ~ 1, data = subset(gehan, treat=="6-MP"),
                      dist = "exp")  # S(t)=e^{-rate*t}
param2 <- flexsurvreg(Surv(time, cens) ~ 1, data = subset(gehan, treat=="6-MP"),
                     dist = "weibull") # S(t)=e^{-(t/scale)^shape}
param2 # Weibull parameter estimation and CI
summary(param2) # survival function estimation and CI
plot(param2, xlab="Months", ylab="Survival Probability", main="6MP (KM and Parametric Est with 95% CI)", cex.lab=1.5, cex.main=1.5)



# KM survival function
###############################################
KM=survfit(Surv(time,cens)~1, data = subset(gehan, treat=="6-MP"), conf.type='log')
plot(KM, conf.int = FALSE, mark.time = TRUE,xlab="Months", ylab="Survival Probability", main="6MP K-M curve", cex.lab=1.5, cex.main=1.5)
plot(KM)
plot(KM,fun='cumhaz') # cumulative hazard fun
# estimate cumulative hazard rates
cbind(KM$time,-log(KM$surv),  cumsum(KM$n.event/KM$n.risk)) # time, KM est, Nelson-Aalen Estimator
#
# obtain survival rate at given time, with CI
summary(KM,time=c(5,10,12.5, 15)) # note: n.event is the cumulative num of events since last listed time 
summary(KM, censored = TRUE)# (if not specify time, then n.event is the # event at each time point)
# median survival time, with CI
print(KM)



# Log Rank test
#############################################
survdiff(Surv(time,cens)~treat, data=gehan) # log rank test
plot(survfit(Surv(time,cens)~treat, data = gehan)) 
library(survminer)
ggsurvplot( survfit(Surv(time, cens) ~ treat, data = gehan), conf.int=TRUE)
