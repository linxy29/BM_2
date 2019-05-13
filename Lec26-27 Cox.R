# this file contains code for Cox PH
# Standard PH and Stratification (time-dependent omitted)
# By Gen Li, 1/8/2018

library(survival)
library(KMsurv) # contains many interesting data sets
?coxph # note: ties, strata, tt



# standard cox ph
data(btrial)
fit=coxph(Surv(time,death)~im, data=btrial)
summary(fit)
# num of events is num of uncensored obs
# CI of HR exp(beta) is based on CI of beta
# interpretation of beta
# predict.coxph is with respect to mean covariate
#
# get survival function (including baseline hazard function) of a patient
S0=survfit(fit,newdata=data.frame(im=1)) # ?survfit.coxph,  estimate baseline h(t) to get the survival function
S1=survfit(fit,newdata=data.frame(im=2)) 
plot(S0$time,S0$surv,xlab='time',ylim=c(0,1),ylab='Estimated Survival Rate',type='s',lty=1) # evaluated at all times (censored or uncensored)
lines(S1$time,S1$surv,type='s',col='blue',lty=2)
text(x=100,y=0.7,'group 1')
text(x=100,y=0.15,'group 2') # with higher hazard/poorer survival






#########################################################
# cox ph with stratification
data(larynx)
# treat stage as factor (stage=1 as reference)
fit0=coxph(Surv(time,delta)~factor(stage)+age,data=larynx,ties='breslow') # with tied uncensored data; factor(stage)
# estimate survival function for a 60 year old with stage 2
Sfit.ph=survfit(fit0,newdata=data.frame(age=60,stage=c(1,2,3,4))) # survival rate evaluated at every uncensored time point in the original data
summary(Sfit.ph)
plot(Sfit.ph,col=c('black','blue','magenta', 'red'),xlab='time',ylab='survival rate',main='Surv Rate of Patient at Age 60, from PH Model')
legend('topright',c('Stage I','Stage II','Stage III','Stage IV'),lty=1,col=c('black','blue','magenta', 'red'))
#
#
# fit cox with stratification
fit1=coxph(Surv(time,delta)~age+strata(stage),data=larynx,ties='breslow')   #strate(stage)
summary(fit1)
Sfit.strata=survfit(fit1,newdata=data.frame(age=60,stage=c(1,2,3,4))) # write out model
plot(Sfit.strata,col=c('black','blue','magenta', 'red'),xlab='time',ylab='survival rate',main='Surv Rate of Patient at Age 60, from PH Model with Strata')
legend('topright',c('Stage I','Stage II','Stage III','Stage IV'),lty=1,col=c('black','blue','magenta', 'red'))


