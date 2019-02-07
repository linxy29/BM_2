# Fit logistic models


# EX1: Show/No-show
data1 = read.table('Lagtime_1.csv',header=TRUE,sep=',')
table(data1[,1]) # 521 show, 125 no-show
dim(data1) # 645*2
names(data1)
levels(data1$Internal.Status) # if resp is a factor, the first level is treated as 0 and second level is treated as 1
data1$Internal.Status1 <- relevel(data1$Internal.Status, "NOSHOW") # switch level order; if we need to test arrival, then set NOSHOW as referrence

fit=glm(Internal.Status~Appointment.Lag,family=binomial(link='logit'),data=data1)   # family: specify link function to be logit
summary(fit) # residual deviance is deviance of fitted model,but the data is ungroup data, comparing residual deviance with chi-square cannot be used
exp(fit$coefficients)[2] # odds ratio of no-show with one day increase in lag 
vcov(fit)
sqrt(vcov(fit)[1,1])   # should be equal to se of intercept
# GoF (ungrouped data) check Hosmer-Lemeshow for ungrouped data
library(ResourceSelection)
hoslem.test(fit$y, fitted(fit), g=10)  # fitted: returns \hat{pi}
# fails to reject, fit is ok
# 95% CI for beta
CI1=fit$coefficients + kronecker(t(c(0,qnorm(0.025),-qnorm(0.025))),t(t(sqrt(diag(vcov(fit))))))
#out=cbind(exp(CI1)[-1,,drop=FALSE],coef(summary(fit))[-1,4,drop=FALSE])
#colnames(out)=c('OR','95% CI','95% CI','p-value')
#rownames(out)=c('Lag')
#out


# EX2: Pub
data2 = read.table('MedEd Stats.csv',header=TRUE,sep=',')
sum(data2$timeoff) # 26 timeoff, 162 no timeoff
# fit model
# (#pr, #total-#pr) ~ timeoff
resp=cbind(data2$Urology.Publication,data2$Total.Publications-data2$Urology.Publication)
pred=data2$timeoff
fit=glm(resp~pred,family=binomial(link='logit'))
summary(fit)
exp(fit$coefficients)[2] # odds ratio of urology pub with timeoff vs no-timeoff
# GoF
#sum(residuals(fit,type='pearson')^2) # pearson chisq 
#dev=sum(residuals(fit,type='deviance')^2);dev # deviance (or obtain from summary(glm_logit)) 
## compare with chisq(188-2)
#pval=1-pchisq(dev,186);pval # fit is not good, later will see why (over dispersion)


# equivalently
uropub=xtabs(data2$Urology.Publication~data2$timeoff)
allpub=xtabs(data2$Total.Publications~data2$timeoff)
resp1=cbind(uropub,allpub-uropub)
fit=glm(resp1~c(0,1),family=binomial(link='logit'))
summary(fit) 
exp(fit$coefficients)[2] 



# 95% CI
#CI1=fit$coefficients + kronecker(t(c(0,qnorm(0.025),-qnorm(0.025))),t(t(sqrt(diag(vcov(fit))))))
#out=cbind(exp(CI1)[-1,,drop=FALSE],coef(summary(fit))[-1,4,drop=FALSE])
#colnames(out)=c('OR','95% CI','95% CI','p-value')
#rownames(out)=c('timeoff')
#out # p=0.008, OR=1.31(1.07,1.59) (odds of PR pub in timeoff vs no timeoff)

