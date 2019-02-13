# Analyze Beetle Mortality Data
# By Gen Li 8/19/2016


# load data
Dose=c(49.05689, 52.99074, 56.91150, 60.84151, 64.75898, 68.69103, 72.61060, 76.54203)
log10Dose=log10(Dose)
num=c(59, 60, 62, 56, 63, 59, 62, 60)
killed=c(6, 13, 18, 28, 52, 53, 61, 60)
data=data.frame(log10Dose,num,killed)

# data visualization
plot(data$log10Dose,data$killed/data$num,xlab='log10(Dose)',ylab='Proportion Killed',cex=1.5,pch=19,cex.lab=1.6,cex.axis = 1.5)

# data preparation
x=data$log10Dose
y=data$killed  # number of cases
m=data$num     # number of total cases
resp=cbind(y,m-y)   # response should include two factors, number of cases and remaining case
#resp=cbind(m-y,y) # surv=1, death=0   # for predicting number of survival


### fit logistic model
glm_logit=glm(resp~x, family=binomial(link='logit'))
#glm_probit=glm(resp~x, family=binomial(link='probit'))
#glm_clog=glm(resp~x, family=binomial(link='cloglog')) # asymmetric
summary(glm_logit) # wald test of coefficients

### logistic model with ungrouped data
new.x=c(rep(x,y),rep(x,m-y))
new.resp=c(rep(1,sum(y)),rep(0,sum(m-y)))
newdata=data.frame(new.resp,new.x)
head(newdata)
glm_logit1=glm(new.resp~new.x, family=binomial(link='logit'),data=newdata)
summary(glm_logit1)







######################################
# Goodness of fit/residuals (check ?residuals.glm to learn more about residuals)
beta0=glm_logit$coefficients[1]
beta1=glm_logit$coefficients[2]
beta0+x*beta1
pihat=fitted(glm_logit) # \hat{pi}
G.res=(y-m*pihat)/sqrt(m*pihat*(1-pihat))  # Pearson Chisq residual
residuals(glm_logit, type = "pearson")
sum(residuals(glm_logit,type='pearson')^2) # pearson chisq 
dev=sum(residuals(glm_logit,type='deviance')^2);dev # deviance (or obtain from summary(glm_logit)) 
# compare with chisq(8-2)
pval=1-pchisq(dev,6);pval # fit is ok, fails to reject
# 
# check Hosmer-Lemeshow for ungrouped data
library(ResourceSelection)
hl <- hoslem.test(glm_logit1$y, fitted(glm_logit1), g=10)  # fitted: returns \hat{pi}
hl # again, fit is ok, fails to reject


 
# Confidence interval 
# CI for beta
vcov(glm_logit) # covariance of beta MLE (fisher information inverse)
beta=glm_logit$coefficients[2]
se=sqrt(vcov(glm_logit)[2,2]) # (same as in summary(glm_logit))
beta+c(qnorm(0.025),0,-qnorm(0.025))*se

# CI for odds ratio: exp(beta)
exp(beta+c(qnorm(0.025),0,-qnorm(0.025))*se)
 
# CI for x\beta
out=predict(glm_logit, se.fit=TRUE);out # predict: returns x\hat{beta}
#?predict.glm # check other options; can also return \hat{pi}
out=predict(glm_logit, data.frame(x=c(1.7)), se.fit=TRUE)
CIval=out$fit[1]+c(qnorm(0.025),0,-qnorm(0.025))*out$se.fit[1] # 95% CI  for eta, based on delta method

# CI for pi
predict(glm_logit, data.frame(x=1.7), se.fit=TRUE,type='response') # predict: \hat{pi}
exp(CIval)/(1+exp(CIval)) # 95% CI for pi
#
# # Note: do NOT use the following to get CI for pi (b/c g^-1(eta) deviates from Gaussian)
# out=predict(glm_logit, data.frame(x=c(1.7)), se.fit=TRUE, type='response')
# wrong.CIval=out$fit+c(qnorm(0.025),0,-qnorm(0.025))*out$se.fit
# wrong.CIval



# LD50 est and CI
beta0=glm_logit$coefficients[1]
beta1=glm_logit$coefficients[2]
betacov=vcov(glm_logit) # inverse fisher information
x0fit=-beta0/beta1
10^x0fit # point estimate of LD50
varx0=betacov[1,1]/(beta1^2)+betacov[2,2]*(beta0^2)/(beta1^4)-2*betacov[1,2]*beta0/(beta1^3)
c(x0fit,sqrt(varx0)) # point est and se
10^(x0fit+c(qnorm(0.025),-qnorm(0.025))*sqrt(varx0)) # 95% CI for LD50




