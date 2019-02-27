# Compare binomial model (ignoring dispersion), constant over-dispersion
# By Gen Li 8/27/2016

#setwd("C:\\Users\\gl2521\\Dropbox\\Courses\\Teaching\\Columbia - Bios Method II\\2018 Spring\\Code")
library(aod) # for wald.test
#library(dispmod) # for glm.binomial.disp

# input data
x=c(rep(0,16),rep(1,16))
y=c(13,12,9,9,8,8,12,11,9,9,8,11,4,5,7,7,12,11,10,9,10,9,9,8,8,4,7,4,5,3,3,0) # survive=1
m=c(13,12,9,9,8,8,13,12,10,10,9,13,5,7,10,10,12,11,10,9,11,10,10,9,9,5,9,7,10,6,10,7)
data=data.frame(x,y,m)


# fit binomial (logistic) without dispersion
none.disp=glm(cbind(y,m-y)~x, family=binomial(link='logit'))
summary(none.disp)
# goodness of fit
pval=1-pchisq(none.disp$deviance,32-2)
pval # bad fit

# wald test of beta_2=0 (equiv to z-test)
wald.test(b = coef(none.disp), Sigma = vcov(none.disp), Terms=2)
# deviance analysis of beta_2=0 (nested model)
test.stat=none.disp$null.deviance-none.disp$deviance
pval=1-pchisq(test.stat,df=1)
pval # rej the smaller model; go with the larger model


# calc dispersion param
G.stat=sum(residuals(none.disp,type='pearson')^2) # pearson chisq 
G.stat
phi=G.stat/(32-2)
phi
tilde.phi=none.disp$deviance/none.disp$df.residual
tilde.phi # similar to the one estimated from pearson chisq 



#######################
# test over-dispersion (half normal plot)
res=residuals(none.disp,type='pearson')
plot(qnorm((32+1:32+0.5)/(2*32+1.125)),sort(abs(res)),xlab='Expected Half-Normal Order Stats',ylab='Ordered Abs Pearson Residuals')
abline(a=0,b=1)
abline(a=0,b=sqrt(phi),lty=2)

########################



# fit model with constant over-dispersion
summary(none.disp,dispersion=phi)

