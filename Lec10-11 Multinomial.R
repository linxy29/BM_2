# Analyze Car Preference Data (Nominal Logistic regression)
# By Gen Li 9/2/2016

car <- data.frame(res.unim=c(26, 9, 5, 40, 17, 8), res.im=c(12, 21, 14, 17, 15, 15), res.veim=c(7, 15, 41, 8, 12, 18), sex=c(rep("F", 3), rep("M",3)),age=rep(c("18-23", "24-40", ">40"), 2)) 
car$age=factor(car$age, levels=c("18-23", "24-40", ">40"))
car

#####################################
#### nominal logistic regression ####
#####################################
library(nnet)
car.mult <- multinom(cbind(res.unim, res.im, res.veim)~sex+age, data=car)
summary(car.mult)

# # an equivalent way to fit the model
# freq=c(car$res.unim,car$res.im,car$res.veim)
# res=c(rep(c("unim","im","veim"),c(6,6,6)))
# res=factor(res,levels=c("unim","im","veim"))
# car.ord<-data.frame(res=res,sex=rep(car$sex,3),age=rep(car$age,3),freq=freq)
# car.mult1 <- multinom(res~sex+age,data=car.ord,weights=freq)
# summary(car.mult1)


# goodness of fit
pihat=predict(car.mult,type='probs') 
m=rowSums(car[,1:3])
res.pearson=(car[,1:3]-pihat*m)/sqrt(pihat*m) # pearson residuals 
# NOTE: do not use residual function directly for multinomial models
# residuals(car.mult, type='pearson') gives you car[,1:3]/m-pihat, which is WRONG!
# do not use deviance from the function either
# It's easy to derive the deviance 
G.stat=sum(res.pearson^2) # Generalized Pearson Chisq Stat
G.stat
pval=1-pchisq(G.stat,df=(6-4)*(3-1)) # fit is good
# deviance
D.stat=sum(2*car[,1:3]*log(car[,1:3]/(m*pihat)))
D.stat





# alternative model is beta_age3=2*beta_age2
car$age.new=factor(car$age, levels=c("18-23", "24-40", ">40"))
levels(car$age.new)=c(0,1,2)
car$age.new=as.numeric(car$age.new)-1
car
car.alt <- multinom(cbind(res.unim, res.im, res.veim)~sex+age.new, data=car)
summary(car.alt)
# gof
pihat1=predict(car.alt,type='probs') 
res.pearson1=(car[,1:3]-pihat1*m)/sqrt(pihat1*m) # pearson residuals
G.stat1=sum(res.pearson1^2) # Generalized Pearson Chisq Stat
G.stat1 # df=(3-1)*(6-3)=6
D.stat1=sum(2*car[,1:3]*log(car[,1:3]/(m*pihat1)))
D.stat1

# deviance analysis
D.stat1-D.stat # chisq(2)
pval=1-pchisq(D.stat1-D.stat ,df=2)
pval # not rejected, go with the smaller model
 







#####################################
#### ordinal logistic regression ####
#####################################
library(MASS)
freq=c(car$res.unim,car$res.im,car$res.veim)
res=c(rep(c("unim","im","veim"),c(6,6,6)))
res=factor(res,levels=c("unim","im","veim"),ordered=T)
car.ord<-data.frame(res=res,sex=rep(car$sex,3),age=rep(car$age,3),freq=freq)
car.ord

# fit proportional odds model
car.polr=polr(res~sex+age,data=car.ord,weights=freq)
summary(car.polr) # pay attention to sign, read help doc of polr (-eta)

# model prediction
predict(car.polr,car,type='p')
predict(car.polr,data.frame(sex='M',age='18-23'),type='p') # predict prob for each category (another option is 'class')

# residuals
pihat=predict(car.polr,car,type='p')
m=rowSums(cbind(car$res.im,car$res.unim,car$res.veim))
res.pearson=(car[,1:3]-pihat*m)/sqrt(pihat*m)
G=sum(res.pearson^2)
G
numsamp=(3-1)*6 # degree of freedom for grouped data
numparam=2+3 # total num of param
pval=1-pchisq(G ,df=numsamp-numparam)
pval # fits well

# *** exercise: use what you have learned from this class to test beta_age3=2beta_age2 in this model


# combine categories and use logistic regression (coeff other than intercept should be similar to the prop odds model)
res1=car.ord$res
res2=car.ord$res
levels(res1)=c('unim','im','im')
levels(res2)=c('unim','unim','veim')
car.ord1=data.frame(res=res,res1=res1,res2=res2,sex=rep(car$sex,3),age=rep(car$age,3),freq=freq)
car.ord1

car.logit1=glm(res1~sex+age, data=car.ord1, weights=freq,family=binomial(link='logit'))
summary(car.logit1) 
car.logit2=glm(res2~sex+age, data=car.ord1, weights=freq,family=binomial(link='logit'))
summary(car.logit2)

