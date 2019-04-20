# this file contains code for linear mixed effects models of the orthodontic example
# Note: we use the R pacakge nlme here. An alternative is the lme4 package.
#       The lmer function in lme4 is more suitable for modeling multiple non-nested random effects
# By Gen Li, 3/19/2018
#setwd("C:\\Users\\gl2521\\Dropbox\\Courses\\Teaching\\Columbia - Bios Method II\\2018 Spring\\Code")
library(nlme)
library(ggplot2)
library (lattice)
Orthodont<- read.table ("./orthodontic.dat",header=TRUE)
head(Orthodont)

Orth.new <- groupedData (distance ~ age | child,data = as.data.frame (Orthodont)) # not necessary for lme
head(Orth.new)
OrthFem <- subset(Orth.new, male==0)
OrthFem[1:5,]
plot(OrthFem) # default: ordered by max resp of each child
ggplot(Orth.new, aes(age, distance, group=child)) + geom_line()


######################################
# fit a random intercept model
LMM1 <- lme (distance ~ age, random = ~1 | child,  data = OrthFem, method='REML') 
# group by "child"
summary (LMM1) # pay attention to: random effects, fixed effects, 
# correlation is about fisher information
VarCorr(LMM1) # covariance estimates for random effects and variance for residuals
LMM1$sigma # std for residuals
vcov(LMM1) # covariance for fixed effects estimates (inv fisher info)
#
fixed.effects(LMM1) # fixed effects coeff 
random.effects(LMM1) # ordered random effects, BLUP (in this case, just b_i)
fitted(LMM1) # fixed+random for each subj in each visit
OrthFem$distance-fitted(LMM1) # residuals
LMM1$residuals
#
plot(OrthFem$age[1:4],OrthFem$distance[1:4],type='b',ylim=c(19,26),col=26) # original data for subj1
points(OrthFem$age[1:4],fixed.effects(LMM1)[1]+fixed.effects(LMM1)[2]*OrthFem$age[1:4], pch=4,col=26)  # fix effect
lines(OrthFem$age[1:4],fitted(LMM1)[1:4], lty=5,type='b',pch=4,col=26)  # fixed effect add radom effect
random.effects(LMM1) # check random effects b_i for the first subj: -1.229
#
lines(OrthFem$age[5:8],OrthFem$distance[5:8],type='b',col=20) # original data for subj2
points(OrthFem$age[5:8],fixed.effects(LMM1)[1]+fixed.effects(LMM1)[2]*OrthFem$age[5:8], pch=4,col=20)
lines(OrthFem$age[5:8],fitted(LMM1)[5:8], lty=5,type='b',pch=4,col=20)
random.effects(LMM1) #check random effects b_i for the first subj: 0.340


# check equivalence to marginal model with compound symmetry correlation structure
summary(gls(distance~age, OrthFem, correlation=corCompSymm(form = ~ 1 |child), method="REML"))
# check rho==sigma_b^2/(sigma_b^2+sigma^2)


###########################################
# compare models (likelihood ratio test)
LMM.1 <- lme (distance ~ age, random = ~ 1 | child,  data = OrthFem, method='ML') # do NOT use REML for likelihood ratio among nested models 
LMM.2 <- lme (distance ~ 1, random = ~ 1 | child,  data = OrthFem, method='ML')
anova(LMM.2,LMM.1) 





#############################################
# fit a random intercept and slope model
LMM2 <- lme (distance ~ age, random = ~ 1+ age | child, data = OrthFem)
summary (LMM2) 
#
plot(OrthFem$age[1:4],OrthFem$distance[1:4],type='b',ylim=c(19,26),col=26) # original data for child1
points(OrthFem$age[1:4],fixed.effects(LMM2)[1]+fixed.effects(LMM2)[2]*OrthFem$age[1:4], pch=4,col=26)
lines(OrthFem$age[1:4],fitted(LMM2)[1:4], lty=5,type='b',pch=4,col=26)
random.effects(LMM2) # check BLUP (Q: why the diff between fixed est and mixed est not equal to random intercept?)
lines(OrthFem$age[5:8],OrthFem$distance[5:8],type='b',col=20) 
points(OrthFem$age[5:8],fixed.effects(LMM2)[1]+fixed.effects(LMM2)[2]*OrthFem$age[5:8], pch=4,col=20)
lines(OrthFem$age[5:8],fitted(LMM2)[5:8], lty=5,type='b',pch=4,col=20)







# ##############################
# # # check diff fit
# resp=c(3.1,2.1,1.0, 6.2,5.1,4.0,  8.3,7.1,6.2, 10.5,9.3,8.4, 13.2,12.1,11.2)
# pred=c(1,2,3, 2,4,6 , 1,4,5, 0,2,4, 4,6,7)
# #pred=c(1,2,5, 1,2,5, 1,2,5, 1,2,5, 1,2,5)
# group=rep(c(1,2,3,4,5),c(3,3,3,3,3))
# test=groupedData(resp~pred|group)
# plot(test)
# summary(lm(resp~pred)) # fixed effect est
# summary(lme (resp ~ pred, random = ~ 1 | group, data =test,method='REML')) # random intercept
# summary(lme (resp ~ pred, random = ~ 1+pred | group, data =test,method='REML')) # random intercept and slope
# summary(gls(resp~pred ,test, correlation=corAR1(form = ~ 1 |group),method="REML")) # marginal with autoregressive cov
# summary(gls(resp~pred ,test, correlation=corCompSymm(form = ~ 1 |group),method="REML")) # marginal with compound symmetry cov (conceptually, this is equiv to random intercept model)
# summary(gls(resp~pred ,test, correlation=corSymm(form = ~ 1 |group),method="REML")) # marginal with no constraint cov

