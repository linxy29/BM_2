# Compare retrospective model and prospective model
# By Gen Li 8/23/2016

# case control study data (i.e., retrospective data)
#################################
#              case      control
# exposed      n11=1707  n10=2455     n1
# unexposed    n01=1293  n00=3321     n0
#              m1=3000    m0=5776        
#################################
n11=1707
n10=2455
n01=1293
n00=3321
m0=n10+n00 # for retro model
m1=n11+n01
n1=n11+n10 # for prosp model
n0=n01+n00



# prospective model
x=c(1,0)
y=cbind(c(n11,n01),c(n10,n00))
# fit logit model
logit.prosp=glm(y~x, family=binomial(link='logit'))
summary(logit.prosp)

# retrospective model
x=c(1,0)
y=cbind(c(n11,n10),c(n01,n00))
# fit logit model
logit.retro=glm(y~x, family=binomial(link='logit'))
summary(logit.retro)

# This is retro data, so we cannot use it to calc RR (even with prosp model)
