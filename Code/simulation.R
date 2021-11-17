library(tidyverse)
library(nlme)

set.seed(1)

# Parameters
N = 100
I = 24
TT = 5 # 4 randomization steps
N.total = 2400

mu = 0.05 # Prevalence
tau2 = 0.000225
sigma2 = mu*(1-mu)

RR = 0.8 # Risk ratio
theta = mu*(RR-1)

sims = 100


# WLS power analysis
# 0.05 --> 0.032 (-36%) should have power of ~0.8
U = I*TT/2
W = sum(((1:4)*6)^2)
V = sum(((1:4)^2)*6)
wc.var = I*(sigma2/N)*((sigma2/N)+TT*tau2) / ( (I*U-W)*(sigma2/N) + (U^2 + I*TT*U - TT*W - I*V)*tau2 )
pnorm(-theta/sqrt(wc.var)-qnorm(0.975))

# Simulate data (equal case)
# --------------------------

rejects = rep(0, sims)

for (n in 1:sims)
{

# Sample crossover times and simulate cluster effects
times = sample(rep(1:4, I/(TT-1))) + 1
alphas = rnorm(I, mean=0, sd=sqrt(tau2))

dat = data.frame(cluster=integer(),
                 time=integer(),
                 crossover=integer(),
                 case=integer())
for (j in 1:TT)
{
  p = rep(pmax(0, mu + alphas + theta*(j >= times)), each=N)
  dat = rbind(dat, data.frame(cluster=rep(1:I, each=N),
                              time=j,
                              crossover=rep(times, each=N),
                              case=rbinom(I*N, 1, p)))
}

dat.cluster = dat %>%
  group_by(cluster, time, crossover) %>%
  summarize(preval=mean(case))
# dat.cluster %>%
#   arrange(crossover, cluster) %>%
#   group_by(time) %>%
#   mutate(sorted.cluster=row_number()) %>%
#   ggplot(aes(time, sorted.cluster, fill=preval)) +
#   geom_tile()

# Analysis
dat.cluster$ept = factor((dat.cluster$time >= dat.cluster$crossover)*1)
dat.cluster$cluster = factor(dat.cluster$cluster)
dat.cluster$time = factor(dat.cluster$time)

# LME
fit = lme(preval~1+ept, random=~1|cluster, data=dat.cluster)
#summary(fit)
#fit = lme(preval~1+time+ept, random=~time|cluster, data=dat.cluster)
#summary(fit)
#plot(fit, resid(.,type="p")~fitted(.)|cluster)

#rejects[n] = (summary(fit)$tTable[2,"p-value"] < 0.05)*1
ept = summary(fit)$tTable[2,]
#rejects[n] = (pnorm(ept[1]/ept[2]) < 0.025)*1
# rejects[n] = (pnorm(ept[1]/sqrt(sigma2/N)) < 0.025)*1
rejects[n] = pnorm(-theta/sqrt(sigma2/N)-qnorm(0.975))
}
table(rejects)