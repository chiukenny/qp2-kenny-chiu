library(tidyverse)
library(nlme)
library(MASS)
library(gee)
# library(geeM)
library(DirichletReg)

set.seed(1)

# Run settings
# ---------------------------
save.files = T
run.LMM = F
run.GLMM = F
run.GEE = T

RR = 0.7 # Risk ratio
equal.size = F
sims = 100

link_f = "logit" #"identity"

# Parameters
# ---------------------------
N = 100
I = 24
TT = 5 # 4 randomization steps
N.total = 2400

mu = 0.05 # Prevalence
tau2 = 0.000225
theta = mu*(RR-1)

#sigma2 = mu*(1-mu)
#sigma2 = (mu+theta)*(1-mu-theta)

#rejects = matrix(0, sims, 3)
thetas = as.data.frame(matrix(0,sims,3))
colnames(thetas) = c("LMM","GLMM","GEE")
SEs = as.data.frame(matrix(0,sims,4))
colnames(SEs) = c("LMM","GLMM","GEE.naive","GEE.robust")

fails = as.data.frame(matrix(0,sims,3))
colnames(fails) = c("LMM","GLMM","GEE")
#powers = matrix(0, sims, 3)

ILMM = 1
IGLMM = 2
IGEE = 3

# Simulation
# ---------------------------
sink("NUL")
for (n in 1:sims)
{
  message(paste("Starting simulation", n))
  # Sample crossover times and simulate cluster effects
  times = sample(rep(1:4, I/(TT-1))) + 1
  alphas = rnorm(I, mean=0, sd=sqrt(tau2))
  
  dat = data.frame(cluster=integer(),
                   time=integer(),
                   crossover=integer(),
                   case=numeric())
  
  if (equal.size)
  {
    Ni = rep(N, I)
  } else {
    Ni = rmultinom(1, N.total-I, rdirichlet(1,rep(1,I))) + 1
  }
  
  # for (j in 1:TT)
  # {
  #   p = rep(pmax(0, mu + alphas + theta*(j >= times)), each=N)
  #   #p = rep(mu + alphas + theta*(j >= times), each=N)
  #   dat = rbind(dat, data.frame(cluster=rep(1:I, each=N),
  #                               time=j,
  #                               crossover=rep(times, each=N),
  #                               case=rbinom(I*N, 1, p)))
  #                               #case=rnorm(I*N, mean=p, sd=sqrt(sigma2))))
  # }
  for (i in 1:I)
  {
    p = rep(pmax(0, mu+alphas[i]+theta*(1:TT>=times[i])), each=Ni[i])
    #p = rep(mu + alphas[i] + theta*(1:TT >= times[i]), each=Ni[i])
    dat = rbind(dat, data.frame(cluster=i,
                                time=rep(1:TT, each=Ni[i]),
                                crossover=times[i],
                                case=rbinom(TT*Ni[i], 1, p)))
    #case=rnorm(TT*Ni[i], mean=p, sd=sqrt(sigma2))))
  }
  dat = arrange(dat, cluster, time)
  dat$cluster = factor(dat$cluster)
  
  # Analysis
  # -------------------------
  
  # LMM
  if (run.LMM)
  {
    dat.cluster = dat %>%
      group_by(cluster, time, crossover) %>%
      summarize(preval=mean(case))
    dat.cluster$ept = factor((dat.cluster$time >= dat.cluster$crossover)*1)
    dat.cluster$cluster = factor(dat.cluster$cluster)
    dat.cluster$time = factor(dat.cluster$time)
    
    fails[n,ILMM] = tryCatch(
      {
        fit.lmm = lme(preval~ept+time, random=~1|cluster, data=dat.cluster)
        ept.lmm = summary(fit.lmm)$tTable[2,]
        thetas[n,ILMM] = ept.lmm[1]
        SEs[n,ILMM] = ept.lmm[2]
        
        # DELETE THIS
        # thetas[n,IGLMM] = ept.lmm[1]
        # 
        # jk.thetas = rep(0, I)
        # for (i in 1:I)
        # {
        #   M = sum(Ni)
        #   jk.lmm = lme(preval~ept+time, random=~1|cluster, data=dat.cluster[which(dat.cluster$cluster!=i),])
        #   jk.thetas[i] = summary(jk.lmm)$tTable[2,1]
        # }
        # jk.pseudo = (M*ept.lmm[1] - (M-Ni)*jk.thetas) / Ni
        # jk.mean = sum(Ni*jk.pseudo) / M
        # SEs[n,IGLMM] = sqrt(sum(Ni^2*(jk.pseudo-jk.mean)^2)/(M^2))
        
        0 # Return 0
      },
      error = function(err)
      {
        message(paste("LMM failed:", err))
        return(1)
      }
    )
  }
  #rejects[n] = abs(ept[1]/ept[2]) > qnorm(0.975)
  
  # dat.cluster %>%
  #   arrange(crossover, cluster) %>%
  #   group_by(time) %>%
  #   mutate(sorted.cluster=row_number()) %>%
  #   ggplot(aes(time, sorted.cluster, fill=preval)) +
  #   geom_tile()
  
  #rejects[n] = (summary(fit)$tTable[2,"p-value"] < 0.05)*1
  #powers[n] = pnorm(abs(theta/ept[2]) - qnorm(0.975))
  #thetas[n] = ept[1]
  #rejects[n] = (ept[5] < 0.05)*1
  #rejects[n] = (pnorm(ept[1]/ept[2]) < 0.025)*1
  #rejects[n] = (pnorm(ept[1]/sqrt(sigma2/N)) < 0.025)*1
  #rejects[n] = pnorm(-theta/sqrt(sigma2/N)-qnorm(0.975))
  #rejects[n] = (pnorm(ept[1]/sqrt(wls.var)-qnorm(0.975)) < 0.025)*1
  
  # GLMM
  if (run.GLMM)
  {
    fails[n,IGLMM] = tryCatch(
      {
        fit.glmm = glmmPQL(case~I(time >= crossover)+time, random=~1|cluster, family=binomial(link=link_f), data=dat, verbose=F)
        ept.glmm = summary(fit.glmm)$tTable[2,]
        thetas[n,IGLMM] = ept.glmm[1]
        SEs[n,IGLMM] = ept.glmm[2]
        
        # DELETE THIS
        # SEs[n,ILMM] = ept.glmm[2]
        # 
        # # jackknife
        # jk.thetas = rep(0, I)
        # for (i in 1:I)
        # {
        #   M = sum(Ni)
        #   jk.glmm = glmmPQL(case~I(time >= crossover)+time, random=~1|cluster, family=binomial(link=link_f), data=dat[which(dat$cluster!=i),], verbose=F)
        #   jk.thetas[i] = summary(jk.glmm)$tTable[2,1]
        # }
        # jk.pseudo = (M*ept.glmm[1] - (M-Ni)*jk.thetas) / Ni
        # jk.mean = sum(Ni*jk.pseudo) / M
        # SEs[n,IGLMM] = sqrt(sum(Ni^2*(jk.pseudo-jk.mean)^2)/(M^2))
        
        0 # Return 0
      },
      error = function(err)
      {
        message(paste("GLMM failed:", err))
        return(1)
      }
    )
  }
  # rejects[n] = abs(ept[1]/ept[2]) > qnorm(0.975)
  
  #summary(fit.glmm)
  #thetas[n] = coefs[2]
  #rejects[n] = (summary(fit.glmm)$tTable[2,"p-value"] < 0.05)*1
  #rejects[n] = (pnorm(coefs[2]/sqrt(wls.var)-qnorm(0.975)) < 0.025)*1
  
  # GEE
  if (run.GEE)
  {
    fails[n,IGEE] = tryCatch(
      {
        ept.gee <- summary(suppressMessages(gee(case~I(time >= crossover)+factor(time), cluster, data=dat, corstr="exchangeable", family=binomial(link=link_f))))$coef[2,]
        thetas[n,IGEE] = ept.gee[1]
        SEs[n,c(IGEE,IGEE+1)] = c(ept.gee[2], ept.gee[4])
        0 # Return 0
      }, error = function(err)
      {
        message(paste("GEE failed:", err))
        return(1)
      }
    )
    # ept.gee = summary(geem(case~I(time >= crossover)+factor(time), cluster, data=dat, corstr="exchangeable", family="binomial", tol=0.1))$coef[2,]
    # #ept = summary(fit.gee)$coef[2,]
    #rejects[n] = abs(ept[1]/ept[4]) > qnorm(0.975)
    
    #fit.gee = gee(case~I(time >= crossover)+factor(time), factor(cluster), data=dat, corstr="exchangeable")
  }
}
sink()

# Save files
if (save.files)
{
  if (equal.size)
  {
    thetas_f = "thetas_eq.txt"
    SEs_f = "SEs_eq.txt"
  } else {
    thetas_f = "thetas_uneq.txt"
    SEs_f = "SEs_uneq.txt"
  }
  write.table(thetas, file=thetas_f, quote=F, row.names=F)
  write.table(thetas, file=SEs_f, quote=F, row.names=F)
}

# Failures
colSums(fails)

# Compute power
powers = c(mean(abs(thetas[!fails[,1],1])/SEs[!fails[,1],1] > qnorm(0.975)),
           mean(abs(thetas[!fails[,2],2])/SEs[!fails[,2],2] > qnorm(0.975)),
           mean(abs(thetas[!fails[,3],3])/SEs[!fails[,3],3] > qnorm(0.975)),
           mean(abs(thetas[!fails[,3],3])/SEs[!fails[,3],4] > qnorm(0.975)))
names(powers) = c("LMM","GLMM","GEE.naive","GEE.robust")
powers
#powers = colSums((abs(cbind(thetas,thetas[,3]))/SEs) > qnorm(0.975)) / sims

# table(rejects)

# jk.tt = rep(0, sims)
# for (i in 1:sims)
# {
#   jk.tt[i] = mean(thetas[-i])
# }
# jk.var = sum((thetas-mean(jk.tt))^2) / (sims*(sims-1))
# pnorm(-theta/sqrt(jk.var) - qnorm(0.975))

#pnorm(-theta/sd(thetas))

#pnorm(-theta/sd(thetas)-qnorm(0.975))
#pnorm(-theta/sqrt(wls.var)-qnorm(0.975))





# OTHER

# WLS power analysis
# 0.05 --> 0.032 (-36%) should have power of ~0.8
# U = I*TT/2
# W = sum(((1:4)*6)^2)
# V = sum(((1:4)^2)*6)
# wls.var = I*(sigma2/N)*((sigma2/N)+TT*tau2) / ( (I*U-W)*(sigma2/N) + (U^2 + I*TT*U - TT*W - I*V)*tau2 )

#wls.var = I*TT*(sigma2/N)*((sigma2/N)+TT*tau2) / ((I*TT*U-U^2)*(sigma2/N)+I*TT*(U*TT-V)*tau2)

# wls.var.f = function(mu)
# {
#   sigma2 = mu*(1-mu)
#   wls.var = I*(sigma2/N)*((sigma2/N)+TT*tau2) / ( (I*U-W)*(sigma2/N) + (U^2 + I*TT*U - TT*W - I*V)*tau2 )
#   return(pnorm(-theta/sqrt(wls.var)-qnorm(0.975)))
# }

# Power calculation
# eff.sizes = 20:50
# powers = rep(0, length(eff.sizes))
# for (i in 1:length(eff.sizes))
# {
#   tt = -mu*eff.sizes[i]/100
#   powers[i] = pnorm(-tt/sqrt(wls.var)-qnorm(0.975))
# }
#plot(eff.sizes, powers, type="l", ylim=c(0,1))