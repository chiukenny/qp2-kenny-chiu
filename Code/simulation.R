# R script for our attempt to replicate Hussey & Hughes' (2007)
# simulation study.

# Set the parameters below before running.


# Simulation parameters
# ---------------------------------------------------------
save.files = F  # Save estimates and SEs: {T,F}
run.LMM    = T  # Run LMM: {T,F}
run.GLMM   = F  # Run GLMM: {T,F}
run.GEE    = F  # Run GEE: {T,F}

RR           = 0.5      # Risk ratio: {0.5,0.6,0.7,1}
equal.size   = T        # Equal cluster size: {T,F}
sims         = 1000     # Number of simulations: >=1
link_f       = "logit"  # Link function: {"logit","identity"}
jackknife.se = F        # Use jackknife variance: {T,F}




# Initialization
# ---------------------------------------------------------
library(tidyverse)
library(nlme)
library(MASS)
library(gee)
library(DirichletReg)

set.seed(1)


# Set parameters and constants
# ---------------------------------------------------------
N = 100  # Number of units per cluster in equal size case
I = 24   # Number of clusters
TT = 5   # Number of time points
N.total = N*I

mu = 0.05
tau2 = 0.000225
theta = mu*(RR-1)

# Index constants for convenience
ILMM = 1
IGLMM = 2
IGEE = 3


# Run simulation
# ---------------------------------------------------------

# Initialize result data frames
thetas = as.data.frame(matrix(0,sims,3)) # Estimates
colnames(thetas) = c("LMM","GLMM","GEE")

SEs = as.data.frame(matrix(0,sims,4)) # Standard errors
colnames(SEs) = c("LMM","GLMM","GEE.naive","GEE.robust")

fails = as.data.frame(matrix(0,sims,3)) # Number of failures
colnames(fails) = c("LMM","GLMM","GEE")

# Run simulation
sink("NUL")
for (n in 1:sims)
{
  message(paste("Starting simulation", n))
  
  # Sample crossover times and cluster effects
  times = sample(rep(1:4, I/(TT-1))) + 1
  alphas = rnorm(I, mean=0, sd=sqrt(tau2))
  
  # Generate individual-level data
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
  for (i in 1:I)
  {
    p = rep(pmax(0, mu+alphas[i]+theta*(1:TT>=times[i])), each=Ni[i])
    dat = rbind(dat, data.frame(cluster=i,
                                time=rep(1:TT, each=Ni[i]),
                                crossover=times[i],
                                case=rbinom(TT*Ni[i], 1, p)))
  }
  
  # Order data by clusters for GLMM and GEE
  dat = arrange(dat, cluster, time)
  dat$cluster = factor(dat$cluster)
  
  # Fit LMM
  if (run.LMM)
  {
    # Aggregate individual-level data
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
        
        if (jackknife.se)
        {
          # Compute jackknife estimate of variance
          jk.thetas = rep(0, I)
          for (i in 1:I)
          {
            jk.lmm = lme(preval~ept+time, random=~1|cluster,
                         data=dat.cluster[which(dat.cluster$cluster!=i),])
            jk.thetas[i] = summary(jk.lmm)$tTable[2,1]
          }
          jk.pseudo = (N.total*ept.lmm[1] - (N.total-Ni)*jk.thetas) / Ni
          jk.mean = sum(Ni*jk.pseudo) / N.total
          SEs[n,ILMM] = sqrt(sum(Ni^2*(jk.pseudo-jk.mean)^2)/(N.total^2))
        } else {
          # Extract standard error from model output
          SEs[n,ILMM] = ept.lmm[2]
        }
        0 # Return 0 if no error
      },
      error = function(err)
      {
        # Flag error if model fitting failed
        message(paste("LMM failed:", err))
        return(1)
      }
    )
  }
  
  # Fit GLMM
  if (run.GLMM)
  {
    fails[n,IGLMM] = tryCatch(
      {
        fit.glmm = glmmPQL(case~I(time >= crossover)+time, random=~1|cluster,
                           family=binomial(link=link_f), data=dat, verbose=F)
        ept.glmm = summary(fit.glmm)$tTable[2,]
        thetas[n,IGLMM] = ept.glmm[1]
        
        if (jackknife.se)
        {
          # Compute jackknife estimate of variance
          jk.thetas = rep(0, I)
          for (i in 1:I)
          {
            jk.glmm = glmmPQL(case~I(time >= crossover)+time, random=~1|cluster,
                              family=binomial(link=link_f),
                              data=dat[which(dat$cluster!=i),], verbose=F)
            jk.thetas[i] = summary(jk.glmm)$tTable[2,1]
          }
          jk.pseudo = (N.total*ept.glmm[1] - (N.total-Ni)*jk.thetas) / Ni
          jk.mean = sum(Ni*jk.pseudo) / N.total
          SEs[n,IGLMM] = sqrt(sum(Ni^2*(jk.pseudo-jk.mean)^2)/(N.total^2))
        } else {
          # Extract standard error from model output
          SEs[n,IGLMM] = ept.glmm[2]
        }
        0 # Return 0 if not error
      },
      error = function(err)
      {
        # Flag error if model fitting failed
        message(paste("GLMM failed:", err))
        return(1)
      }
    )
  }
  
  # Fit GEE
  if (run.GEE)
  {
    fails[n,IGEE] = tryCatch(
      {
        # gee() prints unwanted messages so suppress them
        ept.gee = summary(suppressMessages(
          gee(case~I(time >= crossover)+factor(time), cluster,
              data=dat, corstr="exchangeable", family=binomial(link=link_f))
          ))$coef[2,]
        thetas[n,IGEE] = ept.gee[1]
        
        # Jackknife estimate of variance for GEE not implemented
        SEs[n,c(IGEE,IGEE+1)] = c(ept.gee[2], ept.gee[4])
        
        0 # Return 0 if no error
      }, error = function(err)
      {
        # Flag error if model fitting failed
        message(paste("GEE failed:", err))
        return(1)
      }
    )
  }
}
sink()


# Save outputs
# ---------------------------------------------------------
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


# Print results and diagnostics
# ---------------------------------------------------------

# Failures
print("Number of failures:")
colSums(fails)

# Compute power
print("Estimated powers:")
powers = c(mean(abs(thetas[!fails[,1],1])/SEs[!fails[,1],1] > qnorm(0.975)),
           mean(abs(thetas[!fails[,2],2])/SEs[!fails[,2],2] > qnorm(0.975)),
           mean(abs(thetas[!fails[,3],3])/SEs[!fails[,3],3] > qnorm(0.975)),
           mean(abs(thetas[!fails[,3],3])/SEs[!fails[,3],4] > qnorm(0.975)))
names(powers) = c("LMM","GLMM","GEE.naive","GEE.robust")
powers

# Plot generated clusters
# dat.cluster %>%
#   arrange(crossover, cluster) %>%
#   group_by(time) %>%
#   mutate(sorted.cluster=row_number()) %>%
#   ggplot(aes(time, sorted.cluster, fill=preval)) +
#   geom_tile() +
#   scale_x_discrete(expand=c(0,0)) +
#   scale_y_continuous(expand=c(0,0)) +
#   scale_fill_gradient(limits=c(0,0.13)) +
#   labs(x="Time", y="Sorted cluster", fill="Prevalence") +
#   theme_bw()