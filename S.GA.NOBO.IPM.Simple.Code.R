###############################################################################
## Code for running an integrated population model estimating juvenile survival
##    of northern bobwhites (Colinus virginianus) in southern Georgia, USA, 
##    1998 - 2022. Juveniles survival is inferred from breeding productivity and
##    post-breeding population count and harvest age/sex ratio data.
##  Code for: William B. Lewis, Chloe R. Nater, Justin A. Rectenwald, D. Clay 
##    Sisson, and James A. Martin. 2024. Use of breeding and post-breeding data
##    to estimate juvenile survival with integrated population models. International
##    Statistical Ecology Conference.
##############################################################################

require(parallel)

load("S.GA.NOBO.IPM.data.simple.gzip")

str(NOBO.IPM.data)

# See metadata on GitHub page for description of gzip contents. Female-specific model.
# Juveniles enter population in breeding months (June - September). Juveniles transition
#   to subadults in October and to adults in April.

# Putting model and initial values into function for running in parallel
NOBO_MCMC_code <- function(seed, mod.data, mod.constants, monitors){
  
  require(nimble)
  require(truncnorm)
  
  NOBO.model.code <- nimbleCode({
    
    # Priors --------------------------------------------------------------------
    
    # Initial population sizes for June population 1st column is surviving adults
    #   2nd column is surviving subadults
    for(r in 1:2){
      S.n.f[r] ~ dunif(0,1000)
    }
    S.nb[1,1] <- round(S.n.f[1])
    S.nb[1,2] <- round(S.n.f[2])
    
    # Informative prior for covey size and availability, both vary by year
    for (q in 1:nyears){
      F.Tot.R[q] ~ dunif(0,1)
      covey.size[q] ~ T(dnorm(csize.mean[q], sd=csize.sd[q]), csize.min, csize.max)
      availability.logit[q] ~ T(dnorm(avail.mean[q], sd=avail.sd[q]), avail.min, avail.max)
      availability[q] <- ilogit(availability.logit[q])
    }
    
    # Informative prior for detection of coveys. Doesn't vary by year
    pdet.logit ~ T(dnorm(pdet_mean, sd=pdet_sd), pdet_min, Inf)
    pdet <- ilogit(pdet.logit)
    
    
    for (d in 1:nyears){ # Looping over years
      
      # Juvenile survival
      phi.j.daily.logit[d] ~ dnorm(phi.j.prior.mean, sd=phi.j.prior.sd)
      phi.j.daily[d] <- ilogit(phi.j.daily.logit[d])
      phi.j[d] <- pow(phi.j.daily[d],30.5)
      J.surv[d,1:n.b.months] <- c(pow(phi.j[d],4),
                                  pow(phi.j[d],3),
                                  pow(phi.j[d],2),
                                  phi.j[d])
      
      # Breeding survival
      phi.b.logit[d] ~ dnorm(phi.b.hyper.logit, sd=phi.b.hyper.logit.sd) # Yearly values come from distribution around hyperparameter
      phi.b[d] <- ilogit(phi.b.logit[d]) # Not doing density-dependence
      
      # Non-breeding survival
      for(f in 1:2){ # Looping over ages (1=adult, 2=subadult)
        phi.nb.logit[d,f] ~ dnorm(phi.nb.hyper.logit[f], sd=phi.nb.hyper.logit.sd[f]) # Yearly values come from distribution around hyperparameter
        phi.nb[d,f] <- ilogit(phi.nb.logit[d,f]) # Not doing density-dependence
      }
      
      # Productivity
      for(g in 1:n.b.months){
        productivity.log[d,g] ~ dnorm(prod.year.hyper.log[g], sd=prod.year.hyper.log.sd[g]) # Yearly values come from distribution around hyperparameter
        productivity[d,g] <- exp(productivity.log[d,g]) # Not doing density dependence
      }
    }  
    
    phi.b.hyper.logit <- logit(Survival.b.hyper)
    Survival.b.hyper ~ dunif(0, 1)
    phi.b.hyper.logit.sd ~ dunif(0, 100)
    for(h in 1:2){ # Looping over ages (1=adults, 2=subadults)
      phi.nb.hyper.logit[h] <- logit(Survival.nb.hyper[h])
      Survival.nb.hyper[h] ~ dunif(0, 1)
      phi.nb.hyper.logit.sd[h] ~ dunif(0, 100)
    }
    for(u in 1:n.b.months){ # Looping over months of chick production (June-September)
      prod.year.hyper.log[u] <- log(Prod.hyper[u])
      Prod.hyper[u] ~ dgamma(1,1)
      prod.year.hyper.log.sd[u] ~ dunif(0, 100)
    }
    
    
    
    # State Process ------------------------------------------------------------
    
    
    ## First Year --------------------------------------------------------------
    N.b[1,1] <- S.nb[1,1] + S.nb[1,2] # Adults alive at start of June
    
    # Adults surviving to subsequent months of the breeding season
    for(m in 2:n.b.months){
      N.b[1,m] ~ dbin(phi.b[1], N.b[1,m-1])
    }
    # Adults surviving from last breeding month (September) to start of November (one month of breeding, one of non-breeding)
    N.Nov[1,1] ~ dbin(phi.b[1] * phi.nb[1,1], N.b[1,n.b.months])
    
    # Productivity and juvenile survival model. Assuming 50/50 age ratio of chicks
    for (k in 1:n.b.months){
      Chicks[1,k] ~ dpois(productivity[1,k] * N.b[1,k] * 0.5) # Female chicks produced
      Recruit[1,k] ~ dbin(J.surv[1,k] * phi.nb[1,2], Chicks[1,k]) # Subadult recruits over juvenile period (hatch month - start of October) and one month of non-breeding to start of November
    }
    
    # Ratios
    N.Nov[1,2] <- sum(Recruit[1,1:n.b.months])
    N.Nov.Tot.F[1] <- N.Nov[1,1] + N.Nov[1,2]
    JF.F.R[1] <- N.Nov[1,2] / N.Nov.Tot.F[1] # Ratio of juvenile females to total females
    N.Nov.Tot[1] ~ dpois(N.Nov.Tot.F[1] / F.Tot.R[1]) # Total population size is total females divided by ratio of females in population
    D.Nov.Tot[1] <- N.Nov.Tot[1] / Area
    
    
    
    ## Subsequent Years --------------------------------------------------------
    
    for(t in 2:nyears){
      
      # Survivors from November of previous year to June (5 months of non-breeding, 2 months of breeding)
      for(d in 1:2){ # Looping over ages
        S.nb[t,d] ~ dbin(pow(phi.nb[t-1,d],5) * pow(phi.b[t],2), N.Nov[t-1,d])
      }
      N.b[t,1] <- S.nb[t,1] + S.nb[t,2] # Adults alive at start of June
      
      # Adults surviving to subsequent months of the breeding season
      for(m in 2:n.b.months){
        N.b[t,m] ~ dbin(phi.b[t], N.b[t,m-1])
      }
      # Adults surviving from last breeding month (September) to start of November (one month of breeding, one of non-breeding)
      N.Nov[t,1] ~ dbin(phi.b[t] * phi.nb[t,1], N.b[t,n.b.months])
      
      # Productivity and juvenile survival model. Assuming 50/50 age ratio of chicks
      for (k in 1:n.b.months){
        Chicks[t,k] ~ dpois(productivity[t,k] * N.b[t,k] * 0.5) # Female chicks produced
        Recruit[t,k] ~ dbin(J.surv[t,k] * phi.nb[t,2], Chicks[t,k]) # Subadult recruits over juvenile period (hatch month - start of October) and one month of non-breeding to start of November
      }
      
      # Ratio of juveniles to adults
      N.Nov[t,2] <- sum(Recruit[t,1:n.b.months])
      N.Nov.Tot.F[t] <- N.Nov[t,1] + N.Nov[t,2]
      JF.F.R[t] <- N.Nov[t,2] / N.Nov.Tot.F[t] # Ratio of juvenile females to total females
      N.Nov.Tot[t] ~ dpois(N.Nov.Tot.F[t] / F.Tot.R[t]) # Total population size is total females divided by ratio of females in population
      D.Nov.Tot[t] <- N.Nov.Tot[t] / Area
    }
    
    
    
    # Likelihoods --------------------------------------------------------------
    
    for(r in 1:nyears){
      
      ## Harvest Ratios --------------------------------------------------------
      Ratio[r,1:3] <- c(F.Tot.R[r]*(1-JF.F.R[r]),F.Tot.R[r]*JF.F.R[r],(1-F.Tot.R[r]))  # Ratios of adult females, juvenile females, and males is function of JF.F.R and F.Tot.R
      Harv.sexage[r,1:3] ~ dmulti(Ratio[r,1:3], Harv.N[r])
      
      
      ## Productivity ----------------------------------------------------------
      for (u in 1:n.b.months){
        chicks.nest.mean[r,u] ~ T(dnorm(chicks.nest.act[r,u], sd=chicks.nest.sd[r,u]), 0, Inf) # Incorporating uncertainty in chick production due to missing data
        chicks.nest.act[r,u] ~ dpois(productivity[r,u] * N.tracked[r,u])
      }
      
      ## Fall covey counts -----------------------------------------------------
      count[r] ~ dbinom(availability[r] * pdet, coveys[r]) # Number of coveys detected is a function of number of coveys, availability, and detection probability
      coveys[r] ~ dpois(Avail[r] / covey.size[r]) # Number of coveys is a function of number available divided by average covey size
      Avail[r] ~ dbinom(effort[r],N.Nov.Tot[r]) # The number of bobwhite available to be sampled is conditional on total population size and degree of area surveyed
    }
    
    ## Telemetry -------------------------------------------------------------
    # Placing survival estimates in matrix, converting from monthly to biweekly survival
    CH.phi[1:nyears,1] <- pow(phi.nb[1:nyears,1],0.5)
    CH.phi[1:nyears,2] <- pow(phi.nb[1:nyears,2],0.5)
    CH.phi[1:nyears,3] <- pow(phi.b[1:nyears],0.5)
    
    for (b in 1:nCH){
      for (c in 1:trackperiods[b]){
        CH.state[b,trackdates[b,c]] ~ dbern(CH.phi[CH.year[trackdates[b,c]-1], CH.age[b,trackdates[b,c]-1]])
      }
    } 
    
    
  })
  
  
  # Setting initial values -----------------------------------------------------
  # Seed ensures that initial values are acceptable: 
  #   Initial Total November population size 1 - 5000
  #   Initial June population size 1 - 3000
  #   Initial harvest ratios never reach 0
  #   Initial covey number never less than observed covey count
  set.seed(374)
  N.Nov.Tot.init <- Avail.init <- coveys.init <- JF.F.R.init <- N.Nov.Tot.F.init <- rep(NA, times=mod.constants$nyears)
  Ratio.init <- matrix(NA, nrow=mod.constants$nyears, ncol=3)
  phi.b.init <- runif(mod.constants$nyears,0.8,1)
  phi.nb.init <- matrix(runif(mod.constants$nyears*2,0.8,1), ncol=2)
  F.Tot.R.init <- runif(mod.constants$nyears, 0.25,0.75)
  J.surv.init <- N.b.init <- Chicks.init <- Recruit.init <- matrix(NA, nrow=mod.constants$nyears, ncol=mod.constants$n.b.months)
  availability.init <- rep(NA, times=mod.constants$nyears)
  phi.j.init <- plogis(runif(mod.constants$nyears,mod.constants$phi.j.prior.mean-mod.constants$phi.j.prior.sd,mod.constants$phi.j.prior.mean+mod.constants$phi.j.prior.sd))
  J.surv.init[,1] <- (phi.j.init^30.5)^4
  J.surv.init[,2] <- (phi.j.init^30.5)^3
  J.surv.init[,3] <- (phi.j.init^30.5)^2
  J.surv.init[,4] <- phi.j.init^30.5
  for(i in 1:mod.constants$nyears){
    availability.init[i] <- plogis(rtruncnorm(1,mod.constants$avail.min,mod.constants$avail.max,mod.constants$avail.mean[i],mod.constants$avail.sd[i]))
  }
  productivity.init <- matrix(runif(mod.constants$nyears*mod.constants$n.b.months, 0, 4), nrow=mod.constants$nyears, ncol=mod.constants$n.b.months)
  covey.size.init <- runif(mod.constants$nyears,mod.constants$csize.min,mod.constants$csize.max)
  S.n.f.init <- runif(2,200,600)
  S.nb.init <- N.Nov.init <- matrix(NA, nrow=mod.constants$nyears, ncol=2)
  S.nb.init[1,1] <- round(S.n.f.init[1],0)
  S.nb.init[1,2] <- round(S.n.f.init[2],0)
  
  # First year
  N.b.init[1,1] <- S.nb.init[1,1] + S.nb.init[1,2]
  for(m in 2:mod.constants$n.b.months){
    N.b.init[1,m] <- rbinom(1, N.b.init[1,m-1], phi.b.init[1])
  }
  N.Nov.init[1,1] <- rbinom(1, N.b.init[1,mod.constants$n.b.months], phi.b.init[1] * phi.nb.init[1,1])
  for (k in 1:mod.constants$n.b.months){
    Chicks.init[1,k] <- rpois(1, productivity.init[1,k] * N.b.init[1,k] * 0.5)
    Recruit.init[1,k] <- rbinom(1, Chicks.init[1,k], J.surv.init[1,k] * phi.nb.init[1,2])
  }
  N.Nov.init[1,2] <- sum(Recruit.init[1,1:mod.constants$n.b.months])
  N.Nov.Tot.F.init[1] <- N.Nov.init[1,1] + N.Nov.init[1,2]
  JF.F.R.init[1] <- N.Nov.init[1,2] / N.Nov.Tot.F.init[1]
  N.Nov.Tot.init[1] <- rpois(1, N.Nov.Tot.F.init[1] / F.Tot.R.init[1])
  Ratio.init[1,1:3] <- c(F.Tot.R.init[1]*(1-JF.F.R.init[1]),F.Tot.R.init[1]*JF.F.R.init[1],(1-F.Tot.R.init[1]))
  Avail.init[1] <- rbinom(1, N.Nov.Tot.init[1], mod.constants$effort[1])
  coveys.init[1] <- rpois(1, Avail.init[1] / covey.size.init[1])
  
  # Subsequent years
  for(t in 2:mod.constants$nyears){
    for(d in 1:2){ # Looping over ages
      S.nb.init[t,d] <- rbinom(1, N.Nov.init[t-1,d], phi.nb.init[t-1,d]^5 * phi.b.init[t]^2)
    }
    N.b.init[t,1] <- S.nb.init[t,1] + S.nb.init[t,2]
    for(m in 2:mod.constants$n.b.months){
      N.b.init[t,m] <- rbinom(1, N.b.init[t,m-1], phi.b.init[t])
    }
    N.Nov.init[t,1] <- rbinom(1, N.b.init[t,mod.constants$n.b.months], phi.b.init[t] * phi.nb.init[t,1])
    for (k in 1:mod.constants$n.b.months){
      Chicks.init[t,k] <- rpois(1, productivity.init[t,k] * N.b.init[t,k] * 0.5)
      Recruit.init[t,k] <- rbinom(1, Chicks.init[t,k], J.surv.init[t,k] * phi.nb.init[t,2])
    }
    N.Nov.init[t,2] <- sum(Recruit.init[t,1:mod.constants$n.b.months])
    N.Nov.Tot.F.init[t] <- N.Nov.init[t,1] + N.Nov.init[t,2]
    JF.F.R.init[t] <- N.Nov.init[t,2] / N.Nov.Tot.F.init[t]
    N.Nov.Tot.init[t] <- rpois(1, N.Nov.Tot.F.init[t] / F.Tot.R.init[t])
    Ratio.init[t,1:3] <- c(F.Tot.R.init[t]*(1-JF.F.R.init[t]),F.Tot.R.init[t]*JF.F.R.init[t],(1-F.Tot.R.init[t]))
    Avail.init[t] <- rbinom(1, N.Nov.Tot.init[t], mod.constants$effort[t])
    coveys.init[t] <- rpois(1, Avail.init[t] / covey.size.init[t])
  }

  inits.fun <- function() list(S.n.f = S.n.f.init,
                               covey.size = runif(mod.constants$nyears,mod.constants$csize.min,mod.constants$csize.max),
                               availability.logit = qlogis(availability.init),
                               pdet.logit = rtruncnorm(1,mod.constants$pdet_min,10,mod.constants$pdet_mean,mod.constants$pdet_sd),
                               phi.j.daily.logit = runif(mod.constants$nyears,mod.constants$phi.j.prior.mean-2*mod.constants$phi.j.prior.sd,mod.constants$phi.j.prior.mean+2*mod.constants$phi.j.prior.sd),
                               phi.b.logit = qlogis(runif(mod.constants$nyears,0.5,1)),
                               phi.nb.logit = qlogis(matrix(runif(mod.constants$nyears*2,0.5,1), ncol=2)),
                               productivity.log = matrix(runif(mod.constants$nyears*mod.constants$n.b.months,-1,1), nrow=mod.constants$nyears),
                               Survival.b.hyper = runif(1,0,1),
                               phi.b.hyper.logit.sd = runif(1,0,1),
                               Survival.nb.hyper = runif(2,0,1),
                               phi.nb.hyper.logit.sd = runif(2,0,1),
                               Prod.hyper = runif(mod.constants$n.b.months,0.00001,3),
                               prod.year.hyper.log.sd = runif(mod.constants$n.b.months,0,1),
                               N.b = N.b.init,
                               Chicks = Chicks.init,
                               Recruit = Recruit.init,
                               N.Nov = N.Nov.init,
                               N.Nov.Tot = N.Nov.Tot.init,
                               S.nb = S.nb.init,
                               chicks.nest.act = ceiling(mod.data$chicks.nest.mean),
                               coveys = coveys.init,
                               Avail = Avail.init,
                               F.Tot.R = runif(mod.constants$nyears,0,1))
  
  NOBO.model <- nimbleModel(code = NOBO.model.code,
                         data = mod.data,
                         constants = mod.constants,
                         inits = inits.fun())
  NOBO.mcmc.out <- nimbleMCMC(model = NOBO.model,
                         niter = 120000, nchains = 1, nburnin = 30000, thin=25,
                         samplesAsCodaMCMC=TRUE, monitor=monitors, setSeed = seed)
  return(NOBO.mcmc.out)
}


mod.data <- list(Harv.sexage = NOBO.IPM.data$Harv.sexage,
                 chicks.nest.mean = NOBO.IPM.data$chicks.nest.mean,
                 count = NOBO.IPM.data$count,
                 CH.state = NOBO.IPM.data$CH.state)

mod.constants <- list(N.tracked = NOBO.IPM.data$N.tracked,
                      effort= NOBO.IPM.data$effort,
                      Harv.N = NOBO.IPM.data$Harv.N,
                      CH.age = NOBO.IPM.data$CH.age,
                      CH.year = NOBO.IPM.data$CH.year,
                      nyears = NOBO.IPM.data$nyears,
                      n.b.months = NOBO.IPM.data$n.breeding.months, # June-September
                      nCH = NOBO.IPM.data$nCH,
                      trackdates = NOBO.IPM.data$trackdates,
                      trackperiods = NOBO.IPM.data$trackperiods,
                      csize.mean = NOBO.IPM.data$csize.mean,
                      csize.sd = NOBO.IPM.data$csize.sd,
                      csize.min = NOBO.IPM.data$csize.min,
                      csize.max = NOBO.IPM.data$csize.max,
                      Area = NOBO.IPM.data$Area,
                      phi.j.prior.mean = NOBO.IPM.data$phi.j.prior.mean,
                      phi.j.prior.sd = NOBO.IPM.data$phi.j.prior.sd,
                      chicks.nest.sd = NOBO.IPM.data$chicks.nest.sd,
                      avail.mean = NOBO.IPM.data$avail.mean,
                      avail.sd = NOBO.IPM.data$avail.sd,
                      avail.min = NOBO.IPM.data$avail.min,
                      avail.max = NOBO.IPM.data$avail.max,
                      pdet_mean = NOBO.IPM.data$pdet_mean,
                      pdet_sd = NOBO.IPM.data$pdet_sd,
                      pdet_min = NOBO.IPM.data$pdet_min)

monitors <- c("S.nb", "covey.size", "availability", "pdet", "phi.j.daily", "J.surv",
              "phi.b", "phi.nb", "productivity", "Survival.b.hyper", "phi.b.hyper.logit.sd",
              "Survival.nb.hyper", "phi.nb.hyper.logit.sd", "Prod.hyper", "prod.year.hyper.log.sd",
              "N.b", "N.Nov", "Chicks", "Recruit", "N.Nov.Tot", "D.Nov.Tot",
              "chicks.nest.act", "coveys", "Avail","Ratio","F.Tot.R","JF.F.R")

this_cluster <- makeCluster(3) # Creating cluster

NOBO_IPM <- parLapply(cl = this_cluster, X = 1:3, 
                      fun = NOBO_MCMC_code, 
                      mod.data = mod.data,
                      mod.constants = mod.constants,
                      monitors = monitors)

stopCluster(this_cluster)