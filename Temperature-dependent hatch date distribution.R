1)	R script
#----------------------
# 1: Loading data and libraries
#----------------------

remove(list=ls())
library(R2jags)
load("DB2.RData")

#This file contains the monthly population-averaged GSI database (Table 1).

load("Malta.RData")
load("Balears.RData")

#These files contain the DB1 databases either for Malta and Balearic Islands cohorts. #DB3 database have been incorporated in these files as the sinusoidal function have been #fitted to the thermal regime for each cohort. 

Load(“data.jags.RData”)

#data.jags is the input list for the model. It contains the databases DB1, 2 and 3 ordered #in the following elements:

$n.pops.gsi Number of populations of DB2
$n is the number of data available from the GSI database
$ID.pop are the IDs for the above-mentioned populations
$GSI contain the GSI values (DB2) for each available population (7 populations)
$day.gsi This element contains the mean Julian days for each month of the year, repeated #7 times, one for each population of DB2.
$pops.birth.from First number for DB1 ID populations
$capture Julian capture dates for each DB1 individual
$age.min Julian days when each individual is vulnerable to the fishery (attained >= 20 cm)
$y Determined Julian birth date for each DB1 individual
$ pops.birth.to Last ID number for DB1 ID populations
$n.birthdates Number of birthdates data (length $y element)
$n.pops Total number of populations (cohorts) considered for the model
$sin.min
$sin.amp
$sin.phi 
# These last 3 elements contain the temperature sinusoidal parameters for each cohort #analyzed (18 considering all DB2 and DB1 populations).
$m Mortality rate
$pi pi number

#----------------------
# 2: Model
#----------------------

sink("model.txt")
cat(
  "model {#linea 1
  # Birthday submodel
  for (i in pops.birth.from:pops.birth.to){           # population level
    mu.m[i] <- mu[i] + sd[i]^2*m                      # mortality bias (mu is the true mean; mu.m is the oberved mean after mortality)                  
  }
  for (i in 1:n.birthdates){
    ones[i] ~ dbern(p.birth[i])
    p.birth[i]<-dnorm(y[i],mu.m[ID.pop.birth[i]],tol[ID.pop.birth[i]])
         / pnorm(capture[i]-age.min[i],mu.m[ID.pop.birth[i]],tol[ID.pop.birth[i]])
  }

  # GSI submodel
  for (i in 1:n){ #observation level (ID.pop[i] ranges from 1 to 7)
    GSI[i] ~ dnorm(GSI.hat[i],tol.gsi[ID.pop[i]])
    p[i] <-  (1/sqrt(2*pi)/sd[ID.pop[i]]) *
       exp(-((day.gsi[i] - mu[ID.pop[i]])^2/(2*sd[ID.pop[i]]^2))^pow)  # super gaussian distribution 
    GSI.hat[i] <- GSI.scale[ID.pop[i]]*p[i]
  }
  pow~dgamma(0.01,0.01)
  for (i in 1:n.pops.gsi){
    tol.gsi[i] ~ dgamma(0.01,0.01)    # population-specific tolerance
    GSI.scale[i]~dnorm(0,1E-10)T(0,)  # population-specific scale
  }

  # Population submodel (all pops) Linear combination from sinusoidal parameters
  for (i in 1:n.pops){  
    mu.hat[i] <- exp(p.mu[1]+
                 #p.mu[2]*sin.mean[i]+
                 #p.mu[3]*sin.amp[i]+
                 p.mu[4]*sin.phi[i])
    sd.hat[i] <- exp(p.sd[1] +
                 p.sd[2]*sin.mean[i])#+
                 #p.sd[3]*sin.amp[i]+
                 #p.sd[4]*sin.phi[i])
    mu[i] ~ dnorm(mu.hat[i],tol.mu)
    sd[i] ~ dgamma(sd.hat[i]*rate.sd,rate.sd)
    tol[i] <- 1/sd[i]^2
  }
  rate.sd ~ dgamma(0.01,0.01)
  tol.mu~dgamma(1,0.01)
  for (i in 1:4){
    p.mu[i] ~ dnorm(0,1E-10)
    p.sd[i] ~ dnorm(0,1E-10)
  }

}",fill = TRUE)
sink()

#----------------------
# 3: Settings and running JAGS
#----------------------
sd.ini=rep(10,18)
sd.ini[c(1,4,5,7)]=100
inits <- function() list(sd=sd.ini,
                         mu=rep(365/2,data.jags$n.pops),
                         p.mu=c(5,0,0,-5))
# MCMC  settings
ni <- 100
nt <- 1
nb <- 10
nc <- 3
# Parameters monitored
params <- c("p.mu","p.sd",
            "mu","sd","tol.mu","rate.sd",
            "mu.hat","sd.hat",
            "GSI.hat","tol.gsi",
            "GSI.scale","pow")
# Running
out <- jags(data.jags, inits, params, "model.txt", n.chains = nc, 
            n.thin = nt, n.iter = ni, n.burnin = nb)
time0=Sys.time()
out=update(out,n.iter =  10000, n.thin = 10)
out=update(out,n.iter =  100000, n.thin = 100)
Sys.time()-time0
n.thin = nt, n.iter = ni, n.burnin = nb, jags.seed = 123)
time0=Sys.time()
out=update(out, n.iter = 100000, n.thin = 10)
Sys.time()-time0

# save(out, file = "out.hatch.RData")
# Evaluation of parameters convergence

t(apply(out$BUGSoutput$sims.list$p.mu,2,quantile,c(0.025,0.5,0.975)))
t(apply(out$BUGSoutput$sims.list$p.sd,2,quantile,c(0.025,0.5,0.975)))

traceplot(out, varname="p.mu")
traceplot(out, varname="mu")
traceplot(out, varname="mu.hat")
traceplot(out, varname="tol.mu")
traceplot(out, varname="p.sd")
traceplot(out, varname="sd.hat")
traceplot(out, varname="pow")
traceplot(out, varname="GSI.scale")
