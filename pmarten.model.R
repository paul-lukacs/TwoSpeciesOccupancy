#
# Occupancy JAGS code
#

model{
  # priors
  b0.psi ~ dnorm(0,0.001)T(-10,10)
  b1.psi ~ dnorm(0,0.001)T(-10,10)
  b2.psi ~ dnorm(0,0.001)T(-10,10)
  b3.psi ~ dnorm(0,0.001)T(-10,10)
  b4.psi ~ dnorm(0,0.001)T(-10,10)
  b5.psi ~ dnorm(0,0.001)T(-10,10)
  b6.psi ~ dnorm(0,0.001)T(-10,10)
  b0.p ~ dnorm(0,0.001)T(-10,10)
  b1.p ~ dnorm(0,0.001)T(-10,10)
  sigma ~  dunif(0,4)
  tau <- 1/(sigma*sigma) 
  
  for(i in 1:narea){
    error[i] ~ dnorm(0, tau)
  }
  
  
  # tranformations
  for( i in 1:nSites ){
    logit(psi[i]) <- b0.psi + b1.psi*EVI[i] + b2.psi*lat[i] + b3.psi*rod.d[i] + b4.psi*rod.ts[i] + b5.psi*d.sett[i] + b6.psi*season[i] + error[areaID[i]]
    for( j in 1:nocc ){
      logit(p[i,j]) <- b0.p + b1.p*trail[i]
    }
  }
  
  # likelihood
  for( i in 1:nSites ){
    # process model
    z[i] ~ dbern( psi[i] )
    
    # observation model
    for( j in 1:nocc ){
      y[i,j] ~ dbern( z[i]*p[i,j] )
      
      #Computation of fit statistic (for Bayesian p-value)
      Presi[i,j] <- abs(y[i,j]-p[i,j])	 # Absolute residual
      y.new[i,j]~dbern(z[i]*p[i,j])
      Presi.new[i,j] <- abs(y.new[i,j]-p[i,j])
    }
  }
  
  fit <- sum(Presi[,])# Discrepancy for actual data set
  fit.new <- sum(Presi.new[,]) 		# Discrepancy for replicate data set
  
  # derived parameters
  totOcc <- sum(z[])
  logit(mean.psi) <- b0.psi
  logit(mean.p) <- b0.p
  
}
