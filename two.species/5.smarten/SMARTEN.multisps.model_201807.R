#
# Occupancy JAGS code
#

model{
  #Smarten priors
  b0.psi ~ dnorm(-0.655378282,0.4536812)
  b0.p ~ dnorm(-2.849681465,1/0.119641457)
  sigma ~  dunif(0,4)
  tau <- 1/(sigma*sigma) 
  
  for(i in 1:narea){
    error[i] ~ dnorm(0, tau)
  }
  
  #Genet priors
  b0.psis1D ~ dnorm(0,0.001)T(-10,10)
  b0.psis1 ~ dnorm(0,0.001)T(-10,10)
  b0.ps1 ~ dnorm(0,0.001)T(-10,10)
  b1.ps1 ~ dnorm(0,0.001)T(-10,10)
  
  # tranformations
  
  for( i in 1:nSites ){
    #smarten
    logit(psi[i]) <- b0.psi + error[areaID[i]]
    
    
    #genet
    logit(psis1D[i]) <- b0.psis1D  
    logit(psis1[i]) <- b0.psis1
    
    for( j in 1:nocc ){
      #smarten
      logit(p[i,j]) <- b0.p
      #genet
      logit(ps1[i,j]) <- b0.ps1 + b1.ps1*season[i]
    }
  }
  
  # likelihood
  for( i in 1:nSites ){
    # process model
    z[i] ~ dbern( psi[i] )
    zs1[i] ~ dbern( z[i] * psis1D[i] + (1-z[i]) * psis1[i] )
    
    # observation model
    for( j in 1:nocc ){
      y[i,j] ~ dbern( z[i]*p[i,j] )
      ys1[i,j] ~ dbern( zs1[i] * ps1 [i,j] )

    }
  }

  # derived parameters
  phi1 <- psi[1]*psis1D[1]/(psi[1]*(psi[1]*psis1D[1]+(1-psi[1])*psis1[1]))

}