#
# Occupancy JAGS code
#

model{
  #Mongoose priors
  b0.psi ~ dnorm(-2.151260588,1/0.868145338)#T(-10,10)
  b0.p ~ dnorm(-3.581240066,1/0.428667907)#T(-10,10)
  sigma ~  dunif(0,4)
  tau <- 1/(sigma*sigma) 
  
  for(i in 1:narea){
    error[i] ~ dnorm(0, tau)
  }
  
  #Smarten priors
  b0.psis1D ~ dnorm(0,0.001)T(-10,10)
  b0.psis1 ~ dnorm(0,0.001)T(-10,10)
  b0.ps1 ~ dnorm(0,0.001)T(-10,10)
  
  #Genet priors
  b0.psis2D ~ dnorm(0,0.001)T(-10,10)
  b0.psis2 ~ dnorm(0,0.001)T(-10,10)
  b0.ps2 ~ dnorm(0,0.001)T(-10,10)
  b1.ps2 ~ dnorm(0,0.001)T(-10,10)
  
  # tranformations
  
  for( i in 1:nSites ){
    #mongoose
    logit(psi[i]) <- b0.psi + error[areaID[i]]
    
    #smarten
    logit(psis1D[i]) <- b0.psis1D  
    logit(psis1[i]) <- b0.psis1
    
    
    #genet
    logit(psis2D[i]) <- b0.psis2D  
    logit(psis2[i]) <- b0.psis2
    
    for( j in 1:nocc ){
      #mongoose
      logit(p[i,j]) <- b0.p
      #smarten
      logit(ps1[i,j]) <- b0.ps1
      #genet
      logit(ps2[i,j]) <- b0.ps2 + b1.ps2*season[i]
    }
  }
  
  # likelihood
  for( i in 1:nSites ){
    # process model
    z[i] ~ dbern( psi[i] )
    zs1[i] ~ dbern( z[i] * psis1D[i] + (1-z[i]) * psis1[i] )
    zs2[i] ~ dbern( z[i] * psis2D[i] + (1-z[i]) * psis2[i] )
    
    # observation model
    for( j in 1:nocc ){
      y[i,j] ~ dbern( z[i]*p[i,j] )
      ys1[i,j] ~ dbern( zs1[i] * ps1 [i,j] )
      ys2[i,j] ~ dbern( zs2[i] * ps2 [i,j] )

    }
  }

  # derived parameters
  phi1 <- psi[1]*psis1D[1]/(psi[1]*(psi[1]*psis1D[1]+(1-psi[1])*psis1[1]))
  phi2 <- psi[1]*psis2D[1]/(psi[1]*(psi[1]*psis2D[1]+(1-psi[1])*psis2[1]))

}