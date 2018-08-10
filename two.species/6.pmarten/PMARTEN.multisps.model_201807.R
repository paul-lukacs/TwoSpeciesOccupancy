#
# Occupancy JAGS code
#

model{
  #Pmarten priors
  b0.psi ~ dnorm(-1.528000001,1/0.503237071)#T(-10,10)
  b0.p ~ dnorm(-3.157240677,1/0.20306752)#T(-10,10)
  b1.p ~ dnorm(0.970305022,1/0.435174098)#T(-10,10)
  sigma ~  dunif(0,4)
  tau <- 1/(sigma*sigma) 
  
  for(i in 1:narea){
    error[i] ~ dnorm(0, tau)
  }
  
  #Smarten priors
  b0.psis1D ~ dnorm(0,0.001)T(-10,10)
  b0.psis1 ~ dnorm(0,0.001)T(-10,10)
  b0.ps1 ~ dnorm(0,0.001)T(-10,10)
  
  #Mongoose priors
  b0.psis2D ~ dnorm(0,0.001)T(-10,10)
  b0.psis2 ~ dnorm(0,0.001)T(-10,10)
  b0.ps2 ~ dnorm(0,0.001)T(-10,10)
  
  #Genet priors
  b0.psis3D ~ dnorm(0,0.001)T(-10,10)
  b0.psis3 ~ dnorm(0,0.001)T(-10,10)
  b0.ps3 ~ dnorm(0,0.001)T(-10,10)
  b1.ps3 ~ dnorm(0,0.001)T(-10,10)
  
  # tranformations
  
  for( i in 1:nSites ){
    #pmarten
    logit(psi[i]) <- b0.psi + error[areaID[i]]
    
    #smarten
    logit(psis1D[i]) <- b0.psis1D  
    logit(psis1[i]) <- b0.psis1
    
    #mongoose
    logit(psis2D[i]) <- b0.psis2D  
    logit(psis2[i]) <- b0.psis2
    
    #genet
    logit(psis3D[i]) <- b0.psis3D  
    logit(psis3[i]) <- b0.psis3
    
    for( j in 1:nocc ){
      #pmarten
      logit(p[i,j]) <- b0.p + b1.p*trail[i]
      #smarten
      logit(ps1[i,j]) <- b0.ps1
      #mongoose
      logit(ps2[i,j]) <- b0.ps2
      #genet
      logit(ps3[i,j]) <- b0.ps3 + b1.ps3*season[i]
    }
  }
  
  # likelihood
  for( i in 1:nSites ){
    # process model
    z[i] ~ dbern( psi[i] )
    zs1[i] ~ dbern( z[i] * psis1D[i] + (1-z[i]) * psis1[i] )
    zs2[i] ~ dbern( z[i] * psis2D[i] + (1-z[i]) * psis2[i] )
    zs3[i] ~ dbern( z[i] * psis3D[i] + (1-z[i]) * psis3[i] )
    
    # observation model
    for( j in 1:nocc ){
      y[i,j] ~ dbern( z[i]*p[i,j] )
      ys1[i,j] ~ dbern( zs1[i] * ps1 [i,j] )
      ys2[i,j] ~ dbern( zs2[i] * ps2 [i,j] )
      ys3[i,j] ~ dbern( zs3[i] * ps3 [i,j] )

    }
  }

  # derived parameters
  phi1 <- psi[1]*psis1D[1]/(psi[1]*(psi[1]*psis1D[1]+(1-psi[1])*psis1[1]))
  phi2 <- psi[1]*psis2D[1]/(psi[1]*(psi[1]*psis2D[1]+(1-psi[1])*psis2[1]))
  phi3 <- psi[1]*psis3D[1]/(psi[1]*(psi[1]*psis3D[1]+(1-psi[1])*psis3[1]))

}