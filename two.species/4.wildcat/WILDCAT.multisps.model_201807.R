#
# Occupancy JAGS code
#

model{
  
    #Wildcat priors
  b0.psi ~ dnorm(-2.109986971,1/0.973610734)#T(-10,10)
  b0.p ~ dnorm(-3.867320893,1/0.344432181)#T(-10,10)
  sigma ~  dunif(0,4)
  tau <- 1/(sigma*sigma) 
  
  for(i in 1:narea){
    error[i] ~ dnorm(0, tau)
  }
  
  #Pmarten priors
  b0.psis1D ~ dnorm(0,0.001)T(-10,10)
  b0.psis1 ~ dnorm(0,0.001)T(-10,10)
  b0.ps1 ~ dnorm(0,0.001)T(-10,10)
  b1.ps1 ~ dnorm(0,0.001)T(-10,10)
  
  #Smarten priors
  b0.psis2D ~ dnorm(0,0.001)T(-10,10)
  b0.psis2 ~ dnorm(0,0.001)T(-10,10)
  b0.ps2 ~ dnorm(0,0.001)T(-10,10)
  
  #Mongoose priors
  b0.psis3D ~ dnorm(0,0.001)T(-10,10)
  b0.psis3 ~ dnorm(0,0.001)T(-10,10)
  b0.ps3 ~ dnorm(0,0.001)T(-10,10)
  
  #Genet priors
  b0.psis4D ~ dnorm(0,0.001)T(-10,10)
  b0.psis4 ~ dnorm(0,0.001)T(-10,10)
  b0.ps4 ~ dnorm(0,0.001)T(-10,10)
  b1.ps4 ~ dnorm(0,0.001)T(-10,10)
  
  # tranformations
  
  for( i in 1:nSites ){
    
    #wildcat
    logit(psi[i]) <- b0.psi + error[areaID[i]]
    
    #pmarten
    logit(psis1D[i]) <- b0.psis1D  
    logit(psis1[i]) <- b0.psis1
    
    #smarten
    logit(psis2D[i]) <- b0.psis2D  
    logit(psis2[i]) <- b0.psis2
    
    #mongoose
    logit(psis3D[i]) <- b0.psis3D  
    logit(psis3[i]) <- b0.psis3
    
    #genet
    logit(psis4D[i]) <- b0.psis4D  
    logit(psis4[i]) <- b0.psis4
    
    for( j in 1:nocc ){
      #wildcat
      logit(p[i,j]) <- b0.p
      #pmarten
      logit(ps1[i,j]) <- b0.ps1 + b1.ps1*trail[i]
      #smarten
      logit(ps2[i,j]) <- b0.ps2
      #mongoose
      logit(ps3[i,j]) <- b0.ps3
      #genet
      logit(ps4[i,j]) <- b0.ps4 + b1.ps4*season[i]
    }
  }
  
  # likelihood
  for( i in 1:nSites ){
    # process model
    z[i] ~ dbern( psi[i] )
    zs1[i] ~ dbern( z[i] * psis1D[i] + (1-z[i]) * psis1[i] )
    zs2[i] ~ dbern( z[i] * psis2D[i] + (1-z[i]) * psis2[i] )
    zs3[i] ~ dbern( z[i] * psis3D[i] + (1-z[i]) * psis3[i] )
    zs4[i] ~ dbern( z[i] * psis4D[i] + (1-z[i]) * psis4[i] )
    
    # observation model
    for( j in 1:nocc ){
      y[i,j] ~ dbern( z[i]*p[i,j] )
      ys1[i,j] ~ dbern( zs1[i] * ps1 [i,j] )
      ys2[i,j] ~ dbern( zs2[i] * ps2 [i,j] )
      ys3[i,j] ~ dbern( zs3[i] * ps3 [i,j] )
      ys4[i,j] ~ dbern( zs4[i] * ps4 [i,j] )

    }
  }

  # derived parameters
  phi1 <- psi[1]*psis1D[1]/(psi[1]*(psi[1]*psis1D[1]+(1-psi[1])*psis1[1]))
  phi2 <- psi[1]*psis2D[1]/(psi[1]*(psi[1]*psis2D[1]+(1-psi[1])*psis2[1]))
  phi3 <- psi[1]*psis3D[1]/(psi[1]*(psi[1]*psis3D[1]+(1-psi[1])*psis3[1]))
  phi4 <- psi[1]*psis4D[1]/(psi[1]*(psi[1]*psis4D[1]+(1-psi[1])*psis4[1]))

}