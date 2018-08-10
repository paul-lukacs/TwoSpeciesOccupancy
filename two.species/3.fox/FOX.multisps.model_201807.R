#
# Occupancy JAGS code
#

model{
  
  
  #Fox priors
  b0.psi ~ dnorm(0.049654631,1/0.402275665)T(-0.780184595,0.84234527)
  b1.psi ~ dnorm(0.415183357,1/0.22420182)T(-0.027902466,0.857708284)
  b2.psi ~ dnorm(-0.302511439,1/0.106021655)T(-0.513166452,-0.098394659)
  b0.p ~ dnorm(-1.814893453,1/0.057757601)T(-1.929433438,-1.704332405)
  b1.p ~ dnorm(0.173064624,1/0.075714471)T(0.02525486,0.322775335)
  sigma ~  dunif(0,4)
  tau <- 1/(sigma*sigma) 
  
  for(i in 1:narea){
    error[i] ~ dnorm(0, tau)
  }
  #Wildcat priors
  b0.psis1D ~ dnorm(0,0.001)T(-10,10)
  b0.psis1 ~ dnorm(0,0.001)T(-10,10)
  b0.ps1 ~ dnorm(0,0.001)T(-10,10)
  
  #Pmarten priors
  b0.psis2D ~ dnorm(0,0.001)T(-10,10)
  b0.psis2 ~ dnorm(0,0.001)T(-10,10)
  b0.ps2 ~ dnorm(0,0.001)T(-10,10)
  b1.ps2 ~ dnorm(0,0.001)T(-10,10)
  
  #Smarten priors
  b0.psis3D ~ dnorm(0,0.001)T(-10,10)
  b0.psis3 ~ dnorm(0,0.001)T(-10,10)
  b0.ps3 ~ dnorm(0,0.001)T(-10,10)
  
  #Mongoose priors
  b0.psis4D ~ dnorm(0,0.001)T(-10,10)
  b0.psis4 ~ dnorm(0,0.001)T(-10,10)
  b0.ps4 ~ dnorm(0,0.001)T(-10,10)
  
  #Genet priors
  b0.psis5D ~ dnorm(0,0.001)T(-10,10)
  b0.psis5 ~ dnorm(0,0.001)T(-10,10)
  b0.ps5 ~ dnorm(0,0.001)T(-10,10)
  b1.ps5 ~ dnorm(0,0.001)T(-10,10)
  
  # tranformations
  
  for( i in 1:nSites ){
    #fox
    logit(psi[i]) <- b0.psi + b1.psi*EVI[i] + b2.psi*(EVI[i]^2) + error[areaID[i]]
    
    #wildcat
    logit(psis1D[i]) <- b0.psis1D
    logit(psis1[i]) <- b0.psis1
    
    #pmarten
    logit(psis2D[i]) <- b0.psis2D  
    logit(psis2[i]) <- b0.psis2
    
    #smarten
    logit(psis3D[i]) <- b0.psis3D  
    logit(psis3[i]) <- b0.psis3
    
    #mongoose
    logit(psis4D[i]) <- b0.psis4D  
    logit(psis4[i]) <- b0.psis4
    
    #genet
    logit(psis5D[i]) <- b0.psis5D  
    logit(psis5[i]) <- b0.psis5
    
    for( j in 1:nocc ){
      #fox
      logit(p[i,j]) <- b0.p + b1.p*season[i]
      #wildcat
      logit(ps1[i,j]) <- b0.ps1
      #pmarten
      logit(ps2[i,j]) <- b0.ps2 + b1.ps2*trail[i]
      #smarten
      logit(ps3[i,j]) <- b0.ps3
      #mongoose
      logit(ps4[i,j]) <- b0.ps4
      #genet
      logit(ps5[i,j]) <- b0.ps5 + b1.ps5*season[i]
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
    zs5[i] ~ dbern( z[i] * psis5D[i] + (1-z[i]) * psis5[i] )
    
    # observation model
    for( j in 1:nocc ){
      y[i,j] ~ dbern( z[i]*p[i,j] )
      ys1[i,j] ~ dbern( zs1[i] * ps1 [i,j] )
      ys2[i,j] ~ dbern( zs2[i] * ps2 [i,j] )
      ys3[i,j] ~ dbern( zs3[i] * ps3 [i,j] )
      ys4[i,j] ~ dbern( zs4[i] * ps4 [i,j] )
      ys5[i,j] ~ dbern( zs5[i] * ps5 [i,j] )

    }
  }

  # derived parameters
  phi1 <- psi[1]*psis1D[1]/(psi[1]*(psi[1]*psis1D[1]+(1-psi[1])*psis1[1]))
  phi2 <- psi[1]*psis2D[1]/(psi[1]*(psi[1]*psis2D[1]+(1-psi[1])*psis2[1]))
  phi3 <- psi[1]*psis3D[1]/(psi[1]*(psi[1]*psis3D[1]+(1-psi[1])*psis3[1]))
  phi4 <- psi[1]*psis4D[1]/(psi[1]*(psi[1]*psis4D[1]+(1-psi[1])*psis4[1]))
  phi5 <- psi[1]*psis5D[1]/(psi[1]*(psi[1]*psis5D[1]+(1-psi[1])*psis5[1]))

}