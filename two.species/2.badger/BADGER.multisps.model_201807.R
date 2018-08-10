#
# Occupancy JAGS code
#

model{
  
  
  # Badger priors
  b0.psi ~ dnorm(-2.880671092,1/0.537751046)T(-4.023750604,-1.913091218)
  b1.psi ~ dnorm(-1.305668366,1/1.063349834)T(-4.014522854,0.125484235)
  b2.psi ~ dnorm(0.604389517,1/0.270925201)T(0.067088184,1.143198563)
  b3.psi ~ dnorm(1.10123539,1/0.895640845)T(-0.207652853,3.224275245)
  b0.p ~ dnorm(-2.661365988,1/0.207155372)T(-3.080050942,-2.274492958)
  b1.p ~ dnorm(-1.568496913,1/0.54195838)T(-2.725512582,-0.585469764)
  sigma ~  dunif(0,4)
  tau <- 1/(sigma*sigma) 
  
  for(i in 1:narea){
    error[i] ~ dnorm(0, tau)
  }
  
  #Fox priors
  b0.psis1D ~ dnorm(0,0.001)T(-10,10)
  b1.psis1D ~ dnorm(0,0.001)T(-10,10)
  b2.psis1D ~ dnorm(0,0.001)T(-10,10)
  b0.psis1 ~ dnorm(0,0.001)T(-10,10)
  b1.psis1 ~ dnorm(0,0.001)T(-10,10)
  b2.psis1 ~ dnorm(0,0.001)T(-10,10)
  b0.ps1 ~ dnorm(0,0.001)T(-10,10)
  b1.ps1 ~ dnorm(0,0.001)T(-10,10)
  
  #Wildcat priors
  b0.psis2D ~ dnorm(0,0.001)T(-10,10)
  b0.psis2 ~ dnorm(0,0.001)T(-10,10)
  b0.ps2 ~ dnorm(0,0.001)T(-10,10)
  
  #Pmarten priors
  b0.psis3D ~ dnorm(0,0.001)T(-10,10)
  b0.psis3 ~ dnorm(0,0.001)T(-10,10)
  b0.ps3 ~ dnorm(0,0.001)T(-10,10)
  b1.ps3 ~ dnorm(0,0.001)T(-10,10)
  
  #Smarten priors
  b0.psis4D ~ dnorm(0,0.001)T(-10,10)
  b0.psis4 ~ dnorm(0,0.001)T(-10,10)
  b0.ps4 ~ dnorm(0,0.001)T(-10,10)
  
  #Mongoose priors
  b0.psis5D ~ dnorm(0,0.001)T(-10,10)
  b0.psis5 ~ dnorm(0,0.001)T(-10,10)
  b0.ps5 ~ dnorm(0,0.001)T(-10,10)
  
  #Genet priors
  b0.psis6D ~ dnorm(0,0.001)T(-10,10)
  b0.psis6 ~ dnorm(0,0.001)T(-10,10)
  b0.ps6 ~ dnorm(0,0.001)T(-10,10)
  b1.ps6 ~ dnorm(0,0.001)T(-10,10)
  
  # tranformations
  
  for( i in 1:nSites ){
    #badger
    logit(psi[i]) <- b0.psi + b1.psi*rbbt.ts[i] + b2.psi*d.sett[i] + b3.psi*season[i] + error[areaID[i]]
    
    #fox
    logit(psis1D[i]) <- b0.psis1D + b1.psis1D*EVI[i] + b2.psis1D*(EVI[i]^2)
    logit(psis1[i]) <- b0.psis1 + b1.psis1*EVI[i] + b2.psis1*(EVI[i]^2)
    
    #wildcat
    logit(psis2D[i]) <- b0.psis2D
    logit(psis2[i]) <- b0.psis2
    
    #pmarten
    logit(psis3D[i]) <- b0.psis3D  
    logit(psis3[i]) <- b0.psis3
    
    #smarten
    logit(psis4D[i]) <- b0.psis4D  
    logit(psis4[i]) <- b0.psis4
    
    #mongoose
    logit(psis5D[i]) <- b0.psis5D  
    logit(psis5[i]) <- b0.psis5
    
    #genet
    logit(psis6D[i]) <- b0.psis6D  
    logit(psis6[i]) <- b0.psis6
    
    for( j in 1:nocc ){
      #badger
      logit(p[i,j]) <- b0.p + b1.p*season[i]
      #fox
      logit(ps1[i,j]) <- b0.ps1 + b1.ps1*season[i]
      #wildcat
      logit(ps2[i,j]) <- b0.ps2
      #pmarten
      logit(ps3[i,j]) <- b0.ps3 + b1.ps3*trail[i]
      #smarten
      logit(ps4[i,j]) <- b0.ps4
      #mongoose
      logit(ps5[i,j]) <- b0.ps5
      #genet
      logit(ps6[i,j]) <- b0.ps6 + b1.ps6*season[i]
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
    zs6[i] ~ dbern( z[i] * psis6D[i] + (1-z[i]) * psis6[i] )
    
    # observation model
    for( j in 1:nocc ){
      y[i,j] ~ dbern( z[i]*p[i,j] )
      ys1[i,j] ~ dbern( zs1[i] * ps1 [i,j] )
      ys2[i,j] ~ dbern( zs2[i] * ps2 [i,j] )
      ys3[i,j] ~ dbern( zs3[i] * ps3 [i,j] )
      ys4[i,j] ~ dbern( zs4[i] * ps4 [i,j] )
      ys5[i,j] ~ dbern( zs5[i] * ps5 [i,j] )
      ys6[i,j] ~ dbern( zs6[i] * ps6 [i,j] )

    }
  }

  # derived parameters
  phi1 <- psi[1]*psis1D[1]/(psi[1]*(psi[1]*psis1D[1]+(1-psi[1])*psis1[1]))
  phi2 <- psi[1]*psis2D[1]/(psi[1]*(psi[1]*psis2D[1]+(1-psi[1])*psis2[1]))
  phi3 <- psi[1]*psis3D[1]/(psi[1]*(psi[1]*psis3D[1]+(1-psi[1])*psis3[1]))
  phi4 <- psi[1]*psis4D[1]/(psi[1]*(psi[1]*psis4D[1]+(1-psi[1])*psis4[1]))
  phi5 <- psi[1]*psis5D[1]/(psi[1]*(psi[1]*psis5D[1]+(1-psi[1])*psis5[1]))
  phi6 <- psi[1]*psis6D[1]/(psi[1]*(psi[1]*psis6D[1]+(1-psi[1])*psis6[1]))

}