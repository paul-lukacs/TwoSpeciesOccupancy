#
# Occupancy JAGS code
#

model{
  #Lynx priors
  b0.psi ~ dnorm(-2.926760896,1/0.421127666)
  b0.p ~ dnorm(-3.091837555,1/0.39975512)
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
  
  #Badger priors
  b0.psis3D ~ dnorm(0,0.001)T(-10,10)
  b1.psis3D ~ dnorm(0,0.001)T(-10,10)
  b2.psis3D ~ dnorm(0,0.001)T(-10,10)
  b3.psis3D ~ dnorm(0,0.001)T(-10,10)
  b0.psis3 ~ dnorm(0,0.001)T(-10,10)
  b1.psis3 ~ dnorm(0,0.001)T(-10,10)
  b2.psis3 ~ dnorm(0,0.001)T(-10,10)
  b3.psis3 ~ dnorm(0,0.001)T(-10,10)
  b0.ps3 ~ dnorm(0,0.001)T(-10,10)
  b1.ps3 ~ dnorm(0,0.001)T(-10,10)
  
  #Smarten priors
  b0.psis4D ~ dnorm(0,0.001)T(-10,10)
  b0.psis4 ~ dnorm(0,0.001)T(-10,10)
  b0.ps4 ~ dnorm(0,0.001)T(-10,10)
  
  #Genet priors
  b0.psis6D ~ dnorm(0,0.001)T(-10,10)
  b0.psis6 ~ dnorm(0,0.001)T(-10,10)
  b0.ps6 ~ dnorm(0,0.001)T(-10,10)
  b1.ps6 ~ dnorm(0,0.001)T(-10,10)
  
  #Mongoose priors
  b0.psis7D ~ dnorm(0,0.001)T(-10,10)
  b0.psis7 ~ dnorm(0,0.001)T(-10,10)
  b0.ps7 ~ dnorm(0,0.001)T(-10,10)
  
  # tranformations
  
  for( i in 1:nSites ){
    #Lynx
    logit(psi[i]) <- b0.psi + error[areaID[i]]
    
    #fox
    logit(psis1D[i]) <- b0.psis1D + b1.psis1D*EVI[i] + b2.psis1D*(EVI[i]^2)
    logit(psis1[i]) <- b0.psis1 + b1.psis1*EVI[i] + b2.psis1*(EVI[i]^2)
    
    #wildcat
    logit(psis2D[i]) <- b0.psis2D  
    logit(psis2[i]) <- b0.psis2 
    
    #badger
    logit(psis3D[i]) <- b0.psis3D + b1.psis3D*rbbt.ts[i] + b2.psis3D*d.sett[i] + b3.psis3D*season[i]
    logit(psis3[i]) <- b0.psis3 + b1.psis3*rbbt.ts[i] + b2.psis3*d.sett[i] + b3.psis3*season[i]
    
    #smarten
    logit(psis4D[i]) <- b0.psis4D  
    logit(psis4[i]) <- b0.psis4
    
    #genet
    logit(psis6D[i]) <- b0.psis6D  
    logit(psis6[i]) <- b0.psis6
    
    #mongoose
    logit(psis7D[i]) <- b0.psis7D  
    logit(psis7[i]) <- b0.psis7
    
    for( j in 1:nocc ){
      #lynx
      logit(p[i,j]) <- b0.p
      #fox
      logit(ps1[i,j]) <- b0.ps1 + b1.ps1*season[i]
      #wildcat
      logit(ps2[i,j]) <- b0.ps2
      #badger
      logit(ps3[i,j]) <- b0.ps3 + b1.ps3*season[i]
      #smarten
      logit(ps4[i,j]) <- b0.ps4 
      #genet
      logit(ps6[i,j]) <- b0.ps6 + b1.ps6*season[i]
      #mongoose
      logit(ps7[i,j]) <- b0.ps7
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
    zs6[i] ~ dbern( z[i] * psis6D[i] + (1-z[i]) * psis6[i] )
    zs7[i] ~ dbern( z[i] * psis7D[i] + (1-z[i]) * psis7[i] )
    
    # observation model
    for( j in 1:nocc ){
      y[i,j] ~ dbern( z[i]*p[i,j] )
      ys1[i,j] ~ dbern( zs1[i] * ps1 [i,j] )
      ys2[i,j] ~ dbern( zs2[i] * ps2 [i,j] )
      ys3[i,j] ~ dbern( zs3[i] * ps3 [i,j] )
      ys4[i,j] ~ dbern( zs4[i] * ps4 [i,j] )
      ys6[i,j] ~ dbern( zs6[i] * ps6 [i,j] )
      ys7[i,j] ~ dbern( zs7[i] * ps7 [i,j] )

    }
  }

  # derived parameters
  phi1 <- psi[1]*psis1D[1]/(psi[1]*(psi[1]*psis1D[1]+(1-psi[1])*psis1[1]))
  phi2 <- psi[1]*psis2D[1]/(psi[1]*(psi[1]*psis2D[1]+(1-psi[1])*psis2[1]))
  phi3 <- psi[1]*psis3D[1]/(psi[1]*(psi[1]*psis3D[1]+(1-psi[1])*psis3[1]))
  phi4 <- psi[1]*psis4D[1]/(psi[1]*(psi[1]*psis4D[1]+(1-psi[1])*psis4[1]))
  phi6 <- psi[1]*psis6D[1]/(psi[1]*(psi[1]*psis6D[1]+(1-psi[1])*psis6[1]))
  phi7 <- psi[1]*psis7D[1]/(psi[1]*(psi[1]*psis7D[1]+(1-psi[1])*psis7[1]))
}