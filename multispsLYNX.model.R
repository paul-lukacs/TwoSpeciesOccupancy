#
# Occupancy JAGS code
#

## Definitions (Waddle et al. (2010) Ecol. Appl. 20:1467-1475)
#
# zB = occupancy state of dominant species
# zA = occupancy state of subordinate species
# psiB = Pr(zB =1)
# psiAB = Pr(zA=1|zB=1)
# psi A = Pr(zA=1|zB=0)
# pB = Pr(yB=1|zB=1)
# pAB = Pr(yA=1|zA=1,zB=1)
# pA = Pr(yA=1|zA=1,zB=0)

model{
  #missing data
  for(i in 1:nSites){
    rabbit.ts[i] ~ dnorm(0,1)
  }
  
  ## priors
  #Lynx
  b0.psiD ~ dnorm(0,0.01)T(-10,10)
  b0.pD ~ dnorm(0,0.01)T(-10,10)
#  sigma.psiD ~  dunif(0,4)
#  tau.psiD <- 1/(sigma.psiD*sigma.psiD)
#  for(i in 1:narea){error.psiD[i] ~ dnorm(0, tau.psiD)
#  }
  
  #Badger
  b0.psis1D ~ dnorm(0,0.01)T(-10,10)
  b0.psis1 ~ dnorm(0,0.01)T(-10,10)
  b1.psis1 ~ dnorm(0,0.01)T(-10,10)
  b0.ps1D ~ dnorm(0,0.01)T(-10,10)
  b0.ps1 ~ dnorm(0,0.01)T(-10,10)
  sigma.psis1 ~  dunif(0,4)
  tau.psis1 <- 1/(sigma.psis1*sigma.psis1)
  for(i in 1:narea){error.psis1[i] ~ dnorm(0, tau.psis1)
  }
  
  #Fox
  b0.psis2D ~ dnorm(0,0.01)T(-10,10)
  b0.psis2 ~ dnorm(0,0.01)T(-10,10)
  b1.psis2 ~ dnorm(0,0.01)T(-10,10)
  b2.psis2 ~ dnorm(0,0.01)T(-10,10)
  b0.ps2D ~ dnorm(0,0.01)T(-10,10)
  b0.ps2 ~ dnorm(0,0.01)T(-10,10)
  b1.ps2 ~ dnorm(0,0.01)T(-10,10)
  sigma.psis2 ~  dunif(0,4)
  tau.psis2 <- 1/(sigma.psis2*sigma.psis2)
  for(i in 1:narea){error.psis2[i] ~ dnorm(0, tau.psis2)
  }
  
  #Wildcat
  b0.psis3D ~ dnorm(0,0.01)T(-10,10)
  b0.psis3 ~ dnorm(0,0.01)T(-10,10)
  b0.ps3D ~ dnorm(0,0.01)T(-10,10)
  b1.psis3 ~ dnorm(0,0.01)T(-10,10)
  b2.psis3 ~ dnorm(0,0.01)T(-10,10)
  b0.ps3 ~ dnorm(0,0.01)T(-10,10) 
  sigma.psis3 ~  dunif(0,4)
  tau.psis3 <- 1/(sigma.psis3*sigma.psis3)
  for(i in 1:narea){error.psis3[i] ~ dnorm(0, tau.psis3)
  }
  
  #smarten
  b0.psis4D ~ dnorm(0,0.01)T(-10,10)
  b0.psis4 ~ dnorm(0,0.01)T(-10,10)
  b1.psis4 ~ dnorm(0,0.01)T(-10,10)
  b0.ps4D ~ dnorm(0,0.01)T(-10,10)
  b0.ps4 ~ dnorm(0,0.01)T(-10,10) 
  b1.ps4 ~ dnorm(0,0.01)T(-10,10) 
  sigma.psis4 ~  dunif(0,4)
  tau.psis4 <- 1/(sigma.psis4*sigma.psis4)
  for(i in 1:narea){error.psis4[i] ~ dnorm(0, tau.psis4)
  }
  
  #Genet
  b0.psis5D ~ dnorm(0,0.01)T(-10,10)
  b0.psis5 ~ dnorm(0,0.01)T(-10,10)
  b0.ps5D ~ dnorm(0,0.01)T(-10,10)
  b0.ps5 ~ dnorm(0,0.01)T(-10,10)
  sigma.psis5 ~  dunif(0,4)
  tau.psis5 <- 1/(sigma.psis5*sigma.psis5)
  for(i in 1:narea){error.psis5[i] ~ dnorm(0, tau.psis5)
  }
  
  ## tranformations
  # Occupancy
  for(i in 1:nSites){
    #lynx
    logit(psiD[i]) <- b0.psiD
    #badger
    logit(psis1D[i]) <- b0.psis1D + b1.psis1*d.sett[i]
    logit(psis1[i]) <- b0.psis1 + b1.psis1*d.sett[i] + error.psis1[areaID[i]]
    #fox
    logit(psis2D[i]) <- b0.psis2D + b1.psis2*EVI[i] + b1.psis2*(EVI[i]^2)
    logit(psis2[i]) <- b0.psis2 + b1.psis2*EVI[i] + b1.psis2*(EVI[i]^2) + error.psis2[areaID[i]]
    #wildcat
    logit(psis3D[i]) <- b0.psis3D + b1.psis3*rabbit.ts[i] + b2.psis3*d.sett[i]
    logit(psis3[i]) <- b0.psis3 + b1.psis3*rabbit.ts[i] + b2.psis3*d.sett[i] + error.psis3[areaID[i]]
    #smarten
    logit(psis4D[i]) <- b0.psis4D + b1.psis4*rod.ts[i]
    logit(psis4[i]) <- b0.psis4 + b1.psis4*rod.ts[i] + error.psis4[areaID[i]]
    #genet
    logit(psis5D[i]) <- b0.psis5D
    logit(psis5[i]) <- b0.psis5 + error.psis5[areaID[i]]
    
    #Detection
    for( j in 1:nocc ){
      #lynx
      logit(pD[i,j]) <- b0.pD
      #badger
      logit(ps1D[i,j]) <- b0.ps1D
      logit(ps1[i,j]) <- b0.ps1
      #fox
      logit(ps2D[i,j]) <- b0.ps2D + b1.ps2*ctr[i,j]
      logit(ps2[i,j]) <- b0.ps2 + b1.ps2*ctr[i,j]
      #wildcat
      logit(ps3D[i,j]) <- b0.ps3D
      logit(ps3[i,j]) <- b0.ps3
      #smarten
      logit(ps4D[i,j]) <- b0.ps4D + b1.ps4*ctr[i,j]
      logit(ps4[i,j]) <- b0.ps4 + b1.ps4*ctr[i,j]
      #genet
      logit(ps5D[i,j]) <- b0.ps5D
      logit(ps5[i,j]) <- b0.ps5
    }
  }
  
  
  # likelihood
  # process model
  
  for (i in 1:nSites){
    zD[i] ~ dbern( psiD[i] )
    zs1[i] ~ dbern( zD[i] * psis1D[i] + (1-zD[i]) * psis1[i] )
    zs2[i] ~ dbern( zD[i] * psis2D[i] + (1-zD[i]) * psis2[i] )
    zs3[i] ~ dbern( zD[i] * psis3D[i] + (1-zD[i]) * psis3[i] )
    zs4[i] ~ dbern( zD[i] * psis4D[i] + (1-zD[i]) * psis4[i] )
    zs5[i] ~ dbern( zD[i] * psis5D[i] + (1-zD[i]) * psis5[i] )
    
    # observation model
    for(j in 1:nocc){
      yD[i,j] ~ dbern( zD[i] * pD[i,j] )
      ys1[i,j] ~ dbern( zs1[i] * ( zD[i] * ps1D[i,j] + (1-zD[i]) * ps1[i,j] ) )
      ys2[i,j] ~ dbern( zs2[i] * ( zD[i] * ps2D[i,j] + (1-zD[i]) * ps2[i,j] ) )
      ys3[i,j] ~ dbern( zs3[i] * ( zD[i] * ps3D[i,j] + (1-zD[i]) * ps3[i,j] ) )
      ys4[i,j] ~ dbern( zs4[i] * ( zD[i] * ps4D[i,j] + (1-zD[i]) * ps4[i,j] ) )
      ys5[i,j] ~ dbern( zs5[i] * ( zD[i] * ps5D[i,j] + (1-zD[i]) * ps5[i,j] ) )
    }
  }

  # derived parameters
  phis1D <- psiD[1]*psis1D[1]/(psiD[1]*(psiD[1]*psis1D[1]+(1-psiD[1])*psis1[1]))
  phis2D <- psiD[1]*psis2D[1]/(psiD[1]*(psiD[1]*psis2D[1]+(1-psiD[1])*psis2[1]))
  phis3D <- psiD[1]*psis3D[1]/(psiD[1]*(psiD[1]*psis3D[1]+(1-psiD[1])*psis3[1]))
  phis4D <- psiD[1]*psis4D[1]/(psiD[1]*(psiD[1]*psis4D[1]+(1-psiD[1])*psis4[1]))
  phis5D <- psiD[1]*psis5D[1]/(psiD[1]*(psiD[1]*psis5D[1]+(1-psiD[1])*psis5[1]))

}