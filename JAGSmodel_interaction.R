#
# Occupancy JAGS code
#

## Definitions (adapted from Waddle et al. (2010) Ecol. Appl. 20:1467-1475)
#
# zD = occupancy state of dominant species D
# zs1 = occupancy state of subordinate species s1
# zs2 = occupancy state of subordinate species s2
# psiD = Pr(zD =1)
# psis1D = Pr(zs1=1|zD=1)
# psis1 = Pr(zs1=1|zD=0)
# psis2D = Pr(zs2=1|zD=1)
# psis2 = Pr(zs2=1|zD=0)
# pD = Pr(yD=1|zD=1)
# ps1D = Pr(ys1=1|zs1=1,zD=1)
# ps1 = Pr(ys1=1|zs1=1,zD=0)
# ps2D = Pr(ys2=1|zs2=1,zD=1)
# ps2 = Pr(ys2=1|zs2=1,zD=0)

model{

  # priors
  b0.psiD ~ dnorm(0,0.001)T(-10,10)
  b0.psis1D ~ dnorm(0,0.001)T(-10,10)
  b0.psis1 ~ dnorm(0,0.001)T(-10,10)
  b0.psis2D ~ dnorm(0,0.001)T(-10,10)
  b0.psis2 ~ dnorm(0,0.001)T(-10,10)
  b0.pD ~ dnorm(0,0.001)T(-10,10)
  b0.ps1D ~ dnorm(0,0.001)T(-10,10)
  b0.ps1 ~ dnorm(0,0.001)T(-10,10)
  b0.ps2D ~ dnorm(0,0.001)T(-10,10)
  b0.ps2 ~ dnorm(0,0.001)T(-10,10)
  sigma.psi ~  dunif(0,4)
  tau.psi <- 1/(sigma.psi*sigma.psi)
  
  for(i in 1:narea){
    error.psi[i] ~ dnorm(0, tau.psi)
  }
  
  # tranformations
  for(i in 1:nSites){
    logit(psiD[i]) <- b0.psiD + error.psi[areaID[i]] # include random effects of study area on dominant species occupancy
    logit(psis1D[i]) <- b0.psis1D                        
    logit(psis1[i]) <- b0.psis1 + error.psi[areaID[i]] # include random effects of study area on subordinate species occupancy
    logit(psis2D[i]) <- b0.psis2D                        
    logit(psis2[i]) <- b0.psis2 + error.psi[areaID[i]] # include random effects of study area on subordinate species occupancy
    
    for( j in 1:nocc ){
      logit(pD[i,j]) <- b0.pD
      logit(ps1D[i,j]) <- b0.ps1D
      logit(ps1[i,j]) <- b0.ps1
      logit(ps2D[i,j]) <- b0.ps2D
      logit(ps2[i,j]) <- b0.ps2
    }
  }
  
  
  # likelihood
  # process model
  
  for (i in 1:nSites){
    zD[i] ~ dbern( psiD[i] )
    zs1[i] ~ dbern( zD[i] * psis1D[i] + (1-zD[i]) * psis1[i] )
    zs2[i] ~ dbern( zD[i] * psis2D[i] + (1-zD[i]) * psis2[i] )
    
    # observation model
    for(j in 1:nocc){
      yD[i,j] ~ dbern( zD[i] * pD[i,j] )
      ys1[i,j] ~ dbern( zs1[i] * ( zD[i] * ps1D[i,j] + 
                              (1-zD[i]) * ps1[i,j] )
                      )
    }
  }
  
  
  #fit <- sum(Presi[,])# Discrepancy for actual data set
  #fit.new <- sum(Presi.new[,]) 		# Discrepancy for replicate data set
  #ZZBa <- sum(zBa[])
  # derived parameters
  phis1D <- psiD[1]*psis1D[1]/(psiD[1]*(psiD[1]*psis1D[1]+(1-psiD[1])*psis1[1]))
  phis2D <- psiD[1]*psis2D[1]/(psiD[1]*(psiD[1]*psis2D[1]+(1-psiD[1])*psis2[1]))
  
}
