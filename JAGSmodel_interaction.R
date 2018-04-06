#
# Occupancy JAGS code
#

## Definitions (Waddle et al. (2010) Ecol. Appl. 20:1467-1475)
#
# zB = occupancy state of dominant species
# zA = occupancy state of subordinate species
# psiB = Pr(zB =1)
# psiAB = Pr(zA=1|zB=1)
# psi Ab = Pr(zA=1|zB=0)
# pB = Pr(yB=1|zB=1)
# pAB = Pr(yA=1|zA=1,zB=1)
# pAb = Pr(yA=1|zA=1,zB=0)

model{

  # priors
  b0.psiB ~ dnorm(0,0.001)T(-10,10)
  b0.psiAB ~ dnorm(0,0.001)T(-10,10)
  b0.psiAb ~ dnorm(0,0.001)T(-10,10)
  b0.pB ~ dnorm(0,0.001)T(-10,10)
  b0.pAB ~ dnorm(0,0.001)T(-10,10)
  b0.pAb ~ dnorm(0,0.001)T(-10,10)  
  
  
  # tranformations
  for(i in 1:nSites){
    logit(psiB[i]) <- b0.psiB
    logit(psiAB[i]) <- b0.psiAB
    logit(psiAb[i]) <- b0.psiAb
    
    for( j in 1:nocc ){
      logit(pB[i,j]) <- b0.pB
      logit(pAB[i,j]) <- b0.pAB
      logit(pAb[i,j]) <- b0.pAb
    }
  }
  
  
  # likelihood
  # process model
  
  for (i in 1:nSites){
    zB[i] ~ dbern( psiB[i] )
    zA[i] ~ dbern( zB[i] * psiAB[i] + (1-zB[i]) * psiAb[i] )
    
    # observation model
    for(j in 1:nocc){
      yB[i,j] ~ dbern( zB[i] * pB[i,j] )
      yA[i,j] ~ dbern( zA[i] * ( zB[i] * pAB[i,j] + 
                              (1-zB[i]) * pAb[i,j] )
                      )
    }
  }
  
  
  #fit <- sum(Presi[,])# Discrepancy for actual data set
  #fit.new <- sum(Presi.new[,]) 		# Discrepancy for replicate data set
  ZZBa <- sum(zBa[])
  # derived parameters
  phi <- psiB[1]*psiAB[1]/(psiB[1]*(psiB[1]*psiAB[1]+(1-psiB[1])*psiAb[1]))
  
}
