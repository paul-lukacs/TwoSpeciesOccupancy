#
# Occupancy JAGS code
#

model{
  # missing data
  #for(i in 1:nSites){
  #rabbt.ts[i] ~ dnorm(0,1)
  #}
  
  # priors
  b0.psiA ~ dnorm(0,0.001)T(-10,10)
  b0.psiBA ~ dnorm(0,0.001)T(-10,10)
  b0.psiBa ~ dnorm(0,0.001)T(-10,10)
  b0.pA ~ dnorm(0,0.001)T(-10,10)
  b0.pB ~ dnorm(0,0.001)T(-10,10)
  b0.rBA ~ dnorm(0,0.001)T(-10,10)
  b0.rBa ~ dnorm(0,0.001)T(-10,10)
  
  
  
  # tranformations
  for(i in 1:nSites){
    logit(psiA[i]) <- b0.psiA
    logit(psiBA[i]) <- b0.psiBA
    logit(psiBa[i]) <- b0.psiBa
    
    for( j in 1:nocc ){
      logit(pA[i,j]) <- b0.pA
      logit(pB[i,j]) <- b0.pB
      logit(rBA[i,j]) <- b0.rBA
      logit(rBa[i,j]) <- b0.rBa
    }
  }
  
  
  
  # likelihood
  # process model
  
  for (i in 1:nSites){
    zA[i] ~ dbern( psiA[i] )
    zBA[i] ~ dbern( zA[i]*psiBA[i] )
    zBa[i] ~ dbern( (1-zA[i])*psiBa[i] )
    
    # observation model
    for(j in 1:nocc){
      yA[i,j] ~ dbern(zA[i]*pA[i,j])
      yB[i,j] ~ dbern( 
                        zBa[i]* pB[i,j] + 
                        zBA[i]* ( rBA[i,j]*yA[i,j] + rBa[i,j]*(1-yA[i,j]) )
                      )
    }
  }
  
  
  #fit <- sum(Presi[,])# Discrepancy for actual data set
  #fit.new <- sum(Presi.new[,]) 		# Discrepancy for replicate data set
  
  # derived parameters
  phi <- psiA[1]*psiBA[1]/(psiA[1]*(psiA[1]*psiBA[1]+(1-psiA[1])*psiBa[1]))
  
}
