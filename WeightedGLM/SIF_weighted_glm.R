#
# SIF weighted GLM
#
#
  
  dat <- read.csv( here::here( "WeightedGLM", "SIF-DIET_data.csv" ) )
  
  sif <- dat$mean.SIF
  weights <- dat$sd.SIF
  pianka <- ( dat$obs.pianka - mean( dat$obs.pianka) )/(sd( dat$obs.pianka ) )
  bmass <- (dat$biomass.ratio - mean( dat$biomass.ratio) )/(sd( dat$biomass.ratio ) )
  
  sif.parms <- c( "b0", "b1", "b2", "b3", "b4", "sigma", "tau" )
  sif.data <- list(
    sif = sif,
    weights = weights,
    pianka = pianka,
    bmass = bmass,
    nobs = nrow(dat)
  )
  sif.inits <- function(){
    list(
      sigma = runif( 1, 0, 1 ),
      b0 = runif( 1, -3, 3 ),
      b1 = runif( 1, -3, 3 ),
      b2 = runif( 1, -3, 3 ),
      b3 = runif( 1, -3, 3 ),
      b4 = runif( 1, -3, 3 )
    )
  }
  ni <- 100000
  nb <- 10000
  nc <- 3
  nt <- 5
  
  sif.out <- R2jags::jags( 
                            sif.data, 
                            sif.inits,
                            sif.parms,
                            "WeightedGLM/SIF_glm.txt",
                            n.chains=nc, 
                            n.iter=ni, 
                            n.burnin=nb,
                            n.thin=nt
                          )
  sif.out