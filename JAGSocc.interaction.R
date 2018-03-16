#
# Occupancy analysis
#

library(mcmcplots)
library(R2jags)

# Check to see who's running the code
# Set the path according to the user.
if( Sys.info()["user"] == "paul.lukacs" ){
  setwd( "c:/Users/paul.lukacs/Documents/GitHub/TwoSpeciesOccupancy" )
} else {
  setwd( "/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/1.rcode/occ/fox-wildcat" )
}


#EHw <- yfs[,4:33]
EHmf <- ymf[,3:32]
EHf <- yvv[,3:32]
nocc <- ncol(EHw)
nSites <- nrow(EHw)

# random effect

areaID <- as.numeric(yfs[,3])
narea <- max(areaID)

# detection covariates

val <- detCovs.nb$v.nb[,3:32]
trail <- detCovs.nb$trail.nb[,3]
cam <- detCovs.nb$lr.nb[,3:32]


# site covariates

#nTC500 <- siteCovs.nb$nTC500
#rabbt.d <- siteCovs.nb$rbbt.ha
#rabbt.ts <- siteCovs.nb$rbbt.ts
#rod.d <- siteCovs.nb$rdt.ha
#rod.ts <- siteCovs.nb$rod.ts

occ.data <- list( yA=EHf, 
                  yB=EHmf,
                  nSites=nSites,
                  nocc=nocc
)

occ.inits <- function(){
  list(
    zA = rep(1,nSites), #apply( EHf, 1, max, na.rm=T ),
    zBA = rep(1,nSites), #apply( EHmf, 1, max, na.rm=T ), 
    zBa = rep(0,nSites), #apply( EHmf, 1, max, na.rm=T ) -
     # apply( EHf, 1, max, na.rm=T ) * apply( EHmf, 1, max, na.rm=T ), 
    b0.psiA = runif(1, -3, 3),
    b0.psiBA = runif(1, -3, 3),
    b0.psiBa = runif(1, -3, 3),
    b0.pA = runif(1, -3, 3),
    b0.pB = runif(1, -3, 3),
    b0.rBA = runif(1, -3, 3),
    b0.rBa = runif(1, -3, 3)
  )
}

occ.parm <- c( "b0.psiA",
               "b0.psiBA",
               "b0.psiBa",
               "b0.pA",
               "b0.pB",
               "b0.rBA",
               "b0.rBa",
               "phi",
               "ZZBa"
)

# set up for MCMC run
ni <- 5000
nt <- 1
nb <- 200
nc <- 3

# run the MCMC chain in JAGS
occ.fox.wildc <- jags( occ.data, 
                    occ.inits,
                    occ.parm,
                    "JAGSmodel_interaction.R",
                    n.chains=nc, 
                    n.iter=ni, 
                    n.burnin=nb,
                    n.thin=nt
)

mcmcplot(occ.fox.wildc)
# assess model fit
plot(occ.result$BUGSoutput$sims.list$fit, occ.result$BUGSoutput$sims.list$fit.new, main = "", xlab = 
       "Discrepancy for actual data set", ylab = "Discrepancy for perfect 
     data sets", las = 1, bty='n', pch=16, cex=0.8)
abline(0,1, lwd = 2)

mean(occ.result$BUGSoutput$sims.list$fit.new > occ.result$BUGSoutput$sims.list$fit) # acceptable values are those that fall between 0.1 and 0.9. Extremely low (<0.1) and high (>0.9) reveal structural failure of model fit
