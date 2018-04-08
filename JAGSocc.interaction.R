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
  setwd( ""/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/1.rcode/occ/two.species"" )
}


# set detection histories and covariate data

source("/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/1.rcode/occ/single.species/data.nb.R")

# data
yD <- EHvv.nb # dominant species (red fox)
ys1 <- EHgg.nb # subordinate species (common genet)
ys2 <- EHmf.nb # subordinate species (stone marten)

# random effect

areaID.nb <- as.numeric(y.nb[,3])
narea <- max(areaID)

# detection covariates

#trail <- detCovs.nb$trail.nb[,3]

# site covariates

#nTC500 <- siteCovs.nb$nTC500
#rabbt.d <- siteCovs.nb$rbbt.ha
#rabbt.ts <- siteCovs.nb$rbbt.ts
#rod.d <- siteCovs.nb$rdt.ha
#rod.ts <- siteCovs.nb$rod.ts

occ.data <- list( yD=yD, 
                  ys1=ys1,
                  ys2=ys2,
                  nSites=nSites,
                  nocc=nocc,
                  areaID=areaID.nb,
                  narea=narea
)

occ.inits <- function(){
  list(
    zD = rep(1,nSites.nb), 
    zs1 = rep(1,nSites.nb), 
    zs2 = rep(1,nSites.nb), 
    b0.psiD = runif(1, -3, 3),
    b0.psis1D = runif(1, -3, 3),
    b0.psis1 = runif(1, -3, 3),
    b0.psis2D = runif(1, -3, 3),
    b0.psis2 = runif(1, -3, 3),
    b0.pD = runif(1, -3, 3),
    b0.ps1D = runif(1, -3, 3),
    b0.ps1 = runif(1, -3, 3),
    b0.ps2D = runif(1, -3, 3),
    b0.ps2 = runif(1, -3, 3),
    sigma.psi = runif(1, 0, 0.1)
  )
}

occ.parm <- c( "b0.psiD",
               "b0.psis1D",
               "b0.psis1",
               "b0.psis2D",
               "b0.psis2",
               "b0.pD",
               "b0.ps1D",
               "b0.ps1",
               "b0.ps2D",
               "b0.ps2",
               "sigma.psi",
               "phis1D",
               "phis2D"
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
