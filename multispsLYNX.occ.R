#
# Occupancy analysis
#

library(mcmcplots)
library(R2jags)

setwd( "/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/1.rcode/occ/two.species" )

# set detection histories and covariate data
source("/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/1.rcode/occ/single.species/data.R")

EHlynx <- y$ylp[,c(5:68)]
EHbadger <- y$ymm[,c(5:68)]
EHfox <- y$yvv[,c(5:68)]
EHwildcat <- y$yfs[,c(5:68)]
EHsmarten <- y$ymf[,c(5:68)]
EHgenet <- y$ygg[,c(5:68)]

nocc <- ncol(EHlynx)
nSites <- nrow(EHlynx)

#naive occupancy estimates
naive.est<- data.frame(c(sum(as.integer(apply(EHlynx,1,sum, na.rm=T)>0))/nrow(EHlynx),
                         sum(as.integer(apply(EHbadger,1,sum, na.rm=T)>0))/nrow(EHbadger),
                         sum(as.integer(apply(EHfox,1,sum, na.rm=T)>0))/nrow(EHfox),
                         sum(as.integer(apply(EHwildcat,1,sum, na.rm=T)>0))/nrow(EHwildcat),
                         sum(as.integer(apply(EHsmarten,1,sum, na.rm=T)>0))/nrow(EHsmarten),
                         sum(as.integer(apply(EHgenet,1,sum, na.rm=T)>0))/nrow(EHgenet)))
colnames(naive.est) <- 'naive.psi'
naive.est$species <- c('lynx','badger','fox','wildcat','smarten','genet')

# random effect
areaID <- as.numeric(siteCovs[,3])
narea <- max(areaID)

# detection covariates
trail <- detCovs$trail[,5]
ctr <- detCovs$ctr[,c(5:68)]

# site covariates
EVI <- siteCovs$EVInb
rbbt.ts <- siteCovs$rbbt.ts
rod.ts <- siteCovs$rod.ts
d.sett <- siteCovs$d.sett
season <- as.numeric(as.factor(siteCovs$season))-1


# data
occ.data <- list( yD=EHlynx, # dominant species
                  ys1=EHbadger,# subordinate species 1
                  ys2=EHfox,# subordinate species 2
                  ys3=EHwildcat,# subordinate species 3
                  ys4=EHsmarten,# subordinate species 4
                  ys5=EHgenet,# subordinate species 5
                  nSites=nSites,
                  EVI=EVI,
                  rbbt.ts=rbbt.ts,
                  rod.ts=rod.ts,
                  d.sett=d.sett,
                  ctr=ctr,
                  nocc=nocc,
                  areaID=areaID,
                  narea=narea
)

# initial values
occ.inits <- function(){
  list(
    zD = as.integer(apply(EHlynx,1,sum, na.rm=T)>0), 
    zs1 = as.integer(apply(EHbadger,1,sum, na.rm=T)>0), 
    zs2 = as.integer(apply(EHfox,1,sum, na.rm=T)>0), 
    zs3 = as.integer(apply(EHwildcat,1,sum, na.rm=T)>0),
    zs4 = as.integer(apply(EHsmarten,1,sum, na.rm=T)>0),
    zs5 = as.integer(apply(EHgenet,1,sum, na.rm=T)>0),
    #occupancy
    b0.psiD = runif(1, -3, 3),
    b0.psis1D = runif(1, -3, 3),
    b0.psis1 = runif(1, -3, 3),
    b1.psis1 = runif(1, -3, 3),
    b0.psis2D = runif(1, -3, 3),
    b0.psis2 = runif(1, -3, 3),
    b1.psis2 = runif(1, -3, 3),
    b2.psis2 = runif(1, -3, 3),
    b0.psis3D = runif(1, -3, 3),
    b0.psis3 = runif(1, -3, 3),
    b1.psis3 = runif(1, -3, 3),
    b2.psis3 = runif(1, -3, 3),
    b0.psis4D = runif(1, -3, 3),
    b0.psis4 = runif(1, -3, 3),
    b1.psis4 = runif(1, -3, 3),
    b0.psis5D = runif(1, -3, 3),
    b0.psis5 = runif(1, -3, 3),
    #detection
    b0.pD = runif(1, -3, 3),
    b0.ps1D = runif(1, -3, 3),
    b0.ps1 = runif(1, -3, 3),
    b0.ps2D = runif(1, -3, 3),
    b0.ps2 = runif(1, -3, 3),
    b1.ps2 = runif(1, -3, 3),
    b0.ps3D = runif(1, -3, 3),
    b0.ps3 = runif(1, -3, 3),
    b0.ps4D = runif(1, -3, 3),
    b0.ps4 = runif(1, -3, 3),
    b1.ps4 = runif(1, -3, 3),
    b0.ps5D = runif(1, -3, 3),
    b0.ps5 = runif(1, -3, 3),
    #random effects
    sigma.psis1 = runif(1, 0, 0.1),
    sigma.psis2 = runif(1, 0, 0.1),
    sigma.psis3 = runif(1, 0, 0.1),
    sigma.psis4 = runif(1, 0, 0.1),
    sigma.psis5 = runif(1, 0, 0.1)
  )
}

occ.parm <- c( "b0.psiD",
               "b0.psis1D",
               "b0.psis1",
               "b1.psis1",
               "b0.psis2D",
               "b0.psis2",
               "b1.psis2",
               "b2.psis2",
               "b0.psis3D",
               "b0.psis3",
               "b1.psis3",
               "b2.psis3",
               "b0.psis4D",
               "b0.psis4",
               "b2.psis4",
               "b0.psis5D",
               "b0.psis5",
               "b0.pD",
               "b0.ps1D",
               "b0.ps1",
               "b0.ps2D",
               "b0.ps2",
               "b1.ps2",
               "b0.ps3D",
               "b0.ps3",
               "b0.ps4D",
               "b0.ps4",
               "b1.ps4",
               "b0.ps5D",
               "b0.ps5",
               "sigma.psis1",
               "sigma.psis2",
               "sigma.psis3",
               "sigma.psis4",
               "sigma.psis5",
               "phis1D",
               "phis2D",
               "phis3D",
               "phis4D",
               "phis5D"
)

# set up for MCMC run
ni <- 50000 #500000
nt <- 10     #3
nb <- 20000 #20000
nc <- 3

# run the MCMC chain in JAGS
multisps <- jags( occ.data, 
                occ.inits,
                occ.parm,
                "multispsLYNX.model.R",
                n.chains=nc, 
                n.iter=ni, 
                n.burnin=nb,
                n.thin=nt
)
mcmcplot(multisps)

#function to obtain the mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
#function to calculate the inverse logit
expit <- function(x){
  exp(x)/(1+exp(x))
}

m <- data.frame(multisps$BUGSoutput$summary)
phis <- m[c("phis1D","phis2D","phis3D","phis4D","phis5D"),]
phis$names <- c("lynx/badger","lynx/fox","lynx/wildcat","lynx/smarten","lynx/genet")

library(easyGgplot2)
p <- ggplot(phis, aes(x=names, y=mean)) + 
  geom_point()+
  geom_errorbar(aes(ymin=X2.5., ymax=X97.5.), width=.1,
                position=position_dodge(0.05)) +
  theme_minimal() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey")  +
  xlab("Species pair") + 
  ylab("Species Interaction Factor") 
    #ggtitle("") +
    #theme(plot.title = element_text(hjust = 0.5))
  #scale_y_continuous(limits = c(0,10))
print(p)  

head(multisps$BUGSoutput$sims.matrix[,"phis1D"])

#calculate percentile for specific value in the posterior distribution
ecdf_fun <- function(x,perc) ecdf(x)(perc)
percentile <- ecdf(multisps$BUGSoutput$sims.matrix[,"phis1D"])
percentile(0.8)

#find specific percentile in the posterior distribution
quantile(multisps$BUGSoutput$sims.matrix[,"phis1D"],0.025)




