#
# Occupancy analysis
#

library(mcmcplots)
library(R2jags)

setwd( "/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/1.rcode/occ/two.species/2.badger" )

# set detection histories and covariate data
source("/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/1.rcode/occ/single.species/data.R")

EHbadger <- y$ymm[,c(5:34)]
EHfox <- y$yvv[,c(5:34)]
EHwildcat <- y$yfs[,c(5:34)]
EHpmarten <- y$ymmt[,c(5:34)]
EHmongoose <- y$yhi[,c(5:34)]
EHsmarten <- y$ymf[,c(5:34)]
EHgenet <- y$ygg[,c(5:34)]

nocc <- ncol(EHbadger)
nSites <- nrow(EHbadger)

# random effect
areaID <- as.numeric(siteCovs[,3])
narea <- max(areaID)

# detection covariates
trail <- detCovs$trail[,5]

# site covariates
EVI <- siteCovs$EVInb
rbbt.ts <- siteCovs$rbbt.ts
d.sett <- siteCovs$d.sett
season <- as.numeric(as.factor(siteCovs$season))-1

# data
occ.data <- list( y=EHbadger, # dominant species
                  ys1=EHfox,# subordinate species
                  ys2=EHwildcat,# subordinate species
                  ys3=EHpmarten,# subordinate species
                  ys4=EHsmarten,# subordinate species
                  ys5=EHmongoose,# subordinate species
                  ys6=EHgenet,# subordinate species
                  nSites=nSites,
                  EVI=EVI,
                  d.sett=d.sett,
                  rbbt.ts=rbbt.ts,
                  season=season,
                  trail=trail,
                  nocc=nocc,
                  areaID=areaID,
                  narea=narea
)

# initial values
occ.inits <- function(){
  list(z = as.integer(apply(EHbadger,1,sum, na.rm=T)>0),
    zs1 = as.integer(apply(EHfox,1,sum, na.rm=T)>0),
    zs2 = as.integer(apply(EHwildcat,1,sum, na.rm=T)>0),
    zs3 = as.integer(apply(EHpmarten,1,sum, na.rm=T)>0),
    zs4 = as.integer(apply(EHsmarten,1,sum, na.rm=T)>0),
    zs5 = as.integer(apply(EHmongoose,1,sum, na.rm=T)>0),
    zs6 = as.integer(apply(EHgenet,1,sum, na.rm=T)>0),
    #badger
    b0.psi = runif(1, -3.418422,-2.34292),
    b1.psi = runif(1, -2.369018,-0.2423185),
    b2.psi = runif(1, 0.3334643,0.8753147),
    b3.psi = runif(1, 0.2055945,1.996876),
    b0.p = runif(1, -2.868521,-2.454211),
    b1.p = runif(1, -2.110455,-1.026539),
    sigma = runif(1, 0, 0.1),
    #fox
    b0.psis1D = runif(1,-4, 4),
    b1.psis1D = runif(1,-4, 4),
    b2.psis1D = runif(1,-4, 4),
    b0.psis1 = runif(1, -4, 4),
    b1.psis1 = runif(1, -4, 4),
    b2.psis1 = runif(1, -4, 4),
    b0.ps1 = runif(1, -4, 4),
    b1.ps1 = runif(1, -4, 4),
    #wildcat
    b0.psis2D = runif(1,-4, 4),
    b0.psis2 = runif(1, -4, 4),
    b0.ps2 = runif(1, -4, 4),
    #pmarten
    b0.psis3D = runif(1,-4, 4),
    b0.psis3 = runif(1, -4, 4),
    b0.ps3 = runif(1, -4, 4),
    b1.ps3 = runif(1, -4, 4),
    #smarten
    b0.psis4D = runif(1,-4, 4),
    b0.psis4 = runif(1, -4, 4),
    b0.ps4 = runif(1, -4, 4),
    #mongoose
    b0.psis5D = runif(1,-4, 4),
    b0.psis5 = runif(1, -4, 4),
    b0.ps5 = runif(1, -4, 4),
    #genet
    b0.psis6D = runif(1,-4, 4),
    b0.psis6 = runif(1, -4, 4),
    b0.ps6 = runif(1, -4, 4),
    b1.ps6 = runif(1, -4, 4)
  )
}

occ.parm <- c( #badger
               "b0.psi",
               "b1.psi",
               "b2.psi",
               "b3.psi",
               "b0.p",
               "b1.p",
               "sigma",
               #fox
               "b0.psis1D",
               "b1.psis1D",
               "b2.psis1D",
               "b0.psis1",
               "b1.psis1",
               "b2.psis1",
               "b0.ps1",
               "b1.ps1",
               #wildcat
               "b0.psis2D",
               "b0.psis2",
               "b0.ps2",
               #pmarten
               "b0.psis3D",
               "b0.psis3",
               "b0.ps3",
               "b1.ps3",
               #smarten
               "b0.psis4D",
               "b0.psis4",
               "b0.ps4",
               #mongoose
               "b0.psis5D",
               "b0.psis5",
               "b0.ps5",
               #genet
               "b0.psis6D",
               "b0.psis6",
               "b0.ps6",
               "b1.ps6",
               #sif
               "phi1",
               "phi2",
               "phi3",
               "phi4",
               "phi5",
               "phi6")

# set up for MCMC run
ni <- 50000
nt <- 10
nb <- 20000
nc <- 3

# run the MCMC chain in JAGS
multisps <- jags( occ.data, 
                occ.inits,
                occ.parm,
                "BADGER.multisps.model_201807.R",
                n.chains=nc, 
                n.iter=ni, 
                n.burnin=nb,
                n.thin=nt
)

# inpecting model results
multisps
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
phis <- m[c("phi1","phi2","phi3","phi4","phi5","phi6"),]
phis$names <- c("badger/fox","badger/wildcat","badger/pmarten","badger/smarten","badger/mongoose","badger/genet")
x <- data.frame(multisps$BUGSoutput$sims.matrix)
phis$mode <- c(Mode(x$phi1),Mode(x$phi2),Mode(x$phi3),Mode(x$phi4),Mode(x$phi5),Mode(x$phi6))


library(easyGgplot2)
p <- ggplot(phis, aes(x=names, y=mean)) + 
  geom_point()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1,
                position=position_dodge(0.05)) +
  theme_minimal() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey")  +
  xlab("Species pair") + 
  ylab("Species Interaction Factor") 
    #ggtitle("") +
    #theme(plot.title = element_text(hjust = 0.5))
  #scale_y_continuous(limits = c(0,10))
print(p)  

#generate data frame with estimates for each phi
x <- x[,c("phi1","phi2","phi3","phi4","phi5","phi6")]
x <- round(x,2)
x1 <- data.frame(c(x$phi1,x$phi2,x$phi3,x$phi4,x$phi5,x$phi6))
colnames(x1) <- "phi"
x1$names <- c(rep("badger-fox",nrow(x)),
              rep("badger-wildcat",nrow(x)),
              rep("badger-pmarten",nrow(x)),
              rep("badger-smarten",nrow(x)),
              rep("badger-mongoose",nrow(x)),
              rep("badger-genet",nrow(x)))

q <- ggplot(x1, aes(x=phi, fill=names)) +
  geom_density(alpha=0.4) +
  scale_fill_brewer(palette = "Paired") +
  theme_minimal() +
  scale_x_continuous(limits = c(-1, 4)) + 
  geom_vline(xintercept = 1, linetype = "dashed", colour = "black") +
  ggtitle("Species Interaction Factor") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ names, nrow = 3, ncol = 3) +
  theme(legend.position="none")

print(q)

write.csv(phis,"phis.badger_201807.csv")




#calculate percentile for specific value in the posterior distribution
ecdf_fun <- function(x,perc) ecdf(x)(perc)
percentile <- ecdf(multisps$BUGSoutput$sims.matrix[,"phis1D"])
percentile(0.8)

#find specific percentile in the posterior distribution
quantile(multisps$BUGSoutput$sims.matrix[,"phis1D"],0.025)




