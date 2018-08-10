#
# Occupancy analysis
#

library(mcmcplots)
library(R2jags)

setwd( "/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/1.rcode/occ/two.species/3.fox" )

# set detection histories and covariate data
source("/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/1.rcode/occ/single.species/data.R")

EHfox <- y$yvv[,c(5:34)]
EHwildcat <- y$yfs[,c(5:34)]
EHpmarten <- y$ymmt[,c(5:34)]
EHsmarten <- y$ymf[,c(5:34)]
EHmongoose <- y$yhi[,c(5:34)]
EHgenet <- y$ygg[,c(5:34)]

nocc <- ncol(EHfox)
nSites <- nrow(EHfox)

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
occ.data <- list( y=EHfox, # dominant species
                  ys1=EHwildcat,# subordinate species
                  ys2=EHpmarten,# subordinate species
                  ys3=EHsmarten,# subordinate species
                  ys4=EHmongoose,# subordinate species
                  ys5=EHgenet,# subordinate species
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
  list(z = as.integer(apply(EHfox,1,sum, na.rm=T)>0),
    zs1 = as.integer(apply(EHwildcat,1,sum, na.rm=T)>0),
    zs2 = as.integer(apply(EHpmarten,1,sum, na.rm=T)>0),
    zs3 = as.integer(apply(EHsmarten,1,sum, na.rm=T)>0),
    zs4 = as.integer(apply(EHmongoose,1,sum, na.rm=T)>0),
    zs5 = as.integer(apply(EHgenet,1,sum, na.rm=T)>0),
    #fox
    b0.psi = runif(1, -0.780184595,0.84234527),
    b1.psi = runif(1, -0.027902466,0.857708284),
    b2.psi = runif(1, -0.513166452,-0.098394659),
    b0.p = runif(1, -1.929433438,-1.704332405),
    b1.p = runif(1, 0.02525486,0.322775335),
    sigma = runif(1, 0, 0.1),
    #wildcat
    b0.psis1D = runif(1,-4, 4),
    b0.psis1 = runif(1, -4, 4),
    b0.ps1 = runif(1, -4, 4),
    #pmarten
    b0.psis2D = runif(1,-4, 4),
    b0.psis2 = runif(1, -4, 4),
    b0.ps2 = runif(1, -4, 4),
    b1.ps2 = runif(1, -4, 4),
    #smarten
    b0.psis3D = runif(1,-4, 4),
    b0.psis3 = runif(1, -4, 4),
    b0.ps3 = runif(1, -4, 4),
    #mongoose
    b0.psis4D = runif(1,-4, 4),
    b0.psis4 = runif(1, -4, 4),
    b0.ps4 = runif(1, -4, 4),
    #genet
    b0.psis5D = runif(1,-4, 4),
    b0.psis5 = runif(1, -4, 4),
    b0.ps5 = runif(1, -4, 4)
  )
}

occ.parm <- c( #fox
               "b0.psi",
               "b1.psi",
               "b2.psi",
               "b0.p",
               "b1.p",
               #wildcat
               "b0.psis1D",
               "b0.psis1",
               "b0.ps1",
               #pmarten
               "b0.psis2D",
               "b0.psis2",
               "b0.ps2",
               "b1.ps2",
               #smarten
               "b0.psis3D",
               "b0.psis3",
               "b0.ps3",
               #mongoose
               "b0.psis4D",
               "b0.psis4",
               "b0.ps4",
               #genet
               "b0.psis5D",
               "b0.psis5",
               "b0.ps5",
               #sif
               "phi1",
               "phi2",
               "phi3",
               "phi4",
               "phi5")

# set up for MCMC run
ni <- 50000
nt <- 10
nb <- 20000
nc <- 3

# run the MCMC chain in JAGS
multisps <- jags( occ.data, 
                occ.inits,
                occ.parm,
                "FOX.multisps.model_201807.R",
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
phis <- m[c("phi1","phi2","phi3","phi4","phi5"),]
phis$names <- c("fox/wildcat","fox/pmarten","fox/smarten","fox/mongoose","fox/genet")
x <- data.frame(multisps$BUGSoutput$sims.matrix)
phis$mode <- c(Mode(x$phi1),Mode(x$phi2),Mode(x$phi3),Mode(x$phi4),Mode(x$phi5))




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
x <- x[,c("phi1","phi2","phi3","phi4","phi5")]
x <- round(x,2)
x1 <- data.frame(c(x$phi1,x$phi2,x$phi3,x$phi4,x$phi5))
colnames(x1) <- "phi"
x1$names <- c(rep("fox-wildcat",nrow(x)),
              rep("fox-pmarten",nrow(x)),
              rep("fox-smarten",nrow(x)),
              rep("fox-mongoose",nrow(x)),
              rep("fox-genet",nrow(x)))

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

write.csv(phis,"phis.fox_201807.csv")




#calculate percentile for specific value in the posterior distribution
ecdf_fun <- function(x,perc) ecdf(x)(perc)
percentile <- ecdf(multisps$BUGSoutput$sims.matrix[,"phis1D"])
percentile(0.8)

#find specific percentile in the posterior distribution
quantile(multisps$BUGSoutput$sims.matrix[,"phis1D"],0.025)




