#
# Occupancy analysis
#

library(mcmcplots)
library(R2jags)

setwd( "/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/1.rcode/occ/two.species/6.pmarten" )

# set detection histories and covariate data
source("/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/1.rcode/occ/single.species/data.R")

EHpmarten <- y$ymmt[,c(5:34)]
EHsmarten <- y$ymf[,c(5:34)]
EHmongoose <- y$yhi[,c(5:34)]
EHgenet <- y$ygg[,c(5:34)]

nocc <- ncol(EHpmarten)
nSites <- nrow(EHpmarten)

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
occ.data <- list( y=EHpmarten, # dominant species
                  ys1=EHsmarten,# subordinate species
                  ys2=EHmongoose,# subordinate species
                  ys3=EHgenet,# subordinate species
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
  list(z = as.integer(apply(EHpmarten,1,sum, na.rm=T)>0),
    zs1 = as.integer(apply(EHsmarten,1,sum, na.rm=T)>0),
    zs2 = as.integer(apply(EHmongoose,1,sum, na.rm=T)>0),
    zs3 = as.integer(apply(EHgenet,1,sum, na.rm=T)>0),
    #pmarten
    b0.psi = runif(1, -2.543709287, -0.516059995),
    b0.p = runif(1, -3.575116711, -2.782373413),
    b1.p = runif(1, 0.05567195, 1.759690078),
    sigma = runif(1, 0, 0.1),
    #smarten
    b0.psis1D = runif(1,-4, 4),
    b0.psis1 = runif(1, -4, 4),
    b0.ps1 = runif(1, -4, 4),
    #mongoose
    b0.psis2D = runif(1,-4, 4),
    b0.psis2 = runif(1, -4, 4),
    b0.ps2 = runif(1, -4, 4),
    #genet
    b0.psis3D = runif(1,-4, 4),
    b0.psis3 = runif(1, -4, 4),
    b0.ps3 = runif(1, -4, 4),
    b1.ps3 = runif(1, -4, 4)
  )
}

occ.parm <- c( #pmarten
               "b0.psi",
               "b0.p",
               "b1.p",
               "sigma",
               #smarten
               "b0.psis1D",
               "b0.psis1",
               "b0.ps1",
               #mongoose
               "b0.psis2D",
               "b0.psis2",
               "b0.ps2",
               #genet
               "b0.psis3D",
               "b0.psis3",
               "b0.ps3",
               "b1.ps3",
               #sif
               "phi1",
               "phi2",
               "phi3")

# set up for MCMC run
ni <- 50000
nt <- 10
nb <- 20000
nc <- 3

# run the MCMC chain in JAGS
multisps <- jags( occ.data, 
                occ.inits,
                occ.parm,
                "PMARTEN.multisps.model_201807.R",
                n.chains=nc, 
                n.iter=ni, 
                n.burnin=nb,
                n.thin=nt
)

# inspecting model results
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
phis <- m[c("phi1","phi2","phi3"),]
phis$names <- c("pmarten/smarten","pmarten/mongoose","pmarten/genet")
x <- data.frame(multisps$BUGSoutput$sims.matrix)
phis$mode <- c(Mode(x$phi1),Mode(x$phi2),Mode(x$phi3))
write.csv(phis,"phis.pmarten_201807.csv")

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
x <- x[,c("phi1","phi2","phi3")]
x <- round(x,2)
x1 <- data.frame(c(x$phi1,x$phi2,x$phi3))
colnames(x1) <- "phi"
x1$names <- c(rep("pmarten-smarten",nrow(x)),
              rep("pmarten-mongoose",nrow(x)),
              rep("pmarten-genet",nrow(x)))

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






#calculate percentile for specific value in the posterior distribution
ecdf_fun <- function(x,perc) ecdf(x)(perc)
percentile <- ecdf(multisps$BUGSoutput$sims.matrix[,"phis1D"])
percentile(0.8)

#find specific percentile in the posterior distribution
quantile(multisps$BUGSoutput$sims.matrix[,"phis1D"],0.025)




