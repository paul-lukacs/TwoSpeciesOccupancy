#
# Occupancy analysis
#

library(mcmcplots)
library(R2jags)

setwd( "/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/1.rcode/occ/single.species/5.smarten" )


# set detection histories and covariate data
source("/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/1.rcode/occ/single.species/data.R")

EH <- y$ymf[,c(5:34)]
nocc <- ncol(EH)
nSites <- nrow(EH)
naive.mf <- sum(as.integer(apply(EH,1,sum, na.rm=T)>0))/nrow(EH)
sum(as.integer(apply(EH,1,sum, na.rm=T)))

# random effect
areaID <- as.numeric(siteCovs[,3])
narea <- max(areaID)

# detection covariates
trail <- detCovs$trail[,5]

# site covariates
EVI <- siteCovs$EVInb
rod.d <- siteCovs$rod.ha
rod.ts <- siteCovs$rod.ts
d.sett <- siteCovs$d.sett
season <- as.numeric(as.factor(siteCovs$season))-1

# eliminate data where control treatment (no attractant) was used
#for(i in 1:nrow(EH)){
#  for(j in 1:ncol(EH)){
#    if(ctr[i,j]==1){
#      EH[i,j]<-NA
#    }
#  }
#}

sum(as.integer(apply(EH,1,sum, na.rm=T)))

occ.data <- list( y=EH, 
                  nSites=nSites,
                  nocc=nocc,
                  EVI=EVI,
                  rod.ts=rod.ts,
                  season=season,
                  trail=trail,
                  areaID=areaID,
                  narea=narea
)

occ.inits <- function(){
  list(
    z = as.integer(apply(EH,1,sum, na.rm=T)>0),
    b0.psi = runif(1, -3, 3),
    b1.psi = runif(1, -3, 3),
    b2.psi = runif(1, -3, 3),
    b3.psi = runif(1, -3, 3),
    b4.psi = runif(1, -3, 3),
    b0.p = runif(1, -3, 3),
    b1.p = runif(1, -3, 3),
    b2.p = runif(1, -3, 3),
    sigma = runif(1, 0, 0.1)
  )
}


occ.parm <- c( "b0.psi",
               "b1.psi",
               "b2.psi",
               "b3.psi",
               "b4.psi",
               "b0.p",
               "b1.p",
               "b2.p",
               "mean.psi",
               "mean.p",
               "totOcc",
               "fit",
               "fit.new",
               "sigma"
)
# set up for MCMC run
ni <- 50000
nt <- 10
nb <- 20000
nc <- 3


# run the MCMC chain in JAGS
smarten <- jags( occ.data, 
                 occ.inits,
                 occ.parm,
                 "smarten.model.20180610.R",
                 n.chains=nc, 
                 n.iter=ni, 
                 n.burnin=nb,
                 n.thin=nt
)

# verify model results
smarten
mcmcplot(smarten)


# assess model fit
plot(smarten$BUGSoutput$sims.list$fit, smarten$BUGSoutput$sims.list$fit.new, main = "", xlab = 
       "Discrepancy for actual data set", ylab = "Discrepancy for perfect data sets", las = 1, bty='n', pch=16, cex=0.8, 
     xlim=c(min(c(smarten$BUGSoutput$sims.list$fit, smarten$BUGSoutput$sims.list$fit.new)),max(c(smarten$BUGSoutput$sims.list$fit, smarten$BUGSoutput$sims.list$fit.new))),
     ylim=c(min(c(smarten$BUGSoutput$sims.list$fit, smarten$BUGSoutput$sims.list$fit.new)),max(c(smarten$BUGSoutput$sims.list$fit, smarten$BUGSoutput$sims.list$fit.new))))
abline(0,1, lwd = 2, col='red')

mean(smarten$BUGSoutput$sims.list$fit.new > smarten$BUGSoutput$sims.list$fit) # acceptable values are those that fall between 0.1 and 0.9. Extremely low (<0.1) and high (>0.9) reveal structural failure of model fit

#function to calculate the inverse logit
expit <- function(x){
  exp(x)/(1+exp(x))
}

#function to obtain the mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

## function to obtain the percentile for specific value 
ecdf_fun <- function(x,perc) ecdf(x)(perc)

# retain variables relevant at the 85% credible interval
m <- data.frame(smarten$BUGSoutput$summary)
for(i in 1:nrow(m)){
  m$zero.perctl[i] <- ecdf_fun(smarten$BUGSoutput$sims.matrix[,i],0)
  if(m$zero.perctl[i]<0.025|m$zero.perctl[i]>0.975){
    m$retain[i]="yes"}else{
      m$retain[i]="no"
    }
}
m <- m[c(1,3,5,2,4,6:8,12:13),]  
m$covariate <- c("p.intercept",
                 "p.trail",
                 "p.season",
                 "psi.intercept",
                 "EVI",
                 "EVI^2",
                 "rodent.ts",
                 "season",
                 "mean.p",
                 "mean.psi"
)
write.csv(m,"smarten_20180610.csv")

#generate data frame with estimates for each psi beta
x <- data.frame(smarten$BUGSoutput$sims.matrix)[,c("b1.psi","b2.psi","b3.psi","b4.psi")]
x <- round(x,2)
x1 <- data.frame(c(x$b1.psi,x$b2.psi,x$b3.psi,x$b4.psi))
colnames(x1) <- "Estimate"
x1$names <- c(rep("EVI",nrow(x)),
              rep("EVI^2",nrow(x)),
              rep("rodent ts",nrow(x)),
              rep("season",nrow(x)))

#install.packages("devtools")
#library(devtools)
#install_github("easyGgplot2", "kassambara")
library(easyGgplot2)

p <- ggplot(x1, aes(x=Estimate, fill=names)) +
  geom_density(alpha=0.4) +
  scale_fill_brewer(palette = "Paired") +
  theme_minimal() +
  scale_x_continuous(limits = c(-2.5, 2.5)) + 
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
  ggtitle("Stone marten") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ names, ncol = 3, nrow = 3) +
  theme(legend.position="none")

print(p)
pdf("smarten_20180610.pdf")
print(p)
dev.off()

# Get codes for colours in a palette
library("RColorBrewer")
brewer.pal(12, "Paired")
# to display that palette:
display.brewer.pal(12, "Paired")

