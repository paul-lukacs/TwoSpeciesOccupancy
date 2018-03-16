
########################
# SET OBSERVATION DATA #
########################

yfs <- read.csv("/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/1.data/occ.data/yfs.nb.csv",sep=",",header = T)
yfs <- yfs[,-1]

yvv <- read.csv("/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/1.data/occ.data/yvv.nb.csv",sep=",",header = T)
yvv <- yvv[,-1]

ymf <- read.csv("/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/1.data/occ.data/ymf.nb.csv",sep=",",header = T)
ymf <- ymf[,-1]

###########################
# PREPARE SITE COVARIATES #
###########################
sCovs.nb <- read.csv("/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/1.data/occ.data/temp/siteCovs.nb.csv",sep=",",header = T)
sCovs.nb <- sCovs.nb[,-1]
#sCovs.nb <- nb[,c(2,36:66)]

siteCovs.nb<-NULL
siteCovs.nb$station <- sCovs.nb$station
siteCovs.nb<-data.frame(siteCovs.nb)

## TRANSFORM STANDARDIZE RAW COVARIATES (CONTINUOUS TO Z-SCORES; PROPORTIONAL TO ARCSINE)  ##

# LATITUDE
attach(sCovs.nb)
#latitude
siteCovs.nb$lat <- (utm_y-mean(utm_y))/sd(utm_y)
#water availability
siteCovs.nb$d.wb <- (dst.waterbodies-mean(dst.waterbodies))/sd(dst.waterbodies)
siteCovs.nb$d.rivers <- (dst_rivers-mean(dst_rivers))/sd(dst_rivers)
siteCovs.nb$d.streams <- (dst_streams-mean(dst_streams))/sd(dst_streams)
siteCovs.nb$d.wtr <- (dst.water-mean(dst.water))/sd(dst.water)
#human disturbance
siteCovs.nb$d.sett <- (dst_settlement-mean(dst_settlement))/sd(dst_settlement)
siteCovs.nb$d.road <- (dst_roads-mean(dst_roads))/sd(dst_roads)
siteCovs.nb$d.hmn <- (dst.human-mean(dst.human))/sd(dst.human)
siteCovs.nb$rds1000 <- (roads1000m-mean(roads1000m))/sd(roads1000m)
siteCovs.nb$rds500 <- (roads500m-mean(roads500m))/sd(roads500m)
siteCovs.nb$rds200 <- (roads200m-mean(roads200m))/sd(roads200m)
#prey
siteCovs.nb$rbbt.ts <- (rabbit.ts-mean(rabbit.ts, na.rm = T))/sd(rabbit.ts, na.rm = T)
siteCovs.nb$rod.ts <- (rodent.ts-mean(rodent.ts, na.rm = T))/sd(rodent.ts, na.rm = T)
siteCovs.nb$rbbt.ha <- (rbbt.ha-mean(rbbt.ha, na.rm = T))/sd(rbbt.ha, na.rm = T)
siteCovs.nb$rdt.ha <- (rdt.ha-mean(rdt.ha, na.rm = T))/sd(rdt.ha, na.rm = T)
siteCovs.nb$As.ha <- (As.ha-mean(As.ha, na.rm = T))/sd(As.ha, na.rm = T)
#lynx
siteCovs.nb$lynx.ts <- (lynx.ts-mean(lynx.ts, na.rm = T))/sd(lynx.ts, na.rm = T)
#tree cover
siteCovs.nb$TC1000 <- (TC1000-mean(TC1000, na.rm = T))/sd(TC1000, na.rm = T)
siteCovs.nb$TC500 <- (TC500-mean(TC500, na.rm = T))/sd(TC500, na.rm = T)
siteCovs.nb$TC200 <- (TC200-mean(TC200, na.rm = T))/sd(TC200, na.rm = T)
#non-tree cover
siteCovs.nb$nTC1000 <- (nTC1000-mean(nTC1000, na.rm = T))/sd(nTC1000, na.rm = T)
siteCovs.nb$nTC500 <- (nTC500-mean(nTC500, na.rm = T))/sd(nTC500, na.rm = T)
siteCovs.nb$nTC200 <- (nTC200-mean(nTC200, na.rm = T))/sd(nTC200, na.rm = T)
#open areas
siteCovs.nb$nV1000 <- (nV1000-mean(nV1000, na.rm = T))/sd(nV1000, na.rm = T)
siteCovs.nb$nV500 <- (nV500-mean(nV500, na.rm = T))/sd(nV500, na.rm = T)
siteCovs.nb$nV200 <- (nV200-mean(nV200, na.rm = T))/sd(nV200, na.rm = T)
#EVI
siteCovs.nb$EVI1000 <- (EVI1000-mean(EVI1000, na.rm = T))/sd(EVI1000, na.rm = T)
siteCovs.nb$EVI500 <- (EVI500-mean(EVI500, na.rm = T))/sd(EVI500, na.rm = T)
siteCovs.nb$EVI200 <- (EVI200-mean(EVI200, na.rm = T))/sd(EVI200, na.rm = T)
#sdEVI
siteCovs.nb$sdEVI1000 <- (sdEVI1000-mean(sdEVI1000, na.rm = T))/sd(sdEVI1000, na.rm = T)
siteCovs.nb$sdEVI500 <- (sdEVI500-mean(sdEVI500, na.rm = T))/sd(sdEVI500, na.rm = T)
siteCovs.nb$sdEVI200 <- (sdEVI200-mean(sdEVI200, na.rm = T))/sd(sdEVI200, na.rm = T)

rm(sCovs.nb)

###########################################
# LOAD AND PREPARE OBSERVATION COVARIATES #
###########################################

ynb <- read.csv('/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/1.data/occ.data/temp/siteCovs.nb.csv', sep = ',', header = T)
ynb <- ynb[,-1]

yb <- read.csv('/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/1.data/occ.data/temp/siteCovs.b.csv', sep = ',', header = T)
yb <- yb[,-1]

ct.nb <- data.frame(ynb[,c(2,1)])
colnames(ct.nb) <- c("station","study.area")

ct.b <- data.frame(yb[,c(2,1)])
colnames(ct.b) <- c("station","study.area")

write.csv(ct.nb, '/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/5.detection covariates/ct.nb.csv')

write.csv(ct.b, '/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/5.detection covariates/ct.b.csv')

require(gdata)

# upload detection covariates
### lynx urine
lu.nb<-read.xls(xls = '/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/5.detection covariates/detCovs.nb.xls',sheet = 'lynx.urine')

lu.nb <- merge(ct.nb, lu.nb, all.x = T, by = 'station', all.y = F, incomparables = NA)
lu.nb <- lu.nb[,c(1:2,4:34)]
colnames(lu.nb)[2] <- 'study.area'

lu.b<-read.xls(xls = '/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/5.detection covariates/detCovs.b.xls',sheet = 'lynx.urine')

lu.b <- merge(ct.b, lu.b, all.x = T, by = 'station', all.y = F, incomparables = NA)
lu.b <- lu.b[,c(1:2,4:33)]
colnames(lu.b)[2] <- 'study.area'

### valerian
v.nb<-read.xls(xls = '/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/5.detection covariates/detCovs.nb.xls',sheet = 'valerian')

v.nb <- merge(ct.nb, v.nb, all.x = T, by = 'station', all.y = F, incomparables = NA)
v.nb <- v.nb[,c(1:2,4:33)]
colnames(v.nb)[2] <- 'study.area'

v.b<-read.xls(xls = '/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/5.detection covariates/detCovs.b.xls',sheet = 'valerian')

v.b <- merge(ct.b, v.b, all.x = T, by = 'station', all.y = F, incomparables = NA)
v.b <- v.b[,c(1:2,4:33)]
colnames(v.b)[2] <- 'study.area'

### no.attractant
ctr.nb<-read.xls(xls = '/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/5.detection covariates/detCovs.nb.xls',sheet = 'no.attract')

ctr.nb <- merge(ct.nb, ctr.nb, all.x = T, by = 'station', all.y = F, incomparables = NA)
ctr.nb <- ctr.nb[,c(1:2,4:33)]
colnames(ctr.nb)[2] <- 'study.area'

ctr.b<-read.xls(xls = '/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/5.detection covariates/detCovs.b.xls',sheet = 'no.attract')

ctr.b <- merge(ct.b, ctr.b, all.x = T, by = 'station', all.y = F, incomparables = NA)
ctr.b <- ctr.b[,c(1:2,4:33)]
colnames(ctr.b)[2] <- 'study.area'

### on trail
trail.nb<-read.xls(xls = '/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/5.detection covariates/detCovs.nb.xls',sheet = 'on.trail')

trail.nb <- merge(ct.nb, trail.nb, all.x = T, by = 'station', all.y = F, incomparables = NA)
trail.nb <- trail.nb[,c(1:2,4)]
colnames(trail.nb)[2] <- 'study.area'

rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

r<-data.frame(rep.col(trail.nb$on.trail,30))
trail.nb <- cbind(trail.nb,r)
colnames(trail.nb) <- colnames((lu.nb))

trail.b<-read.xls(xls = '/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/5.detection covariates/detCovs.b.xls',sheet = 'on.trail')

trail.b <- merge(ct.b, trail.b, all.x = T, by = 'station', all.y = F, incomparables = NA)
trail.b <- trail.b[,c(1:2,4)]
colnames(trail.b)[2] <- 'study.area'

r<-data.frame(rep.col(trail.b$on.trail,29))
trail.b <- cbind(trail.b,r)
colnames(trail.b) <- colnames((lu.b))

### camera model: leafriver
lr.nb<-read.xls(xls = '/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/5.detection covariates/detCovs.nb.xls',sheet = 'leafriver')

lr.nb <- merge(ct.nb, lr.nb, all.x = T, by = 'station', all.y = F, incomparables = NA)
lr.nb <- lr.nb[,c(1:2,4:33)]
colnames(lr.nb)[2] <- 'study.area'

lr.b<-read.xls(xls = '/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/5.detection covariates/detCovs.b.xls',sheet = 'leafriver')

lr.b <- merge(ct.b, lr.b, all.x = T, by = 'station', all.y = F, incomparables = NA)
lr.b <- lr.b[,c(1:2,4:33)]
colnames(lr.b)[2] <- 'study.area'

### camera model: HCO
hco.nb<-read.xls(xls = '/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/5.detection covariates/detCovs.nb.xls',sheet = 'HCO')

hco.nb <- merge(ct.nb, hco.nb, all.x = T, by = 'station', all.y = F, incomparables = NA)
hco.nb <- hco.nb[,c(1:2,4:33)]
colnames(hco.nb)[2] <- 'study.area'

hco.b<-read.xls(xls = '/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/5.detection covariates/detCovs.b.xls',sheet = 'HCO')

hco.b <- merge(ct.b, hco.b, all.x = T, by = 'station', all.y = F, incomparables = NA)
hco.b <- hco.b[,c(1:2,4:33)]
colnames(hco.b)[2] <- 'study.area'

detCovs.nb <- list(lu.nb=lu.nb,v.nb=v.nb,ctr.nb=ctr.nb,lr.nb=lr.nb,hco.nb=hco.nb,trail.nb=trail.nb)
detCovs.b <- list(lu.b=lu.b,v.b=v.b,ctr.b=ctr.b,lr.b=lr.b,hco.b=hco.b,trail.b=trail.b)

rm(list=c('lu.nb','v.nb','ctr.nb','lr.nb','hco.nb','trail.nb','lu.b','v.b','ctr.b','lr.b','hco.b','trail.b'))
rm(list = c('ct.nb','ct.b','ynb','yb','r'))
rm('rep.col')