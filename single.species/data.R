
########################
# SET OBSERVATION DATA #
########################

load( here::here("occ.data","occ.RData") )
cams <- y$ylp[,c(1:4)]

###########################
# PREPARE SITE COVARIATES #
###########################
sCovs <- read.csv( here::here("occ.data","siteCovs.csv"),sep=",",header = T) 
sCovs <- sCovs[,-1]

siteCovs<-NULL
siteCovs$station.season <- sCovs$station.season
siteCovs<-data.frame(siteCovs)

## TRANSFORM STANDARDIZE RAW COVARIATES (CONTINUOUS TO Z-SCORES; PROPORTIONAL TO ARCSINE)  ##
# LATITUDE
attach(sCovs)
#latitude
siteCovs$lat <- (utm_y-mean(utm_y))/sd(utm_y)
#water availability
siteCovs$d.wb <- (dst.waterbodies-mean(dst.waterbodies))/sd(dst.waterbodies)
siteCovs$d.wtr <- (dst.water-mean(dst.water))/sd(dst.water)
#human disturbance
siteCovs$d.sett <- (dst_settlement-mean(dst_settlement))/sd(dst_settlement)
siteCovs$d.hmn <- (dst.human-mean(dst.human))/sd(dst.human)
siteCovs$rds500 <- (roads500m-mean(roads500m))/sd(roads500m)
#prey
rbbtts<-NULL
for (i in 1:nrow(sCovs)){
  rbbtts[i] <- ifelse(rabbit.ts[i]==0,0.001,rabbit.ts[i])}
l.rbbtts<-log(rbbtts)

rdtts<-NULL
for (i in 1:nrow(sCovs)) {
  rdtts[i] <- ifelse(rodent.ts[i]==0,0.001,rodent.ts[i])}
l.rdtts<-log(rdtts)

lynxts<-NULL
for (i in 1:nrow(sCovs)) {
  lynxts[i] <- ifelse(lynx.ts[i]==0,0.001,lynx.ts[i])}
l.lynxts<-log(lynxts)

siteCovs$rbbt.ts <- (rabbit.ts-mean(rabbit.ts))/sd(rabbit.ts)
siteCovs$log.rbbt.ts <- (l.rbbtts-mean(l.rbbtts))/sd(l.rbbtts)
siteCovs$rod.ts <- (rodent.ts-mean(rodent.ts))/sd(rodent.ts)
siteCovs$log.rod.ts <- (l.rdtts-mean(l.rdtts))/sd(l.rdtts)
siteCovs$lynx.ts <- (lynx.ts-mean(lynx.ts))/sd(lynx.ts)
siteCovs$log.lynx.ts <- (l.lynxts-mean(l.lynxts))/sd(l.lynxts)

for(i in 1:nrow(sCovs)){
  if(sCovs$season[i]=='breeding'){
    siteCovs$rbbt.ha[i] <- (rbbt.ha.b[i]-mean(rbbt.ha.b))/sd(rbbt.ha.b)
  }else{
    siteCovs$rbbt.ha[i] <- (rbbt.ha.nb[i]-mean(rbbt.ha.nb))/sd(rbbt.ha.nb)
  }
}


for(i in 1:nrow(sCovs)){
  if(sCovs$season[i]=='breeding'){
    siteCovs$rod.ha[i] <- (rdt.ha.b[i]-mean(rdt.ha.b))/sd(rdt.ha.b)
  }else{
    siteCovs$rod.ha[i] <- (rdt.ha.nb[i]-mean(rdt.ha.nb))/sd(rdt.ha.nb)
  }
}

#tree cover
siteCovs$TCnb <- (TC500nb-mean(TC500nb))/sd(TC500nb)
siteCovs$log.TCnb <- (log(TC500nb)-mean(log(TC500nb)))/sd(log(TC500nb))
siteCovs$TCb <- (TC500b-mean(TC500b))/sd(TC500b)
siteCovs$log.TCb <- (log(TC500b)-mean(log(TC500b)))/sd(log(TC500b))
#non-tree cover
siteCovs$nTCnb <- (nTC500nb-mean(nTC500nb))/sd(nTC500nb)
siteCovs$log.nTCnb <- (log(nTC500nb)-mean(log(nTC500nb)))/sd(log(nTC500nb))
siteCovs$nTCb <- (nTC500b-mean(nTC500b))/sd(nTC500b)
siteCovs$log.nTCb <- (log(nTC500b)-mean(log(nTC500b)))/sd(log(nTC500b))
#open areas
siteCovs$nVnb <- (nV500nb-mean(nV500nb))/sd(nV500nb)
siteCovs$log.nVnb <- (log(nV500nb)-mean(log(nV500nb)))/sd(log(nV500nb))
siteCovs$nVb <- (nV500b-mean(nV500b))/sd(nV500b)
siteCovs$log.nVb <- (log(nV500b)-mean(log(nV500b)))/sd(log(nV500b))

#EVI
siteCovs$EVInb <- (EVI500nb-mean(EVI500nb))/sd(EVI500nb)
siteCovs$log.EVInb <- (log(EVI500nb)-mean(log(EVI500nb)))/sd(log(EVI500nb))
siteCovs$EVIb <- (EVI500b-mean(EVI500b))/sd(EVI500b)
siteCovs$log.EVIb <- (log(EVI500b)-mean(log(EVI500b)))/sd(log(EVI500b))

rm(sCovs)

siteCovs <- merge(cams, siteCovs, by='station.season',all.x = T, all.y = F, incomparables = NA)

### ADD LANDCOVER MODIS_IGBP DATASET ###

ct.north <- read.csv( here::here("occ.data", "ctNorth_LCtype1%.csv"), header=T, sep=',') 

ct.south <- read.csv( here::here("occ.data","ctSouth_LCtype1%.csv"), header=T, sep=',') 

ct.north <- ct.north[,c(1:4,7,8,11,12,9,15,10,14,6,13)]
ct.south <- ct.south[,c(1:4,9,13,10,15,7,12,6,11,8,14)]
lc <- rbind(ct.north,ct.south)

lc <- lc[which(lc$study.area=="CNP"|lc$study.area=="GVNP"|lc$study.area=="PGNP"|lc$study.area=="MNR"|lc$study.area=="SANP"|lc$study.area=="MNP"),]

for(i in 1:nrow(lc)){
  lc$lc[i] <- if(lc$study.area[i]=="CNP"|lc$study.area[i]=="GVNP"){lc$LC2009[i]}else{
    if(lc$study.area[i]=="PGNP"|lc$study.area[i]=="MNR"){lc$LC2010[i]}else{
      if(lc$study.area[i]=="SANP"|lc$study.area[i]=="MNP"){lc$LC2012[i]}
    }} 
  lc$Xlc[i] <- if(lc$study.area[i]=="CNP"|lc$study.area[i]=="GVNP"){lc$X.LC2009[i]}else{
    if(lc$study.area[i]=="PGNP"|lc$study.area[i]=="MNR"){lc$X.LC2010[i]}else{
      if(lc$study.area[i]=="SANP"|lc$study.area[i]=="MNP"){lc$X.LC2012[i]}
    }} 
}


## ADD IDENTIFIER FOR BIOGEOGRAPHIC REGION: 0 - MEDITERRANEAN; 1 - TEMPERATE

for(i in 1:nrow(lc)){
  lc$region[i] <- if(lc$study.area[i]=="MNR"|lc$study.area[i]=="PGNP"){1}else{0}
}

for(i in 1:nrow(lc)){
  lc$Frst[i] <- if(lc$lc[i]==1|lc$lc[i]==5){1}else{0}
  lc$Shrub[i] <- if(lc$lc[i]==6|lc$lc[i]==7){1}else{0}
  lc$WoodSavanna[i] <- if(lc$lc[i]==8){1}else{0}
  lc$Savanna[i] <- if(lc$lc[i]==9){1}else{0}
  lc$GrassCrop[i] <- if(lc$lc[i]==10|lc$lc[i]==12|lc$lc[i]==14){1}else{0}
  lc$Open[i] <- if(lc$lc[i]==9|lc$lc[i]==10|lc$lc[i]==12|lc$lc[i]==14){1}else{0}
}
lc <- lc[,c("station","region","Frst","Shrub","WoodSavanna","Savanna","GrassCrop","Open")]


siteCovs<-merge(siteCovs, lc, all.x = T, all.y = F, by = 'station', incomparables=NA)

rm(list = c("lc","ct.north","ct.south"))


###########################################
# LOAD AND PREPARE OBSERVATION COVARIATES #
###########################################

load("occ.data/det.Covs.RData")

rm(cams)
rm(list = c('i','l.rbbtts','l.rdtts','l.lynxts','lynxts','rbbtts','rdtts'))
