
########################
# SET OBSERVATION DATA #
########################

#Red fox
yvv.nb <- read.csv("/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/1.data/occ.data/yvv.nb.csv",sep=",",header = T)
yvv.nb <- yvv.nb[,-1]

#Stone marten
ymf.nb <- read.csv("/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/1.data/occ.data/ymf.nb.csv",sep=",",header = T)
ymf.nb <- ymf.nb[,-1]

#European wildcat
yfs.nb <- read.csv("/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/1.data/occ.data/yfs.nb.csv",sep=",",header = T)
yfs.nb <- yfs.nb[,-1]

#Common genet
ygg.nb <- read.csv("/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/1.data/occ.data/ygg.nb.csv",sep=",",header = T)
ygg.nb <- ygg.nb[,-1]

#Pine marten
ymmt.nb <- read.csv("/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/1.data/occ.data/ymmt.nb.csv",sep=",",header = T)
ymmt.nb <- ymmt.nb[,-1]

#Egyptian mongoose
yhi.nb <- read.csv("/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/1.data/occ.data/yi.nb.csv",sep=",",header = T)
yhi.nb <- yhi.nb[,-1]

#Eurasian badger
ymm.nb <- read.csv("/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/1.data/occ.data/ymm.nb.csv",sep=",",header = T)
ymm.nb <- ymm.nb[,-1]

#Iberian lynx
ylp.nb <- read.csv("/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/1.data/occ.data/ylp.nb.csv",sep=",",header = T)
ylp.nb <- ylp.nb[,-1]

#study areas and staton codes
y.nb <- yvv.nb[,c(1:3)]

###########################
# PREPARE SITE COVARIATES #
###########################
sCovs.nb <- read.csv("/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/1.data/occ.data/temp/siteCovs.nb.csv",sep=",",header = T)
sCovs.nb <- sCovs.nb[,-1]

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
#siteCovs.nb$d.rivers <- (dst_rivers-mean(dst_rivers))/sd(dst_rivers)
#siteCovs.nb$d.streams <- (dst_streams-mean(dst_streams))/sd(dst_streams)
siteCovs.nb$d.wtr <- (dst.water-mean(dst.water))/sd(dst.water)
#human disturbance
siteCovs.nb$d.sett <- (dst_settlement-mean(dst_settlement))/sd(dst_settlement)
#siteCovs.nb$d.road <- (dst_roads-mean(dst_roads))/sd(dst_roads)
siteCovs.nb$d.hmn <- (dst.human-mean(dst.human))/sd(dst.human)
siteCovs.nb$rds1000 <- (roads1000m-mean(roads1000m))/sd(roads1000m)
siteCovs.nb$rds500 <- (roads500m-mean(roads500m))/sd(roads500m)
siteCovs.nb$rds200 <- (roads200m-mean(roads200m))/sd(roads200m)
#prey
rbbtts<-NULL
for (i in 1:nrow(sCovs.nb)) {
  rbbtts[i] <- ifelse(rabbit.ts[i]==0,0.001,rabbit.ts[i])
}
l.rbbtts<-log(rbbtts)

rdtts<-NULL
for (i in 1:nrow(sCovs.nb)) {
  rdtts[i] <- ifelse(rodent.ts[i]==0,0.001,rodent.ts[i])
}
l.rdtts<-log(rdtts)

siteCovs.nb$rbbtts <- (rabbit.ts-mean(rabbit.ts, na.rm = T))/sd(rabbit.ts, na.rm = T)
siteCovs.nb$log.rbbtts <- (l.rbbtts-mean(l.rbbtts, na.rm = T))/sd(l.rbbtts, na.rm = T)

siteCovs.nb$rod.ts <- (rodent.ts-mean(rodent.ts, na.rm = T))/sd(rodent.ts, na.rm = T)
siteCovs.nb$log.rods <- (l.rdtts-mean(l.rdtts, na.rm = T))/sd(l.rdtts, na.rm = T)

siteCovs.nb$rbbt.ha <- (rbbt.ha-mean(rbbt.ha, na.rm = T))/sd(rbbt.ha, na.rm = T)
siteCovs.nb$rdt.ha <- (rdt.ha-mean(rdt.ha, na.rm = T))/sd(rdt.ha, na.rm = T)

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

### ADD LANDCOVER MODIS_IGBP DATASET ###

ct.north <- read.csv("/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/1.data/ctNorth_LCtype1%.csv", header=T, sep=',')

ct.south <- read.csv("/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/1.data/ctSouth_LCtype1%.csv", header=T, sep=',')

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


siteCovs.nb<-merge(siteCovs.nb, lc, all.x = T, all.y = F, by = 'station', incomparables=NA)

rm(list = c("lc","ct.north","ct.south"))



###########################################
# LOAD AND PREPARE OBSERVATION COVARIATES #
###########################################

ynb <- read.csv('/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/1.data/occ.data/temp/siteCovs.nb.csv', sep = ',', header = T)
ynb <- ynb[,-1]

ct.nb <- data.frame(ynb[,c(2,1)])
colnames(ct.nb) <- c("station","study.area")

require(gdata)

# upload detection covariates
### lynx urine
lu.nb<-read.xls(xls = '/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/5.detection covariates/detCovs.nb.xls',sheet = 'lynx.urine')

lu.nb <- merge(ct.nb, lu.nb, all.x = T, by = 'station', all.y = F, incomparables = NA)
lu.nb <- lu.nb[,c(1:2,4:34)]
colnames(lu.nb)[2] <- 'study.area'

### valerian
v.nb<-read.xls(xls = '/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/5.detection covariates/detCovs.nb.xls',sheet = 'valerian')

v.nb <- merge(ct.nb, v.nb, all.x = T, by = 'station', all.y = F, incomparables = NA)
v.nb <- v.nb[,c(1:2,4:33)]
colnames(v.nb)[2] <- 'study.area'

### no.attractant
ctr.nb<-read.xls(xls = '/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/5.detection covariates/detCovs.nb.xls',sheet = 'no.attract')

ctr.nb <- merge(ct.nb, ctr.nb, all.x = T, by = 'station', all.y = F, incomparables = NA)
ctr.nb <- ctr.nb[,c(1:2,4:33)]
colnames(ctr.nb)[2] <- 'study.area'

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

### camera model: leafriver
lr.nb<-read.xls(xls = '/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/5.detection covariates/detCovs.nb.xls',sheet = 'leafriver')

lr.nb <- merge(ct.nb, lr.nb, all.x = T, by = 'station', all.y = F, incomparables = NA)
lr.nb <- lr.nb[,c(1:2,4:33)]
colnames(lr.nb)[2] <- 'study.area'

### camera model: HCO
hco.nb<-read.xls(xls = '/Users/pedromonterroso/Documents/Trabalho/Publications/Papers/In prep/1.spatiotemporal interactions among Iberian carnivores/3.analyses/5.detection covariates/detCovs.nb.xls',sheet = 'HCO')

hco.nb <- merge(ct.nb, hco.nb, all.x = T, by = 'station', all.y = F, incomparables = NA)
hco.nb <- hco.nb[,c(1:2,4:33)]
colnames(hco.nb)[2] <- 'study.area'

### Join all Detection covariates in a single list object
detCovs.nb <- list(lu.nb=lu.nb,v.nb=v.nb,ctr.nb=ctr.nb,lr.nb=lr.nb,hco.nb=hco.nb,trail.nb=trail.nb)


rm(list=c('lu.nb','v.nb','ctr.nb','lr.nb','hco.nb','trail.nb'))
rm(list = c('ct.nb','ynb','r'))
rm('rep.col')

#SPECIES ENCOUNTER HISTORIES
EHvv.nb <- yvv.nb[,4:33]
EHmf.nb <- ymf.nb[,4:33]
EHfs.nb <- yfs.nb[,4:33]
EHgg.nb <- ygg.nb[,4:33]
EHmmt.nb <- ymmt.nb[,4:33]
EHhi.nb <- yhi.nb[,4:33]
EHmm.nb <- ymm.nb[,4:33]
EHlp.nb <- ylp.nb[,4:33]

eff.nb <- yvv.nb[,4:33]
#CONVERT DETECTION HISTORIES FOR 5-DAY SAMPLING OCCASIONS
#Red fox
EHvv.nb[31:36] <- NA
EHvv.nb <- sapply( EHvv.nb, as.numeric )

for(j in 31:36){
  for(i in 1:nrow(EHvv.nb)){
    EHvv.nb[i, j] <- if(sum(is.na(EHvv.nb[i, c((j/j+5*(j-31)):((j/j+5*(j-31))+4))]))==5){NA} else{
      if(sum(EHvv.nb[i, (j/j+5*(j-31))],
             EHvv.nb[i, (j/j+5*(j-31))+1],
             EHvv.nb[i, (j/j+5*(j-31))+2],
             EHvv.nb[i, (j/j+5*(j-31))+3],
             EHvv.nb[i, (j/j+5*(j-31))+4], na.rm=T)>0){1}
      else{0}
    }
  }
}

EHvv.nb <- data.frame(EHvv.nb[,c(31:36)])

#Stone marten
EHmf.nb[31:36] <- NA
EHmf.nb <- sapply( EHmf.nb, as.numeric )

for(j in 31:36){
  for(i in 1:nrow(EHmf.nb)){
    EHmf.nb[i, j] <- if(sum(is.na(EHmf.nb[i, c((j/j+5*(j-31)):((j/j+5*(j-31))+4))]))==5){NA} else{
      if(sum(EHmf.nb[i, (j/j+5*(j-31))],
             EHmf.nb[i, (j/j+5*(j-31))+1],
             EHmf.nb[i, (j/j+5*(j-31))+2],
             EHmf.nb[i, (j/j+5*(j-31))+3],
             EHmf.nb[i, (j/j+5*(j-31))+4], na.rm=T)>0){1}
      else{0}
    }
  }
}

EHmf.nb <- data.frame(EHmf.nb[,c(31:36)])

#European wildcat
EHfs.nb[31:36] <- NA
EHfs.nb <- sapply( EHfs.nb, as.numeric )

for(j in 31:36){
  for(i in 1:nrow(EHfs.nb)){
    EHfs.nb[i, j] <- if(sum(is.na(EHfs.nb[i, c((j/j+5*(j-31)):((j/j+5*(j-31))+4))]))==5){NA} else{
      if(sum(EHfs.nb[i, (j/j+5*(j-31))],
             EHfs.nb[i, (j/j+5*(j-31))+1],
             EHfs.nb[i, (j/j+5*(j-31))+2],
             EHfs.nb[i, (j/j+5*(j-31))+3],
             EHfs.nb[i, (j/j+5*(j-31))+4], na.rm=T)>0){1}
      else{0}
    }
  }
}

EHfs.nb <- data.frame(EHfs.nb[,c(31:36)])

#Common genet
EHgg.nb[31:36] <- NA
EHgg.nb <- sapply( EHgg.nb, as.numeric )

for(j in 31:36){
  for(i in 1:nrow(EHgg.nb)){
    EHgg.nb[i, j] <- if(sum(is.na(EHgg.nb[i, c((j/j+5*(j-31)):((j/j+5*(j-31))+4))]))==5){NA} else{
      if(sum(EHgg.nb[i, (j/j+5*(j-31))],
             EHgg.nb[i, (j/j+5*(j-31))+1],
             EHgg.nb[i, (j/j+5*(j-31))+2],
             EHgg.nb[i, (j/j+5*(j-31))+3],
             EHgg.nb[i, (j/j+5*(j-31))+4], na.rm=T)>0){1}
      else{0}
    }
  }
}

EHgg.nb <- data.frame(EHgg.nb[,c(31:36)])

#Pine marten
EHmmt.nb[31:36] <- NA
EHmmt.nb <- sapply( EHmmt.nb, as.numeric )

for(j in 31:36){
  for(i in 1:nrow(EHmmt.nb)){
    EHmmt.nb[i, j] <- if(sum(is.na(EHmmt.nb[i, c((j/j+5*(j-31)):((j/j+5*(j-31))+4))]))==5){NA} else{
      if(sum(EHmmt.nb[i, (j/j+5*(j-31))],
             EHmmt.nb[i, (j/j+5*(j-31))+1],
             EHmmt.nb[i, (j/j+5*(j-31))+2],
             EHmmt.nb[i, (j/j+5*(j-31))+3],
             EHmmt.nb[i, (j/j+5*(j-31))+4], na.rm=T)>0){1}
      else{0}
    }
  }
}

EHmmt.nb <- data.frame(EHmmt.nb[,c(31:36)])

#Egyptian mongoose
EHhi.nb[31:36] <- NA
EHhi.nb <- sapply( EHhi.nb, as.numeric )

for(j in 31:36){
  for(i in 1:nrow(EHhi.nb)){
    EHhi.nb[i, j] <- if(sum(is.na(EHhi.nb[i, c((j/j+5*(j-31)):((j/j+5*(j-31))+4))]))==5){NA} else{
      if(sum(EHhi.nb[i, (j/j+5*(j-31))],
             EHhi.nb[i, (j/j+5*(j-31))+1],
             EHhi.nb[i, (j/j+5*(j-31))+2],
             EHhi.nb[i, (j/j+5*(j-31))+3],
             EHhi.nb[i, (j/j+5*(j-31))+4], na.rm=T)>0){1}
      else{0}
    }
  }
}

EHhi.nb <- data.frame(EHhi.nb[,c(31:36)])

#Eurasian badger
EHmm.nb[31:36] <- NA
EHmm.nb <- sapply( EHmm.nb, as.numeric )

for(j in 31:36){
  for(i in 1:nrow(EHmm.nb)){
    EHmm.nb[i, j] <- if(sum(is.na(EHmm.nb[i, c((j/j+5*(j-31)):((j/j+5*(j-31))+4))]))==5){NA} else{
      if(sum(EHmm.nb[i, (j/j+5*(j-31))],
             EHmm.nb[i, (j/j+5*(j-31))+1],
             EHmm.nb[i, (j/j+5*(j-31))+2],
             EHmm.nb[i, (j/j+5*(j-31))+3],
             EHmm.nb[i, (j/j+5*(j-31))+4], na.rm=T)>0){1}
      else{0}
    }
  }
}

EHmm.nb <- data.frame(EHmm.nb[,c(31:36)])

#Iberian lynx
EHlp.nb[31:36] <- NA
EHlp.nb <- sapply( EHlp.nb, as.numeric )

for(j in 31:36){
  for(i in 1:nrow(EHlp.nb)){
    EHlp.nb[i, j] <- if(sum(is.na(EHlp.nb[i, c((j/j+5*(j-31)):((j/j+5*(j-31))+4))]))==5){NA} else{
      if(sum(EHlp.nb[i, (j/j+5*(j-31))],
             EHlp.nb[i, (j/j+5*(j-31))+1],
             EHlp.nb[i, (j/j+5*(j-31))+2],
             EHlp.nb[i, (j/j+5*(j-31))+3],
             EHlp.nb[i, (j/j+5*(j-31))+4], na.rm=T)>0){1}
      else{0}
    }
  }
}

EHlp.nb <- data.frame(EHlp.nb[,c(31:36)])

#NAIVE OCCUPANCY ESTIMATES
naive.vv.nb <- sum(apply(EHvv.nb, 1, sum, na.rm=T) > 0) / nrow(EHvv.nb) 
naive.mf.nb <- sum(apply(EHmf.nb, 1, sum, na.rm=T) > 0) / nrow(EHmf.nb) 
naive.fs.nb <- sum(apply(EHfs.nb, 1, sum, na.rm=T) > 0) / nrow(EHfs.nb) 
naive.gg.nb <- sum(apply(EHgg.nb, 1, sum, na.rm=T) > 0) / nrow(EHgg.nb) 
naive.mmt.nb <- sum(apply(EHmmt.nb, 1, sum, na.rm=T) > 0) / nrow(EHmmt.nb) 
naive.hi.nb <- sum(apply(EHhi.nb, 1, sum, na.rm=T) > 0) / nrow(EHhi.nb) 
naive.mm.nb <- sum(apply(EHmm.nb, 1, sum, na.rm=T) > 0) / nrow(EHmm.nb) 
naive.lp.nb <- sum(apply(EHlp.nb, 1, sum, na.rm=T) > 0) / nrow(EHlp.nb) 

#CONVERT DETECTION COVARIATES FOR 5-DAY SAMPLING OCCASIONS
#LEAFRIVER CAMERAS
lr.nb <- detCovs.nb$lr.nb[,c(3:32)]
lr.nb[31:36] <- NA
lr.nb <- sapply( lr.nb, as.numeric )

for(j in 31:36){
  for(i in 1:nrow(lr.nb)){
    lr.nb[i, j] <- sum(lr.nb[i, (j/j+5*(j-31))],
                    lr.nb[i, (j/j+5*(j-31))+1],
                    lr.nb[i, (j/j+5*(j-31))+2],
                    lr.nb[i, (j/j+5*(j-31))+3],
                    lr.nb[i, (j/j+5*(j-31))+4], na.rm=T)/5
  }
}

lr.nb <- data.frame(lr.nb[,c(31:36)])

#HCO CAMERAS
hco.nb <- detCovs.nb$hco.nb[,c(3:32)]
hco.nb[31:36] <- NA
hco.nb <- sapply( hco.nb, as.numeric )

for(j in 31:36){
  for(i in 1:nrow(hco.nb)){
    hco.nb[i, j] <- sum(hco.nb[i, (j/j+5*(j-31))],
                     hco.nb[i, (j/j+5*(j-31))+1],
                     hco.nb[i, (j/j+5*(j-31))+2],
                     hco.nb[i, (j/j+5*(j-31))+3],
                     hco.nb[i, (j/j+5*(j-31))+4], na.rm=T)/5
  }
}

hco.nb <- data.frame(hco.nb[,c(31:36)])


#OBTAIN SAMPLING EFFORT FOR EACH OF THE 5-DAY SAMPLING OCCASIONS
for(j in 31:36){
  for(i in 1:nrow(eff.nb)){
    eff.nb[i, j] <- (5-sum(is.na(eff.nb[i, c((j/j+5*(j-31)):((j/j+5*(j-31))+4))])))/5
  }
}

eff.nb <- data.frame(eff.nb[,c(31:36)])


#NUMBER OF SAMPLING OCCASIONS AND SAMPLING SITES
nocc.nb <- ncol(eff.nb)
nSites.nb <- nrow(eff.nb)

rm(list = c('yvv.nb','ymf.nb','yfs.nb','ygg.nb','ymmt.nb','yhi.nb','ymm.nb','ylp.nb'))
