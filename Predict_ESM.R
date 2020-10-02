setwd("~/OneDrive/SCS_pico_niche")
#Calculate the relative importance of each envr. factor
library(dismo)
library(gbm)
library(foreach)
library(plot3D)
load('np.Rdata')
source('Calc_bath.R')
set.seed(1)

#ESM files
ESMdir <- 'CMIP6_SST_Chl_PAR/'
title  <- 'CESM2'
CESM2  <- paste0(ESMdir,title,'ChlSSTPAR.nc')

Nrep <- 10  #10 randomizations

#Get grid
LONe <- ncread(CESM2, 'LON')
LATe <- ncread(CESM2, 'LAT')

#Create vertical depth
Bom <- -150 #Maximal depth 150 m for prediction
Hz  <- 10   #Vertical interval of each layer
Z_r <- seq(-Hz/2, Bom, by = -Hz)
Z_w <- seq(0,  Bom, by = -Hz)

#Create DOY
Month_len <- c(31, 28, 31, 30, 31, 30, 31, 31, 30,31,30,31)
DOYe      <- numeric(length(Month_len))

acc_days <- 0
for (i in 1:length(Month_len)){
  DOYe[i] <- as.integer(Month_len[i]/2) +  acc_days
  acc_days<- acc_days + Month_len[i]
}

tDOY <- function(x) cos(x/365 * 2*pi) #Transform DOY

#Transform back
btDOY <- function(x)acos(x)*365/2/pi

#Get Year
Yrs <- ncread(CESM2, 'Year')

#Get surface Chl, PAR, and SST
CHLe <- ncread(CESM2, 'Chl')
PARe <- ncread(CESM2, 'PAR')
SSTe <- ncread(CESM2, 'SST')

#Real lon range: 109-121
#Real lat range: 11-23
klon <- which(LONe >= 110 & LONe <= 120)
LONR <- LONe[klon]
klat <- which(LATe >= 16  & LATe <= 23)
LATR <- LATe[klat]

#Extract the real Chl, PAR and SST
CHLe <- CHLe[klat,klon,,]
PARe <- PARe[klat,klon,,]
SSTe <- SSTe[klat,klon,,]

#Construct a dataframe as input
newdat <- expand.grid(lon   = LONR,
                      lat   = LATR,
                      Depth = Z_r)
newdat$Depth <- abs(newdat$Depth)

#Function calculating total abundance in the Whole SCS
#First need to change unit of the biological variable

#2nd, calculate the surface area (m^2) of each grid
SA <- function(latr, dlon=.5, dlat=.5){
  pi   <- 3.1415926535
  lat1 <- latr - dlat/2
  lat2 <- latr + dlat/2
  
  lat1 <- lat1/180 * pi
  lat2 <- lat2/180 * pi
  R    <- 6371e3  #radius of the earth (m)
  return((pi/180)*R^2 * abs(sin(lat1)-sin(lat2)) *abs(dlon))
}

NLON <- length(LONR)
NLAT <- length(LATR)
NZ   <- length(Z_r)
NDOY <- length(DOYe)
NYR  <- length(Yrs)

#Approximately 111 km = 1 degree at equator
#Volume of all grids

#Set up mesh
LATMesh <- matrix(rep(LATR, length(LONR)), nr = length(LONR), byrow=T)

#surface area
sa <- SA(LATMesh)

#Check sa
image2D(sa, LONR, LATR,
        col  = terrain.colors(20),
        cex.lab = 1.4, cex.axis = 1.4, cex.main = 1.4,
        xlab = "Longitude (ºE)",
        ylab = "Latitude (ºN)",
        main = 'Surface area' )

#Caclulate volume
V  <- array(NA, dim = c(NLON, NLAT, NZ))
for (k in 1:length(Z_r)){
    V[,,k] <- Hz * sa 
}

#Need to remove land
Bot_Depth <- get_bath(newdat)
Bot_Depth <- matrix(Bot_Depth,
                    nr = NLON,
                    nc = NLAT,
                    byrow = F)

#Check bathymetry
image2D(Bot_Depth, LONR, LATR,
            zlim = c(-4000,0),
            col  = cm.colors(20),
            cex.lab = 1.4, cex.axis = 1.4, cex.main = 1.4,
            xlab = "Longitude (ºE)",
            ylab = "Latitude (ºN)",
            main = 'Depth' )

for (i in 1:length(LONR)){
  for (j in 1:length(LATR)){
    
    if (Bot_Depth[i,j] > 0) V[i,j,] <- NA
    kz <- which(Z_w < Bot_Depth[i,j] ) #Find the vertical grid indeces below the bottom depth
    
    if (length(kz) >= 1){
      V[i,j,kz-1] <- NA
    }
  }
}

#Volume of coastal areas (bottom depth < 200 m)
Vcoast <- V

#Volume of oceanic areas (bottom depth >= 200 m)
Vocean <- V

for (i in 1:length(LONR)){
  for (j in 1:length(LATR)){
    if( Bot_Depth[i,j] > -200) {
      Vocean[i,j,] <- NA
    }else{
      Vcoast[i,j,] <- NA
    }
  }
}

#Check volume
image2D(V[,,1], LONR, LATR,
        col  = cm.colors(20),
        cex.lab = 1.4, cex.axis = 1.4, cex.main = 1.4,
        xlab = "Longitude (ºE)",
        ylab = "Latitude (ºN)",
        main = 'Volume' )

#Construct the new data for all years
tmp    <- newdat
newdat <- newdat[0,]
for (i in 1:NYR){
  for (j in 1:NDOY){
    tmp$DOY <- tDOY(DOYe[j])
    chl     <- CHLe[,,j,i]
    chl     <- t(chl)
    
    #fill in NaNs
    #Transform chl to dataframe
    fillNaN <- function(chl){
      #Check if there is NA
      chl[Bot_Depth > 0] <- -999999
      
      ff <- matrix(NA, nr = length(LONR)*length(LATR), nc = 3)
      ff <- as.data.frame(ff)
      names(ff) <- c('lon', 'lat', 'chl')
      ff[,c(1,2)] <- expand.grid(LONR, LATR)
      ff[,3]      <- as.vector(chl)
      
      train       <- ff[ff$chl > 0 & !is.na(ff$chl), c('lon', 'lat')]
      trueclass   <- ff[ff$chl > 0 & !is.na(ff$chl), 'chl']
      test        <- ff[is.na(ff$chl), c('lon', 'lat')]
      ff$chl[is.na(ff$chl)] <- as.numeric(as.character(knn(train, 
                                                           test, trueclass)))
      
      #Convert -999999 to NA
      ff$chl[ff$chl < 0] <- NA
      return(ff$chl)
    }
    #convert matrix (chl) to a vector with masks applied and NaNs filled
    chl    <- fillNaN(chl)
    tmp$logChl0 <- log(rep(chl, length(Z_r)))

    sst    <- SSTe[,,j,i]
    sst    <- t(sst)
    
    #remove the NaN in sst
    sst    <- fillNaN(sst)
    tmp$T0 <- rep(sst, length(Z_r))
    
    par    <- PARe[,,j,i]
    par    <- t(par)
    #remove the NaN of par
    par    <- fillNaN(par)
    tmp$I0 <- rep(par, length(Z_r))
    
    newdat <- rbind(newdat, tmp)
  }
}

Nrep    = 10
predChl = array(NA, dim=c(NLON, NLAT, NZ, NDOY,NYR, Nrep))
predPro = predChl
predSyn = predChl
predPeuk= predChl

#Run BRT model
for (p in 1:Nrep){
  x      <- sample(rownames(np), 0.5*nrow(np))
  Train  <- np[x,]   #Data for training
  
  c_brt  <- gbm.step(data  = Train, 
                     gbm.x = 1:7, 
                     gbm.y = 8,  
                     tree.complexity = 15,
                     learning.rate   = .01,
                     family = "gaussian", 
                     silent = T)
  #Predict Chl
  pred.Chl  <- predict.gbm(c_brt, newdat,
                           n.trees=c_brt$gbm.call$best.trees, 
                           type="response")
  
  #Transform the predicted Chl to an array
  pred.Chl  <- array(pred.Chl, dim=c(NLON, NLAT, NZ, NDOY,NYR))
  
  #Transform back to original unit (mg/m3)
  pred.Chl  <- exp(pred.Chl)
  
  p_brt  <- gbm.step(data=Train, gbm.x = 1:7, gbm.y = 9,  
                     tree.complexity = 15,
                     learning.rate   = .01,
                     family = "gaussian", silent = T)
  
  pred.Pro  <- predict.gbm(p_brt, newdat,
                           n.trees=p_brt$gbm.call$best.trees, 
                           type="response")
  
  #Transform the predicted Pro to an array
  pred.Pro  <- array(pred.Pro, dim=c(NLON, NLAT, NZ, NDOY,NYR))
            
  #Transform back to original unit (cells/mL)
  pred.Pro  <- exp(pred.Pro)
  
  s_brt  <- gbm.step(data=Train, gbm.x = 1:7, gbm.y = 10,  
                     tree.complexity = 15,
                     learning.rate   = .01,
                     family = "gaussian", silent = T)
  
  pred.Syn  <- predict.gbm(s_brt, newdat,
                           n.trees=s_brt$gbm.call$best.trees, 
                           type="response")
  
  #Transform the predicted Syn to an array
  pred.Syn  <- array(pred.Syn, dim=c(NLON, NLAT, NZ, NDOY,NYR))
  
  #Transform back to original unit (cells/mL)
  pred.Syn  <- exp(pred.Syn)
 
  e_brt  <- gbm.step(data=Train, gbm.x = 1:7, gbm.y = 11,  
                     tree.complexity = 15,
                     learning.rate   = .01,
                     family = "gaussian", silent = T)
  
  pred.Peuk  <- predict.gbm(e_brt, newdat,
                            n.trees=e_brt$gbm.call$best.trees, 
                            type="response")
  #Transform the predicted Peuk to an array
  pred.Peuk  <- array(pred.Peuk, 
                      dim=c(NLON, NLAT, NZ, NDOY,NYR))
  
  #Transform back to original unit (cells/mL)
  pred.Peuk  <- exp(pred.Peuk)
  
  predChl[,,,,,p] <- pred.Chl
  predPro[,,,,,p] <- pred.Pro
  predSyn[,,,,,p] <- pred.Syn
  predPeuk[,,,,,p]<- pred.Peuk
  
}

fname <- paste0(title,'Pred.Rdata')
save(predChl, predPro, predSyn, predPeuk, file = fname)

#Calculate total Chl, Pro, Syn and Peuk
Calc_AnnT <- function(title){
  fname <- paste0(title,'Pred.Rdata')
  load(fname)
  TChl    <- array(NA, dim=c(NDOY,NYR, Nrep))
  TChlc   <- TChl
  TChlo   <- TChl
  TPro    <- TChl
  TProc   <- TChl
  TProo   <- TChl
  TSyn    <- TChl
  TPeuk   <- TChl
  TPeukc  <- TChl
  TPeuko  <- TChl
  TSync   <- TChl  #Total coastal Syn
  TSyno   <- TChl  #Total oceanic Syn
  
  for (i in 1:NDOY){
    for (j in 1:NYR){
      for (k in 1:Nrep){
        TChl[i,j,k] <- sum(predChl[,,,i,j,k] * V, na.rm=T)/1e12 #unit: Gg (10^9 g)
        TChlc[i,j,k]<- sum(predChl[,,,i,j,k] * Vcoast, na.rm=T)/1e12
        TChlo[i,j,k]<- sum(predChl[,,,i,j,k] * Vocean, na.rm=T)/1e12
        
        TPro[i,j,k] <- sum(predPro[,,,i,j,k] * V, na.rm=T)/1e14 #unit: 10^20 cells
        TProc[i,j,k]<- sum(predPro[,,,i,j,k] * Vcoast, na.rm=T)/1e14
        TProo[i,j,k]<- sum(predPro[,,,i,j,k] * Vocean, na.rm=T)/1e14
        
        TSyn[i,j,k] <- sum(predSyn[,,,i,j,k] * V, na.rm=T)/1e14
        TSync[i,j,k]<- sum(predSyn[,,,i,j,k] * Vcoast, na.rm=T)/1e14
        TSyno[i,j,k]<- sum(predSyn[,,,i,j,k] * Vocean, na.rm=T)/1e14
        
        TPeuk[i,j,k] <- sum(predPeuk[,,,i,j,k]* V, na.rm=T)/1e14 #unit: 10^20 cells
        TPeukc[i,j,k]<- sum(predPeuk[,,,i,j,k]* Vcoast, na.rm=T)/1e14
        TPeuko[i,j,k]<- sum(predPeuk[,,,i,j,k]* Vocean, na.rm=T)/1e14
      }
    }
  }
  TChlavg <- apply(TChl, c(1,2), mean) #unit: Gg 
  TChlSD  <- apply(TChl, c(1,2), sd)
  TChlcavg <- apply(TChlc, c(1,2), mean)
  TChloavg <- apply(TChlo, c(1,2), mean)
  
  TProavg <- apply(TPro, c(1,2), mean)/1e4 #unit: 10^24 (Septillion) cells, 
  #total length 1e24 * 0.6*1e-6 = 60 light years
  TProcavg <- apply(TProc, c(1,2), mean)/1e4
  TProoavg <- apply(TProo, c(1,2), mean)/1e4
  
  TProSD  <- apply(TPro, c(1,2), sd)/1e4
  
  TSynavg <- apply(TSyn, c(1,2), mean)/1e4 
  TSynSD  <- apply(TSyn, c(1,2), sd)/1e4
  
  TSyncavg <- apply(TSync, c(1,2), mean)/1e4 
  TSyncSD  <- apply(TSync, c(1,2), sd)/1e4
  
  TSynoavg <- apply(TSyno, c(1,2), mean)/1e4 
  TSynoSD  <- apply(TSyno, c(1,2), sd)/1e4
  
  TPeukavg <- apply(TPeuk, c(1,2), mean)/1e4 
  TPeukSD  <- apply(TPeuk, c(1,2), sd)/1e4
  TPeukcavg <- apply(TPeukc, c(1,2), mean)/1e4
  TPeukoavg <- apply(TPeuko, c(1,2), mean)/1e4
  
  AnnTChl  <- apply(TChlavg, 2, mean)
  AnnTChlc  <- apply(TChlcavg, 2, mean)
  AnnTChlo  <- apply(TChloavg, 2, mean)
  
  AnnTPro  <- apply(TProavg, 2, mean)
  AnnTProc  <- apply(TProcavg, 2, mean)
  AnnTProo  <- apply(TProoavg, 2, mean)
  
  AnnTSyn  <- apply(TSynavg, 2, mean)
  AnnTSync  <- apply(TSyncavg, 2, mean)
  AnnTSyno  <- apply(TSynoavg, 2, mean)
  
  AnnTPeuk <- apply(TPeukavg, 2, mean)
  AnnTPeukc <- apply(TPeukcavg, 2, mean)
  AnnTPeuko <- apply(TPeukoavg, 2, mean)
  
  return(list(Chl = AnnTChl, Chlc = AnnTChlc, Chlo = AnnTChlo, 
              Pro = AnnTPro, Proc = AnnTProc, Proo = AnnTProo,
              Syn = AnnTSyn, Sync = AnnTSync, Syno = AnnTSyno, 
              Peuk = AnnTPeuk, Peukc = AnnTPeukc, Peuko = AnnTPeuko))
}

CESM2Pred   <- Calc_AnnT('CESM2')

pdf("ESM_pred1Oct.pdf", width=7, height=7)

op   <- par(font.lab=1,
            family ="serif",
            mar    = c(2,4,2,2),
            mgp    = c(2,1,0),
            mfrow  = c(2,2),
            oma    = c(3,3,0,0),
            cex.axis=1.2, 
            cex.lab =1.2)

Ymin <- min(CESM2Pred$Chl, CESM2Pred$Chlc, CESM2Pred$Chlo, na.rm=T)
Ymax <- max(CESM2Pred$Chl, CESM2Pred$Chlc, CESM2Pred$Chlo, na.rm=T)

plot(Yrs, CESM2Pred$Chl, 
     ylim = c(Ymin, Ymax+2),
     type = 'b',
     xlab = '',
     ylab = bquote('Annual Mean Chl '*italic(a)*' (Gg)'))
abline(lm(CESM2Pred$Chl ~ Yrs), lwd = 1.5) #p < 0.05, significantly decreasing
points(Yrs, CESM2Pred$Chlc, type = 'b', col = 2)
points(Yrs, CESM2Pred$Chlo, type = 'b', col = 3)
abline(lm(CESM2Pred$Chlc ~ Yrs), lwd = 1.5, col = 2)  #Significantly decreasing
abline(lm(CESM2Pred$Chlo ~ Yrs), lwd = 1.5, col = 3)  #Significantly decreasing
text(2060, 8, bquote(italic(p)*' < 0.001'), pos = 4)
mtext('A) Chl', adj = -0.1, cex=1.4)
legend('topright', c('Total', 'Coastal', 'Oceanic'), lty = 1, lwd = 1.5, col = 1:3)

Ymin <- min(CESM2Pred$Pro, CESM2Pred$Proc, CESM2Pred$Proo, na.rm=T)
Ymax <- max(CESM2Pred$Pro, CESM2Pred$Proc, CESM2Pred$Proo, na.rm=T)
plot(Yrs, CESM2Pred$Pro, 
     ylim = c(Ymin, Ymax),
     type = 'b',
     xlab = '',
     ylab = bquote('Annual Mean abundance ('*10^24* ' cells)'))
points(Yrs, CESM2Pred$Proc, type = 'b', col = 2)
points(Yrs, CESM2Pred$Proo, type = 'b', col = 3)
abline(lm(CESM2Pred$Proc ~ Yrs), lwd = 1.5, col = 2) #Significantly increasing
abline(lm(CESM2Pred$Proo ~ Yrs), lwd = 1.5, col = 3) 

abline(lm(CESM2Pred$Pro ~ Yrs), lwd = 1.5) #p < 0.05, significantly increasing
text(2060, 3, bquote(italic(p)*' < 0.001'), pos = 4)
mtext('B) Pro', adj = -0.1, cex=1.4)

#Calculate the limits of yaxis of Syn
Ymin <- min(CESM2Pred$Syn, CESM2Pred$Sync, CESM2Pred$Syno, na.rm=T)
Ymax <- max(CESM2Pred$Syn, CESM2Pred$Sync, CESM2Pred$Syno, na.rm=T)
plot(Yrs, CESM2Pred$Syn, 
     ylim = c(Ymin, Ymax),
     type = 'b',
     xlab = '',
     ylab = bquote('Annual Mean abundance ('*10^24* ' cells)'))
points(Yrs, CESM2Pred$Sync, type = 'b', col = 2)
points(Yrs, CESM2Pred$Syno, type = 'b', col = 3)
abline(lm(CESM2Pred$Syn  ~ Yrs), lwd = 1.5)  #p > 0.05
abline(lm(CESM2Pred$Sync ~ Yrs), lwd = 1.5, col = 2)  #p < 0.05, significantly increasing
abline(lm(CESM2Pred$Syno ~ Yrs), lwd = 1.5, col = 3)  #p > 0.05
text(2020, .198, bquote(italic(p)*' < 0.05'), pos = 4, col = 2)
text(2020, .475,   bquote(italic(p)*' > 0.05'), pos = 4, col = 1)
mtext('C) Syn', adj = -0.1, cex=1.4)

Ymin <- min(CESM2Pred$Peuk, CESM2Pred$Peukc, CESM2Pred$Peuko, na.rm=T)
Ymax <- max(CESM2Pred$Peuk, CESM2Pred$Peukc, CESM2Pred$Peuko, na.rm=T)
plot(Yrs, CESM2Pred$Peuk, 
     ylim = c(Ymin, Ymax),
     type = 'b',
     xlab = '',
     ylab = bquote('Annual Mean abundance ('*10^24* ' cells)'))

points(Yrs, CESM2Pred$Peukc, type = 'b', col = 2)
points(Yrs, CESM2Pred$Peuko, type = 'b', col = 3)

abline(lm(CESM2Pred$Peuk ~ Yrs), lwd = 1.5)   #p > 0.05
abline(lm(CESM2Pred$Peukc ~ Yrs), lwd = 1.5, col = 2)  #p > 0.05
abline(lm(CESM2Pred$Peuko ~ Yrs), lwd = 1.5, col = 3)  
text(2060, .1, bquote(italic(p)*' > 0.05'), pos = 4)
mtext('D) Peuk', adj = -0.1, cex=1.4)

mtext('Year', side = 1, outer = T, line = .5, cex = 1.2)
dev.off()

#Make a video
# today  = Sys.Date()
# folder = paste0('Output',today)
# unlink(folder, recursive = T)
# dir.create(folder)
# 
# Months <- c('Jan','Feb', 'Mar','Apr', 'May','Jun', 'Jul','Aug',
#             'Sep','Oct', 'Nov', 'Dec')
# #Calculate avg
# CHLA <- apply(predChl, c(1:5), mean)
# PROA <- apply(predPro, c(1:5), mean)
# SYNA <- apply(predSyn, c(1:5), mean)
# PEUA <- apply(predPeuk,c(1:5), mean)
# 
# Chlrange <- quantile(CHLA, probs=c(0.025,0.975))
# PROrange <- quantile(PROA, probs=c(0.025,0.975))
# SYNrange <- quantile(SYNA, probs=c(0.025,0.975))
# PEUrange <- quantile(PEUA, probs=c(0.025,0.975))
# #plot graph using image2D
# jet.colors   <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
#                                    "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")) 
# 
# LL <- 0
# for (i in 1:NYR){
#   for (j in 1:NDOY){
#     Year  <- sprintf("%02d",  i)
#     Mon   <- sprintf("%02d",  j)
#     LL    <- LL + 1
#     LL    <- sprintf("%04d",  j)
#     Chlname <- paste0(folder,'/Chl',LL,'.png')
#     png(Chlname, width=1800, height=1800, res=300)
#     op <- par(font.lab = 1, 
#               cex    = 1.4,
#               family ="serif",
#               mar    = c(4,4,4,2),
#               mgp    = c(2.3,1,0),
#               mfrow  = c(1,1)   ) 
#     
#     COLS <- jet.colors(24)
#     COLKEY <- T
#     chla   <- CHLA[,,1,j,i]
#     chla[Bot_Depth >= 0]  <- NA  #Apply mask
#     chla[chla > Chlrange[2]] <- Chlrange[2]
#     chla[chla < Chlrange[1]] <- Chlrange[1]
#     image2D(chla, LONR, LATR,
#             zlim = Chlrange,
#             col  = COLS,
#             colkey=COLKEY,
#             cex.lab = 1.4, cex.axis = 1.4, cex.main = 1.4,
#             clab = bquote("Chl (mg "*m^-3*")"),
#             xlab = "Longitude (ºE)", 
#             ylab = "Latitude (ºN)",
#             main = paste0(Yrs[i], ' ', Months[j]) )
#     
#     dev.off()
#     
#     Proname <- paste0(folder,'/Pro',LL,'.png')
#     png(Proname, width=1800, height=1800, res=300)
#     op <- par(font.lab = 1, 
#               cex    = 1.4,
#               family ="serif",
#               mar    = c(4,4,4,2),
#               mgp    = c(2.3,1,0),
#               mfrow  = c(1,1)   ) 
#     
#     pro   <- PROA[,,1,j,i]
#     pro[Bot_Depth >= 0]  <- NA  #Apply mask
#     pro[pro > PROrange[2]] <- PROrange[2]
#     pro[pro < PROrange[1]] <- PROrange[1]
#     image2D(pro, LONR, LATR,
#             zlim = PROrange,
#             col  = COLS,
#             colkey=COLKEY,
#             cex.lab = 1.4, cex.axis = 1.4, cex.main = 1.4,
#             clab = bquote("Pro (cells "*mL^-1*")"),
#             xlab = "Longitude (ºE)", 
#             ylab = "Latitude (ºN)",
#             main = paste0(Yrs[i], ' ', Months[j]) )
#     
#     dev.off()
#   }
# }
#setwd(folder)
#system('ffmpeg -framerate 4 -i Chl%4d.png ChlESM.avi')
#setwd('../')
