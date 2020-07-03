setwd("~/OneDrive/SCS_pico_niche")
rm(list=ls())
library(plot3D)
library(dismo)
library(gbm)
load('np.Rdata')
set.seed(1)

tDOY <- function(x) cos(x/365 * 2*pi) #Transform DOY

#Transform back
btDOY <- function(x)acos(x)*365/2/pi

#Use tree complexity = 15 and learning rate = 0.01 to examine individual effects
SEATSLon <- 116
SEATSLat <- 18
Test_DOY <- tDOY(15) #Test winter

#First, test interaction of lon and lat
#set up new dataframe
lonx    <- seq(min(np$lon), max(np$lon), by = .1)
latx    <- seq(min(np$lat), max(np$lat), by = .1)
lonlatx <- expand.grid(lonx, latx)
names(lonlatx) <- c('lon','lat')
newdat_lonlat <- data.frame(lon   = lonlatx[,1],
                            lat   = lonlatx[,2],
                            Depth = 50,
                            DOY   = Test_DOY, 
                            logChl0 = median(np$logChl0),
                            T0    = median(np$T0),
                            I0    = median(np$I0))

Nrep     <- 10
New_Chl  <- matrix(NA, nr = nrow(newdat_lonlat), nc = Nrep)
New_Pro  <- New_Chl
New_Syn  <- New_Chl
New_Peuk <- New_Chl


#2nd new dataframe (Depth)
depx <- seq(0, -150, by = -1)

newdat_Dep <- data.frame(lon   = SEATSLon, #SEATS
                         lat   = SEATSLat,  #SEATS
                         Depth = abs(depx),
                         DOY   = Test_DOY, #median(np$DOY),
                         logChl0 = median(np$logChl0),
                         T0    = median(np$T0),
                         I0    = median(np$I0))
New_Chl1  <- matrix(NA, 
                    nr = nrow(newdat_Dep),
                    nc = Nrep)

New_Pro1  <- New_Chl1
New_Syn1  <- New_Chl1
New_Peuk1 <- New_Chl1

#3rd new dataframe (DOY)
doy <- seq(1, 365, 1)
newdat_DOY <- data.frame(lon   = SEATSLon, #SEATS
                         lat   = SEATSLat,  #SEATS
                         Depth = 50,
                         DOY   = tDOY(doy),
                         logChl0 = median(np$logChl0),
                         T0    = median(np$T0),
                         I0    = median(np$I0))
New_Chl2  <- matrix(NA, 
                    nr = nrow(newdat_DOY),
                    nc = Nrep)

New_Pro2  <- New_Chl2
New_Syn2  <- New_Chl2
New_Peuk2 <- New_Chl2

#4th dataframe (Chl)
newdat_Chl <- data.frame(lon   = SEATSLon, #SEATS
                         lat   = SEATSLat,  #SEATS
                         Depth = 50,
                         DOY   = Test_DOY,
                         logChl0 = seq(log(.01), log(5), length.out=100),
                         T0    = median(np$T0),
                         I0    = median(np$I0))
New_Chl3  <- matrix(NA, 
                    nr = nrow(newdat_Chl),
                    nc = Nrep)

New_Pro3  <- New_Chl3
New_Syn3  <- New_Chl3
New_Peuk3 <- New_Chl3

#5th dataframe (SST)
newdat_SST <- data.frame(lon   = SEATSLon, #SEATS
                         lat   = SEATSLat,  #SEATS
                         Depth = 50,
                         DOY   = Test_DOY,
                         logChl0 = median(np$logChl0),
                         T0    = seq(quantile(np$T0, probs=.01), 
                                     quantile(np$T0, probs=.99), 
                                     length.out = 100),
                         I0    = median(np$I0))
New_Chl4  <- matrix(NA, 
                    nr = nrow(newdat_SST),
                    nc = Nrep)

New_Pro4  <- New_Chl4
New_Syn4  <- New_Chl4
New_Peuk4 <- New_Chl4

#6th dataframe (I0)
newdat_I0 <- data.frame(lon   = SEATSLon, #SEATS
                        lat   = SEATSLat,  #SEATS
                        Depth = 50,
                        DOY   = Test_DOY,
                        logChl0 = median(np$logChl0),
                        T0    = median(np$T0),
                        I0    = seq(quantile(np$I0,probs=.01), 
                                    quantile(np$I0,probs=.99), 
                                    length.out=100))
New_Chl5  <- matrix(NA, 
                    nr = nrow(newdat_I0),
                    nc = Nrep)

New_Pro5  <- New_Chl5
New_Syn5  <- New_Chl5
New_Peuk5 <- New_Chl5

for(i in 1:Nrep) {
  x      <- sample(rownames(np), 0.5*nrow(np))
  Train  <- np[x,]   #Data for training
  
  c_brt  <- gbm.step(data  = Train, 
                     gbm.x = 1:7, 
                     gbm.y = 8,  
                     tree.complexity = 15,
                     learning.rate   = .01,
                     family = "gaussian", 
                     silent = T)
  
  New_Chl[,i]  <- predict.gbm(c_brt, newdat_lonlat,
                              n.trees=c_brt$gbm.call$best.trees, 
                              type="response")
  
  New_Chl1[,i] <- predict.gbm(c_brt, newdat_Dep,
                              n.trees=c_brt$gbm.call$best.trees, 
                              type="response")
  
  New_Chl2[,i] <- predict.gbm(c_brt, newdat_DOY,
                              n.trees=c_brt$gbm.call$best.trees, 
                              type="response")
  
  New_Chl3[,i] <- predict.gbm(c_brt, newdat_Chl,
                              n.trees=c_brt$gbm.call$best.trees, 
                              type="response") 
  
  New_Chl4[,i] <- predict.gbm(c_brt, newdat_SST,
                              n.trees=c_brt$gbm.call$best.trees, 
                              type="response") 
  
  New_Chl5[,i] <- predict.gbm(c_brt, newdat_I0,
                              n.trees=c_brt$gbm.call$best.trees, 
                              type="response") 
  
  p_brt  <- gbm.step(data=Train, gbm.x = 1:7, gbm.y = 9,  
                     tree.complexity = 15,
                     learning.rate   = .01,
                     family = "gaussian", silent = T)
  
  New_Pro[,i]    <- predict.gbm(p_brt, newdat_lonlat,
                                n.trees=p_brt$gbm.call$best.trees, 
                                type="response")
  
  New_Pro1[,i]   <- predict.gbm(p_brt, newdat_Dep,
                                n.trees=p_brt$gbm.call$best.trees, 
                                type="response")
  
  New_Pro2[,i]   <- predict.gbm(p_brt, newdat_DOY,
                                n.trees=p_brt$gbm.call$best.trees, 
                                type="response")
  
  New_Pro3[,i]   <- predict.gbm(p_brt, newdat_Chl,
                                n.trees=p_brt$gbm.call$best.trees, 
                                type="response")
  
  New_Pro4[,i]   <- predict.gbm(p_brt, newdat_SST,
                                n.trees=p_brt$gbm.call$best.trees, 
                                type="response")
  
  New_Pro5[,i]   <- predict.gbm(p_brt, newdat_I0,
                                n.trees=p_brt$gbm.call$best.trees, 
                                type="response")
  
  s_brt  <- gbm.step(data=Train, gbm.x = 1:7, gbm.y = 10,  
                     tree.complexity = 15,
                     learning.rate   = .01,
                     family = "gaussian", silent = T)
  
  New_Syn[,i] <- predict.gbm(s_brt, newdat_lonlat,
                             n.trees=s_brt$gbm.call$best.trees, 
                             type="response")
  
  New_Syn1[,i] <- predict.gbm(s_brt, newdat_Dep,
                              n.trees=s_brt$gbm.call$best.trees, 
                              type="response")
  
  New_Syn2[,i] <- predict.gbm(s_brt, newdat_DOY,
                              n.trees=s_brt$gbm.call$best.trees, 
                              type="response")
  
  New_Syn3[,i] <- predict.gbm(s_brt, newdat_Chl,
                              n.trees=s_brt$gbm.call$best.trees, 
                              type="response") 
  
  New_Syn4[,i] <- predict.gbm(s_brt, newdat_SST,
                              n.trees=s_brt$gbm.call$best.trees, 
                              type="response") 
  
  New_Syn5[,i] <- predict.gbm(s_brt, newdat_I0,
                              n.trees=s_brt$gbm.call$best.trees, 
                              type="response") 
  
  e_brt  <- gbm.step(data=Train, gbm.x = 1:7, gbm.y = 11,  
                     tree.complexity = 15,
                     learning.rate   = .01,
                     family = "gaussian", silent = T)
  
  New_Peuk[,i] <- predict.gbm(e_brt, newdat_lonlat,
                              n.trees=e_brt$gbm.call$best.trees, 
                              type="response")
  
  New_Peuk1[,i] <- predict.gbm(e_brt, newdat_Dep,
                               n.trees=e_brt$gbm.call$best.trees, 
                               type="response")
  
  
  New_Peuk2[,i] <- predict.gbm(e_brt, newdat_DOY,
                               n.trees=e_brt$gbm.call$best.trees, 
                               type="response")
  
  New_Peuk3[,i] <- predict.gbm(e_brt, newdat_Chl,
                               n.trees=e_brt$gbm.call$best.trees, 
                               type="response")
  New_Peuk4[,i] <- predict.gbm(e_brt, newdat_SST,
                               n.trees=e_brt$gbm.call$best.trees, 
                               type="response")
  New_Peuk5[,i] <- predict.gbm(e_brt, newdat_I0,
                               n.trees=e_brt$gbm.call$best.trees, 
                               type="response")
}

#Take average
NewChl  <- apply(New_Chl,  1, mean)
NewPro  <- apply(New_Pro,  1, mean)
NewSyn  <- apply(New_Syn,  1, mean)
NewPeuk <- apply(New_Peuk, 1, mean)

NewChl1 <- apply(New_Chl1,  1, mean)
NewPro1 <- apply(New_Pro1,  1, mean)
NewSyn1 <- apply(New_Syn1,  1, mean)
NewPeuk1<- apply(New_Peuk1, 1, mean)
NewChl2 <- apply(New_Chl2,  1, mean)
NewPro2 <- apply(New_Pro2,  1, mean)
NewSyn2 <- apply(New_Syn2,  1, mean)
NewPeuk2<- apply(New_Peuk2, 1, mean)
NewChl3 <- apply(New_Chl3,  1, mean)
NewPro3 <- apply(New_Pro3,  1, mean)
NewSyn3 <- apply(New_Syn3,  1, mean)
NewPeuk3<- apply(New_Peuk3, 1, mean)
NewChl4 <- apply(New_Chl4,  1, mean)
NewPro4 <- apply(New_Pro4,  1, mean)
NewSyn4 <- apply(New_Syn4,  1, mean)
NewPeuk4<- apply(New_Peuk4, 1, mean)
NewChl5 <- apply(New_Chl5,  1, mean)
NewPro5 <- apply(New_Pro5,  1, mean)
NewSyn5 <- apply(New_Syn5,  1, mean)
NewPeuk5<- apply(New_Peuk5, 1, mean)

#calculate SE for the last 5 envr. variables
SEChl1  <- apply(New_Chl1,  1, sd)
SEPro1  <- apply(New_Pro1,  1, sd)
SESyn1  <- apply(New_Syn1,  1, sd)
SEPeuk1 <- apply(New_Peuk1, 1, sd)
SEChl2  <- apply(New_Chl2,  1, sd)
SEPro2  <- apply(New_Pro2,  1, sd)
SESyn2  <- apply(New_Syn2,  1, sd)
SEPeuk2 <- apply(New_Peuk2, 1, sd)
SEChl3  <- apply(New_Chl3,  1, sd)
SEPro3  <- apply(New_Pro3,  1, sd)
SESyn3  <- apply(New_Syn3,  1, sd)
SEPeuk3 <- apply(New_Peuk3, 1, sd)
SEChl4  <- apply(New_Chl4,  1, sd)
SEPro4  <- apply(New_Pro4,  1, sd)
SESyn4  <- apply(New_Syn4,  1, sd)
SEPeuk4 <- apply(New_Peuk4, 1, sd)
SEChl5  <- apply(New_Chl5,  1, sd)
SEPro5  <- apply(New_Pro5,  1, sd)
SESyn5  <- apply(New_Syn5,  1, sd)
SEPeuk5 <- apply(New_Peuk5, 1, sd)

#plot graph using image2D
jet.colors   <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")) 
NewChl <- matrix(NewChl,
                 nr = length(lonx),
                 nc = length(latx),
                 byrow = F)

NewPro <- matrix(NewPro,
                 nr = length(lonx),
                 nc = length(latx),
                 byrow = F)

NewSyn <- matrix(NewSyn,
                 nr = length(lonx),
                 nc = length(latx),
                 byrow = F)

NewPeuk <- matrix(NewPeuk,
                  nr = length(lonx),
                  nc = length(latx),
                  byrow = F)
#apply mask
source('Calc_bath.R')
Bot_Depth <- get_bath(newdat_lonlat)
Bot_Depth <- matrix(Bot_Depth,
                    nr = length(lonx),
                    nc = length(latx),
                    byrow = F)

NewChl[Bot_Depth >= 0]  <- NA  #Apply mask
NewPro[Bot_Depth >= 0]  <- NA
NewSyn[Bot_Depth >= 0]  <- NA
NewPeuk[Bot_Depth >= 0] <- NA

pdf("Fig4_ind_effectWinter25June2020.pdf", width=8, height=8)
op   <- par(font.lab=1,
            family ="serif",
            mar    = c(4,4,2,3),
            mgp    = c(2,1,0),
            mfrow  = c(3,3),
            cex.axis=1.2, cex.lab=1.2)

#Plot the partial effect of Depth
x2   <- c(10, 10^2, 10^3, 10^4, 10^5, 10^6)
X2   <- c(10, expression(paste(10^2)), expression(paste(10^3)), expression(paste(10^4)), expression(paste(10^5)), expression(paste(10^6)))
y2   <- log(x2)

YL <- NewPro1 - 1.96*SEPro1
YU <- NewPro1 + 1.96*SEPro1
plot(NewPro1, depx, type = 'l', lwd = 1.5,
     xlim= range(y2),
     xaxt='n',
     xlab= expression(paste('Abundance (cells ' * mL^-1 * ')')),
     ylab= 'Depth (m)')
mtext('A', side = 3, adj = -0.1, cex=1.4)

#add 95% CI
lines(YL, depx, lty = 3, col = 1)
lines(YU, depx, lty = 3, col = 1)
axis(1, at= y2, label = X2)

YL <- NewSyn1 - 1.96*SESyn1
YU <- NewSyn1 + 1.96*SESyn1
lines(NewSyn1, depx, col=2, lwd=1.5)

#add 95% CI
lines(YL, depx, lty = 3, col = 2)
lines(YU, depx, lty = 3, col = 2)
axis(1, at= y2, label = X2)

YL <- NewPeuk1 - 1.96*SEPeuk1
YU <- NewPeuk1 + 1.96*SEPeuk1
lines(NewPeuk1, depx, col = 3, lwd = 1.5)

#add 95% CI
lines(YL, depx, lty = 3, col = 3)
lines(YU, depx, lty = 3, col = 3)

#Compute 95% CI envelopes
YL <- NewChl1 - 1.96*SEChl1
YU <- NewChl1 + 1.96*SEChl1
lines(NewChl1+6, depx, col=4, lwd = 1.5)

#add 95% CI
lines(YL+6, depx, lty = 3, col = 4)
lines(YU+6, depx, lty = 3, col = 4)

#Change scale
X1   <- seq(0.1,1,by=.2)
y1   <- log(X1)
axis(3, at= y1+6, label = X1)
mtext(side=3, bquote('Chl (mg '*m^-3*')'), adj = .8, line=.7, cex=.8)
#Add legend
legend('bottomright', 
       legend=c('Pro','Syn','Peuk','Chl'),
       col=1:4, lwd=1.5, lty=1)

#Plot the partial effect of DOY
YL <- NewPro2 - 1.96*SEPro2
YU <- NewPro2 + 1.96*SEPro2
x2   <- c(10^3, 10^4, 10^5, 10^6)
X2   <- c(expression(paste(10^3)), expression(paste(10^4)), expression(paste(10^5)), expression(paste(10^6)))
y2   <- log(x2)
plot(doy, NewPro2, type = 'l', lwd = 1.5,
     ylim= range(y2),
     xaxt='n',
     yaxt='n',
     ylab= expression(paste('Abundance (cells ' * mL^-1 * ')')),
     xlab= 'Month')
lines(doy, YL, lty = 3)
lines(doy, YU, lty = 3)
mtext('B', side = 3, adj = -0.1, cex=1.4)
DOY_Label <- c(15, 75, 135, 195, 255, 315)
axis(1, at= DOY_Label, 
     label = c('Jan', 'Mar', 'May', 'Jul', 'Sep', 'Nov'))
axis(2, at= y2, label = X2)

YL <- NewSyn2 - 1.96*SESyn2
YU <- NewSyn2 + 1.96*SESyn2
lines(doy, NewSyn2, lwd = 1.5, col=2)
lines(doy, YL, col = 2, lty = 3)
lines(doy, YU, col = 2, lty = 3)

YL <- NewPeuk2 - 1.96*SEPeuk2
YU <- NewPeuk2 + 1.96*SEPeuk2
lines(doy, NewPeuk2, lwd = 1.5, col = 3)
lines(doy, YL, col = 3, lty = 3)
lines(doy, YU, col = 3, lty = 3)

#Add Chl
YL <- NewChl2 - 1.96*SEChl2
YU <- NewChl2 + 1.96*SEChl2
lines(doy, NewChl2 + 13, col = 4, lwd = 1.5)
lines(doy, YL+13, lty = 3, col = 4)
lines(doy, YU+13, lty = 3, col = 4)

#Add Chl axis on the right
X1   <- c(0.1, .2, .4, 0.8)
axis(4, at= log(X1) + 13, label = X1)
mtext(expression(paste('Chl (mg '*m^-3*')')),
      side=4,line=.8, cex=.8, adj=.1) 

#Plot the partial effect of Chl
x2   <- c(10^3, 10^4, 10^5, 10^6)
X2   <- c(expression(paste(10^3)), expression(paste(10^4)), expression(paste(10^5)), expression(paste(10^6)))
y2   <- log(x2)
#End of panel F

YL <- NewPro3 - 1.96*SEPro3
YU <- NewPro3 + 1.96*SEPro3
plot(newdat_Chl$logChl0, NewPro3, type='l',
     ylim= range(y2),
     lwd = 1.5,
     xaxt='n',
     yaxt='n',
     xlab= expression(paste('Surface Chl (mg '*m^-3*')')),
     ylab= expression(paste('Abundance (cells '*mL^-1*')')))
axis(2, at= y2, label = X2)

X1   <- c(0.01, 0.1, 1)
axis(1, at= log(X1), label = X1)

lines(newdat_Chl$logChl0, YL, lty = 3)
lines(newdat_Chl$logChl0, YU, lty = 3)
mtext('C', side = 3, adj = -0.1, cex=1.4)

YL <- NewSyn3 - 1.96*SESyn3
YU <- NewSyn3 + 1.96*SESyn3
lines(newdat_Chl$logChl0, NewSyn3,  col=2,lwd = 1.5)
lines(newdat_Chl$logChl0, YL, col = 2, lty = 3)
lines(newdat_Chl$logChl0, YU, col = 2, lty = 3)

YL <- NewPeuk3 - 1.96*SEPeuk3
YU <- NewPeuk3 + 1.96*SEPeuk3
lines(newdat_Chl$logChl0, NewPeuk3, col=3, lwd = 1.5)
lines(newdat_Chl$logChl0, YL, col = 3, lty = 3)
lines(newdat_Chl$logChl0, YU, col = 3, lty = 3)

YL <- NewChl3 - 1.96*SEChl3
YU <- NewChl3 + 1.96*SEChl3
lines(newdat_Chl$logChl0, NewChl3+13.3, col=4,lwd = 1.5)
lines(newdat_Chl$logChl0, YL+13.3, col = 4, lty = 3)
lines(newdat_Chl$logChl0, YU+13.3, col = 4, lty = 3)
#Add Chl axis on the right
X1   <- c(0.1, .2, .4, 0.8)
axis(4, at= log(X1) + 13.3, label = X1)
mtext(expression(paste('Chl (mg '*m^-3*')')),
      side=4,line=.5, adj=.1, cex=.8) 

#Plot the partial effect of SST
YL <- NewPro4 - 1.96*SEPro4
YU <- NewPro4 + 1.96*SEPro4
plot(newdat_SST$T0, NewPro4, lwd=1.5, type = 'l',
     ylim=c(min(YL), max(YU)),
     yaxt='n',
     xlab= 'Sea surface temperature (ºC)',
     ylab= expression(paste('Abundance (cells '*mL^-1*')')))
axis(2, at= y2, label = X2)
lines(newdat_SST$T0, YL, col = 1, lty = 3)
lines(newdat_SST$T0, YU, col = 1, lty = 3)
mtext('D', side = 3, adj = -0.1, cex=1.4)

YL <- NewSyn4 - 1.96*SESyn4
YU <- NewSyn4 + 1.96*SESyn4
lines(newdat_SST$T0, NewSyn4, col=2, lwd=1.5)
lines(newdat_SST$T0, YL, col = 2, lty = 3)
lines(newdat_SST$T0, YU, col = 2, lty = 3)

YL <- NewPeuk4 - 1.96*SEPeuk4
YU <- NewPeuk4 + 1.96*SEPeuk4
lines(newdat_SST$T0, NewPeuk4, lwd = 1.5,col=3)
lines(newdat_SST$T0, YL, col = 3, lty = 3)
lines(newdat_SST$T0, YU, col = 3, lty = 3)

YL <- NewChl4 - 1.96*SEChl4
YU <- NewChl4 + 1.96*SEChl4
lines(newdat_SST$T0, NewChl4+6, lwd = 1.5, col = 4)
lines(newdat_SST$T0, YL+6, col = 4, lty = 3)
lines(newdat_SST$T0, YU+6, col = 4, lty = 3)

#Add Chl axis on the right
X1   <- c(0.1, .2, .4, 0.8)
axis(4, at= log(X1) + 6, label = X1)
mtext(expression(paste('Chl (mg '*m^-3*')')),
      side=4,line=.6, adj = .8, cex=.8) 

#Partial effect of surface PAR
YL <- NewPro5 - 1.96*SEPro5
YU <- NewPro5 + 1.96*SEPro5

plot(newdat_I0$I0, NewPro5, lwd=1.5, type = 'l',
     ylim=range(y2),
     yaxt='n',
     xlab= expression(paste(PAR[sat]*' (E '*m^-2*' '*d^-1*')')),
     ylab= expression(paste('Abundance (cells '*mL^-1*')')))
axis(2, at= y2, label = X2)
mtext('F', side = 3, adj = -0.1, cex=1.4)
lines(newdat_I0$I0, YL, col = 1, lty = 3)
lines(newdat_I0$I0, YU, col = 1, lty = 3)

YL <- NewSyn5 - 1.96*SESyn5
YU <- NewSyn5 + 1.96*SESyn5
lines(newdat_I0$I0, NewSyn5, lwd = 1.5,col=2)
lines(newdat_I0$I0, YL, col = 2, lty = 3)
lines(newdat_I0$I0, YU, col = 2, lty = 3)

YL <- NewPeuk5 - 1.96*SEPeuk5
YU <- NewPeuk5 + 1.96*SEPeuk5
lines(newdat_I0$I0, NewPeuk5, lwd=1.5, col=3)
lines(newdat_I0$I0, YL, col = 3, lty = 3)
lines(newdat_I0$I0, YU, col = 3, lty = 3)

YL <- NewChl5 - 1.96*SEChl5
YU <- NewChl5 + 1.96*SEChl5
lines(newdat_I0$I0, NewChl5+14, lwd=1.5,type='l', col=4)
lines(newdat_I0$I0, YL + 14, col = 4, lty = 3)
lines(newdat_I0$I0, YU + 14, col = 4, lty = 3)

#Add Chl axis on the right
X1   <- seq(0.1,1,by=.2)
y1   <- log(X1)
axis(4, at= log(X1) + 14, label = X1)
mtext(expression(paste('Chl (mg '*m^-3*')')),
      side=4,line=.6, adj = .1, cex=.8)

image2D(exp(NewChl), lonx, latx, 
        col  = jet.colors(20),
        cex.lab = 1.2, 
        cex.axis = 1.2, 
        xlab = 'Longitude (ºE)', 
        ylab = "Latitude (ºN)")
mtext(bquote('G) Chl (mg '*m^-3*')'), side = 3, adj = 0)

image2D(exp(NewPro)/1e4, lonx, latx, 
        col  = jet.colors(20),
        cex.lab = 1.2, 
        cex.axis = 1.2, 
        xlab = 'Longitude (ºE)', 
        ylab = "Latitude (ºN)")
mtext(bquote('H) Pro ('*10^4*' cells '*mL^-1*')'), 
      side = 3, adj = 0)

image2D(exp(NewSyn)/1e4, lonx, latx, 
        col  = jet.colors(20),
        cex.lab = 1.2, 
        cex.axis = 1.2, 
        xlab = 'Longitude (ºE)', 
        ylab = "Latitude (ºN)")
mtext(bquote('I) Syn ('*10^4*' cells '*mL^-1*')'), 
      side = 3, adj = 0)

image2D(exp(NewPeuk)/1e4, lonx, latx, 
        col  = jet.colors(20),
        cex.lab = 1.2, 
        cex.axis = 1.2, 
        xlab = 'Longitude (ºE)', 
        ylab = "Latitude (ºN)")
mtext(bquote('J) Peuk ('*10^4*' cells '*mL^-1*')'), 
      side = 3, adj = 0)
dev.off()
