setwd("~/OneDrive/SCS_pico_niche")
#Compute monthly climatology in different regions of the South China Sea
#Obtain SST, Chl, and I0 climatology

#Load fitted models
load('Full_BRT_Model.Rdata')

seats   <-  read.csv("seats8d.csv")
seats$sst[seats$sst > 30] <- NA
w <- nchar(as.character(seats$time))
seats$day[w==8] <- str_sub(seats$time[w==8], 1, 1)
seats$day[w==9] <- str_sub(seats$time[w==9], 1, 2)
seats$yr[w==8]  <- as.numeric(paste("20",str_sub(seats$time[w==8], 7, 8),sep=""))
seats$yr[w==9]  <- as.numeric(paste("20",str_sub(seats$time[w==9], 8, 9),sep=""))
seats$mo[w==8]  <- str_sub(seats$time[w==8], 3, 5)
seats$mo[w==9]  <- str_sub(seats$time[w==9], 4, 6)
seats$mo[seats$mo=="Jan"] <- 1
seats$mo[seats$mo=="Feb"] <- 2
seats$mo[seats$mo=="Mar"] <- 3
seats$mo[seats$mo=="Apr"] <- 4
seats$mo[seats$mo=="May"] <- 5
seats$mo[seats$mo=="Jun"] <- 6
seats$mo[seats$mo=="Jul"] <- 7
seats$mo[seats$mo=="Aug"] <- 8
seats$mo[seats$mo=="Sep"] <- 9
seats$mo[seats$mo=="Oct"] <- 10
seats$mo[seats$mo=="Nov"] <- 11
seats$mo[seats$mo=="Dec"] <- 12
seats$mo <- as.numeric(seats$mo)
seats$Date <- as.Date(ISOdate(seats$yr,seats$mo,seats$day))
seats$DOY  <- as.numeric(strftime(seats$Date, format = "%j"))

tDOY         <- function(x) cos(x/365 * 2*pi) #Transform DOY
seats$DOY    <- tDOY(seats$DOY)
seats$logChl <- log(seats$Chl)
seats        <- na.omit(seats)
newx         <- data.frame(lon=rep(116,nrow(seats)*150),
                           lat=rep(18,nrow(seats)*150),
                           Depth=rep(150:1,nrow(seats)),
                           DOY=rep(seats$DOY,each=150),
                           logChl0=rep(seats$logChl,each=150),
                           T0=rep(seats$sst,each=150),
                           I0=rep(seats$par,each=150))
newx$chl    <- exp(predict.gbm(c_brt_full, newx, n.trees=c_brt$gbm.call$best.trees, type="response"))
newx$pro    <- exp(predict.gbm(p_brt_full, newx, n.trees=p_brt$gbm.call$best.trees, type="response"))
newx$syn    <- exp(predict.gbm(s_brt_full, newx, n.trees=s_brt$gbm.call$best.trees, type="response"))
newx$peu    <- exp(predict.gbm(e_brt_full, newx, n.trees=e_brt$gbm.call$best.trees, type="response"))
newx$year   <- rep(seats$yr,each=150)
newx$Day    <- (newx$year-2002)*365+newx$DOY
Chl      <- matrix(newx$chl,nc=150,nr=length(unique(newx$Day)),byrow=TRUE)
Pro      <- matrix(newx$pro,nc=150,nr=length(unique(newx$Day)),byrow=TRUE)
Syn      <- matrix(newx$syn,nc=150,nr=length(unique(newx$Day)),byrow=TRUE)
Peu      <- matrix(newx$peu,nc=150,nr=length(unique(newx$Day)),byrow=TRUE)
jet.colors   <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")) 
x <- 1+365*(0:11)
y <- 2002:2013
pdf("seats.pdf",width=6,height=3)
par(font.lab = 1,
    family  = "serif",
    mar     = c(2,4,.5,2),
    mgp     = c(2.5,1,0),
    mfrow   = c(1,1)) 
filled.contour(unique(newx$Day),-150:-1,Chl,col = jet.colors(15),
               xlab = "", ylab = "Depth (m)",
               plot.axes = {axis(side=1, at=x, labels = y);
                 axis(2)})
filled.contour(unique(newx$Day),-150:-1,Pro,col = jet.colors(21),
               xlab = "", ylab = "Depth (m)",
               plot.axes = {axis(side=1, at=x, labels = y);
                 axis(2)})
filled.contour(unique(newx$Day),-150:-1,Syn,col = jet.colors(19),
               xlab = "", ylab = "Depth (m)",
               plot.axes = {axis(side=1, at=x, labels = y);
                 axis(2)})
filled.contour(unique(newx$Day),-150:-1,Peu,col = jet.colors(21),
               xlab = "", ylab = "Depth (m)",
               plot.axes = {axis(side=1, at=x, labels = y);
                 axis(2)})							
dev.off()
