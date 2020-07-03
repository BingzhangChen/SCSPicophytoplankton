setwd("~/OneDrive/SCS_pico_niche")
load('np.Rdata')
require(ncdf4)

#Draw the station map
#Draw contours of bathymetry

#To get the ncread function

#Simple function to read .nc data:
ncread  <- function(file, VAR, start = NA, count = NA){
  if(file.exists(file)){
    nc    <- nc_open(file)
  }else{
    stop(paste(file,'does not exist!'))
  }
  data    <- ncvar_get(nc, VAR, start = start, count = count)
  nc_close(nc)
  return(data)
}

bathfile <- "etopo05.nc"
depth    <- ncread(bathfile, 'ROSE')
lon      <- ncread(bathfile, 'ETOPO05_X')
lat      <- ncread(bathfile, 'ETOPO05_Y')
i        <- which(lon > min(np$lon) & lon < max(np$lon))
j        <- which(lat > min(np$lat) & lat < max(np$lat))
depth.scs    <- depth[i,j]
depth.scs[depth.scs > 0] <- NA
jet.colors   <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")) 

pdf("Fig1_station1.pdf", width=6, height=6.25)
op   <- par(font.lab=1,
            family ="serif",
            mar    = c(4,4,2,2),
            mgp    = c(2,1,0),
            mfrow  = c(1,1))

ZLEVEL <- c(-5000, -4000, -3000, -2000, -1000, -500, -100, 0)
filled.contour(lon[i], lat[j], depth.scs,
               levels = ZLEVEL,
               col = jet.colors(length(ZLEVEL)-1),
               xlab = "Longitude (°E)", ylab = "Latitude (°N)",
               plot.axes = {axis(1);
                 axis(2);
                 points(np$lon, np$lat, 
                        pch=18, cex=.5, col='white');
                 points(116, 18, cex=1.2, pch=17, col='yellow');
                 text(116, 17.57, 'SEATS', col='yellow')},  #Add SEATS
               key.axes = axis(4, ZLEVEL))

text(109,23.2, pos=4,"Mainland, China")
text(117,15,   pos=4,"Philippines", srt = 90)		
mtext('Depth (m)', side = 3, adj = 0.9)
dev.off()


#Add histograms of different variables
pdf("Fig2_hist.pdf", width=7, height=7)
op   <- par(font.lab=1,
            family ="serif",
            cex.axis = 1.2,
            cex.lab  = 1.2,
            oma    = c(4,4,0,0),
            mar    = c(3.5,4,1.5,.2),
            mgp    = c(2.2,1,0),
            mfrow  = c(3,3))

hist(np$Depth, 
     xlab = 'Sampling depth (m)',
     ylab = '',
     main='')
box()
mtext('A', side = 3, adj = -0.1, cex = 1.2)

hist(np$DOY,
     xaxt = 'n',
     xlab = 'Sampling date',
     ylab = '',
     main='')
DOY_Label <- c(15, 75, 135, 195, 255, 315)
axis(1, at= DOY_Label, 
     label = c('Jan', 'Mar', 'May', 'Jul', 'Sep', 'Nov'))
box()
mtext('B', side = 3, adj = -0.1, cex = 1.2)

hist(np$logChl0,
     xaxt = 'n',
     xlab = bquote('Surface Chl '*italic(a)*' (mg '*m^-3*')'),
     ylab = '',
     main='')
X1   <- c(0.01, 0.1, 1, 10)
y1   <- log(X1)
axis(1, at= y1, label = X1)
box()
mtext('C', side = 3, adj = -0.1, cex = 1.2)

hist(np$T0,
     xlab = 'Sea surface temperature (ºC)',
     ylab = '',
     main='')
box()
mtext('D', side = 3, adj = -0.1, cex = 1.2)

hist(np$I0,
     xlab = expression(paste(PAR[sat]*' (E '*m^-2*' '*d^-1*')')),
     ylab = '',
     main='')
box()
mtext('E', side = 3, adj = -0.1, cex = 1.2)

x2   <- c(1, 10, 10^2, 10^3, 10^4, 10^5, 10^6)
X2   <- c(0, 10, expression(paste(10^2)), expression(paste(10^3)), expression(paste(10^4)), expression(paste(10^5)), expression(paste(10^6)))
y2   <- log(x2)
hist(np$logPro,
     xaxt = 'n',
     xlab = expression(paste('Pro abundance (cells '*mL^-1*')')),
     ylab = '',
     main='')
axis(1, at= y2, label = X2)
box()
mtext('F', side = 3, adj = -0.1, cex = 1.2)

hist(np$logSyn,
     xaxt = 'n',
     xlab = expression(paste('Syn abundance (cells '*mL^-1*')')),
     ylab = '',
     main='')
axis(1, at= y2, label = X2)
box()
mtext('G', side = 3, adj = -0.1, cex = 1.2)

hist(np$logPeuk,
     xaxt = 'n',
     xlab = expression(paste('Peuk abundance (cells '*mL^-1*')')),
     ylab = '',
     main='')
axis(1, at= y2, label = X2)
box()
mtext('H', side = 3, adj = -0.1, cex = 1.2)

hist(np$logChl,
     xaxt = 'n',
     xlab = bquote('Chl '*italic(a)*' (mg '*m^-3*')'),
     ylab = '',
     main = '')
axis(1, at= y1, label = X1)
box()
mtext('I', side = 3, adj = -0.1, cex = 1.2)

mtext('Number of samples', 
      side = 2, adj = .5, cex = 1.4, outer=T, line = .2)
dev.off()

