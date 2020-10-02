#Predict picophytoplankton over surface of the SCS
#Obtain SST for seasonal climatology
setwd("~/OneDrive/SCS_pico_niche")
source('Calc_bath.R')

#MODIS Aqua files
Aquadir <- 'MODISAqua/SeasonClim/'

#Read seasonal climatology data from MODIS-Aqua nc files
File    <- 'A20021722019263.L3m_SCSU_CHL_chlor_a_9km.nc'
File    <- paste0(Aquadir, 'Chl/', File)

#Get grid
LONe <- ncread(File, 'lon')
LATe <- ncread(File, 'lat')

tDOY <- function(x) cos(x/365 * 2*pi) #Transform DOY

#Real lon range: 109-121
#Real lat range: 11-23
klon <- which(LONe >= 110 & LONe <= 120)
LONR <- LONe[klon]
klat <- which(LATe >= 16  & LATe <= 23)
LATR <- LATe[klat]


#Obtain seasonal climatology of SST
surface <- expand.grid(LONR, LATR)
colnames(surface) <- c('lon','lat')

#Get Chl, PAR and SST for each season
surface$Chl_SP <- NA
surface$Chl_SU <- NA
surface$Chl_AU <- NA
surface$Chl_WI <- NA

surface$PAR_SP <- NA
surface$PAR_SU <- NA
surface$PAR_AU <- NA
surface$PAR_WI <- NA

surface$SST_SP <- NA
surface$SST_SU <- NA
surface$SST_AU <- NA
surface$SST_WI <- NA

#Write a loop to find the seasonal data
seasons  <- c('SP', 'SU', 'AU', 'WI')
varnames <- c('chlor_a','par','sst')
folders  <- c('Chl', 'PAR', 'SST')
DOYs     <- c(125, 218, 309, 35)
Fullseason  <- c('Spring', 'Summer', 'Autumn', 'Winter')
for (j in 1:length(varnames)) {
  files <- list.files(paste0(Aquadir, folders[j], '/'))
  for (i in 1:length(seasons)) {
    kf   <- grepl(seasons[i], files, fixed = TRUE)
    file <- paste0(Aquadir, folders[j], '/', files[kf])
    dat  <- ncread(file, varnames[j], 
                   start = c(klon[1], klat[1]), 
                   count = c(length(klon), length(klat) ) )
    #Test map
    #Reverse the matrix
    # dat <- dat[, ncol(dat):1]
    # dat[dat > 1] <- 1
    # image2D(
    #   dat,
    #   LONR,
    #   rev(LATR),
    #   col  = jet.colors(25),
    #   cex.lab  = 1,
    #   cex.axis = 1,
    #   xlab = '',
    #   ylab = ''
    # )
    
    #Copy to surface matrix
    colname <- paste0(folders[j],'_',seasons[i])
    surface[,colname] <- as.vector(dat)
  }
}

#Get bottom depth
surface$Bot_Depth <- get_bath(surface)

#Get brt model
source('PredictSCSPico_example.R')

#-----------------------------------------------------------
jet.colors   <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
pdf("SCSsurface.pdf",width=9,height=7)
par(font.lab = 1,
    family  = "serif",
    oma     = c(4,4,0,2),
    mar     = c(2,2,2,2),
    mgp     = c(2.5,1,0),
    mfrow   = c(4,4)) 

LL <- 0
for (j in 1:length(seasons)){

    colname      <- paste0(folders,'_',seasons[j])
    newx.spr     <- data.frame(lon     =   surface$lon,
                               lat     =   surface$lat,
                               Depth   =   5,
                               DOY     =   tDOY(DOYs[j]),
                               logChl0 =   log(surface[,colname[1]]),
                               T0      =   surface[,colname[3]],
                               I0      =   surface[,colname[2]] )
    
    chl <- exp(predict.gbm(c_brt_full, newx.spr,
                           n.trees=c_brt_full$gbm.call$best.trees, type="response"))
    
    pro <- exp(predict.gbm(p_brt_full, newx.spr, 
                           n.trees=p_brt_full$gbm.call$best.trees, type="response"))
    
    syn <- exp(predict.gbm(s_brt_full, newx.spr, 
                           n.trees=s_brt_full$gbm.call$best.trees, type="response"))
    
    peu <- exp(predict.gbm(e_brt_full, newx.spr, 
                           n.trees=e_brt_full$gbm.call$best.trees, type="response"))
    
    chl[surface$Bot_Depth >=0] = NA
    pro[surface$Bot_Depth >=0] = NA
    syn[surface$Bot_Depth >=0] = NA
    peu[surface$Bot_Depth >=0] = NA	
    
    Chl      <- matrix(chl, nr=length(LONR), nc=length(LATR))
    Pro      <- matrix(pro, nr=length(LONR), nc=length(LATR))/1e4
    Syn      <- matrix(syn, nr=length(LONR), nc=length(LATR))/1e4
    Peu      <- matrix(peu, nr=length(LONR), nc=length(LATR))/1e4
    
    Chl[Chl>1]    <- 1
    Pro[Pro>20]   <- 20
    Syn[Syn>5]    <- 5
    Peu[Peu>0.5]  <- .5000
    
    #Reverse the matrix
    Chl <- Chl[,ncol(Chl):1]
    
    image2D(Chl, LONR, rev(LATR), 
            col  = jet.colors(25),
            cex.lab  = 1, 
            cex.axis = 1, 
            xlab = '', 
            ylab = '')
    LL  <- LL + 1 
    txt <- bquote(.(LETTERS[LL])*") "*.(Fullseason[j]) *" Chl ("*mg*" "*m^-3*")")
    mtext(side=3, adj = 0, cex=.8, txt)
    
    Pro <- Pro[,ncol(Pro):1]
    image2D(Pro, LONR, rev(LATR), 
            col  = jet.colors(25),
            cex.lab  = 1, 
            cex.axis = 1, 
            xlab = '', 
            ylab = '')
    
    LL  <- LL + 1 
    txt <- bquote(.(LETTERS[LL])*") "*.(Fullseason[j]) *" Pro ("*10^4*" cells "*mL^-1*")")
    mtext(side=3, adj = 0, cex=.8, txt)
    
    Syn <- Syn[,ncol(Syn):1]
    image2D(Syn, LONR, rev(LATR), 
            col  = jet.colors(25),
            cex.lab  = 1, 
            cex.axis = 1, 
            xlab = '', 
            ylab = '')
     LL <- LL + 1 
    txt <- bquote(.(LETTERS[LL])*") "*.(Fullseason[j]) *" Syn ("*10^4*" cells "*mL^-1*")")
    mtext(side=3, adj = 0, cex=.8, txt)
    
    Peu <- Peu[,ncol(Peu):1]
    image2D(Peu, LONR, rev(LATR), 
            zlim = c(0, .5),
            col  = jet.colors(25),
            cex.lab  = 1, 
            cex.axis = 1, 
            xlab = '', 
            ylab = '')
     LL <- LL + 1 
    txt <- bquote(.(LETTERS[LL])*") "*.(Fullseason[j]) *" Peuk ("*10^4*" cells "*mL^-1*")")
    mtext(side=3, adj = 0, cex=.8, txt)
}
mtext(side=1, adj = 0.5, line = 0.4, cex=1.2, outer = T, 'Longtitude (ºE)')
mtext(side=2, adj = 0.5, line = 0.4, cex=1.2, outer = T, 'Latitude (ºN)')
dev.off()
