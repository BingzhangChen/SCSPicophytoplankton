setwd("~/OneDrive/SCS_pico_niche")
#------------------------------------------------------------
#This R script reads new data file and generates a new dataframe NP1.
NP           <- read.csv('NewPico.csv')
NP[(NP$Chl==0)&(!is.na(NP$Chl)),]$Chl <- min(NP$Chl[NP$Chl>0], na.rm=T)          #Force zero Chl to 0.01 µg/L

NP$NO3       <- as.numeric(as.character(NP$NO3))
NP$Pro       <- as.numeric(as.character(NP$Pro))       #Change data type
NP$Syn       <- as.numeric(as.character(NP$Syn))
NP$Peuk      <- as.numeric(as.character(NP$Peuk))
NP$HB        <- as.numeric(as.character(NP$HB))
NP$I0        <- NP$I0sat
NP$logChl0   <- log(NP$Chl0)
NP$SSS       <- NA    #Sea surface salinity
kw           <- 0.04  #light attenuation due to seawater
kchl         <- 0.025 #light attenuation due to Chl a
NP1          <- NP[0,]

for (i in unique(NP$ID)){
  A <- NP[NP$ID == i,]
  k <- which(is.na(A$Chl))
  if (length(k) <= nrow(A) - 2){
    A$Chl[k] <-  approx(A$Depth, A$Chl, A$Depth[k])$y    #linear interpolation
  } else{
    A$Chl <- NA
  }
  l <- which(is.na(A$NO3))
  if (length(l) <= nrow(A) - 2){
    A$NO3[l] <-  approx(A$Depth, A$NO3, A$Depth[l])$y
  } else{
    A$NO3 <- NA
  }
  A$Depth1 <-  c(A$Depth[1], A$Depth[-nrow(A)])     #integrate Chl
  A$Chl1   <-  c(A$Chl0[1], A$Chl[-nrow(A)])
  
  #Estimate PAR at each depth
  for (j in 1:nrow(A)){
    A$PAR[j] <- A$I0[j]*exp(-A$Depth[j]*(kw + kchl*sum((A$Depth1[j] - A$Depth1[j])*(A$Chl[j] + A$Chl1[j])/2, na.rm=T)))
  }
  
  #Extract sea surface salinity
  A$SSS <- rep(A[A$Depth <= 5, 'Sal'][1], nrow(A))
  NP1 <- rbind(NP1, A)
}   

NP1$logpar    <- log(NP1$PAR)
NP1$logno3    <- log(NP1$NO3+.01)
NP1$logPro    <- log(NP1$Pro+1)
NP1$logSyn    <- log(NP1$Syn+1)
NP1$logPeuk   <- log(NP1$Peuk+1)
NP1$logChl    <- log(NP1$Chl)

v <- c("lon","lat","Depth","DOY",
            "logChl0","T0","I0","logChl","logPro","logSyn","logPeuk")
np     <- na.omit(NP1[,v])
np$DOY <- cos(np$DOY/365 * 2*pi)
save(np, file = 'np.Rdata')
