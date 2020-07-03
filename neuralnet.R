setwd("~/OneDrive/SCS_pico_niche")
#Try to predict picophytoplankton abundance from remote sensing variables: SST, surface PAR, surface Chl and RRS variables
rm(list=ls())
library(foreach)
library(neuralnet)
load('np.Rdata')

rmse   <- function(x,y) sd(x-y)

#Neural network
normalize <- function(x) (x-min(x))/(max(x)-min(x))
np1       <- np
for (i in 1:ncol(np1)){
  np1[,i]   <- normalize(np1[,i])
}

#Test the effect of different layers
MAX_N      <- 10
ONE_LAYER  <- 1:MAX_N

#Two hidden layers with 5 possibilities
TWO_LAYERS <- matrix(NA, nr = 5, nc = 2)
TWO_LAYERS[1, ] <- c(2,2)
TWO_LAYERS[2, ] <- c(5,5)
TWO_LAYERS[3, ] <- c(5,10)
TWO_LAYERS[4, ] <- c(10,5)
TWO_LAYERS[5, ] <- c(10,10)

#Combine one and two hidden layers
LAYERS <- vector(mode = 'list', length = MAX_N + nrow(TWO_LAYERS))
for (i in 1:MAX_N){
  LAYERS[[i]] <- ONE_LAYER[i]
}
for (i in (1 + MAX_N):(nrow(TWO_LAYERS) + MAX_N)){
  LAYERS[[i]] <- TWO_LAYERS[i-MAX_N,]
}

#Construct two matrices storing mean and sd of R2 and RMSE of each hidden layer
MEANR2 <- matrix(NA, nr = length(LAYERS), nc = 12)
SDR2   <- MEANR2
v      <- c("lon","lat","Depth","DOY","logChl0","T0","I0")
set.seed(1)
STEPMAX <- 1e8
for (k in 1:length(LAYERS)){
  nn.result <- foreach(i = 1:10, .combine = 'rbind') %do% {
    R2     <- numeric(4)
    RMSEv  <- R2
    MB     <- R2  #Mean bias
    x      <- sample(rownames(np), 0.5*nrow(np))
    y      <- setdiff(rownames(np),x)         #Find the elements of rownames(NP1) which did not belong to x
    Train  <- np1[x,]   #Data for training
    Test   <- np1[y,]
    try(c_nn   <- neuralnet(logChl ~ lon+lat+Depth+DOY+logChl0+T0+I0,
                        data    = Train,hidden       =LAYERS[[k]],
                        stepmax = STEPMAX, linear.output=FALSE))
    
    if (!is.null(c_nn)){
      c.p    <- compute(c_nn,Test[,v])$net.result
      c.p    <- c.p*(max(np$logChl)-min(np$logChl))+min(np$logChl)
      Test$logChl1  <- Test$logChl*(max(np$logChl)-min(np$logChl))+min(np$logChl)
      R2[1]    <- cor(c.p,Test$logChl1)**2
      RMSEv[1] <- rmse(c.p,Test$logChl1)
      MB[1]    <- mean(c.p-Test$logChl1)
    }else{
      R2[1]    <- NA
      RMSEv[1] <- NA
      MB[1]    <- NA
    }
    
    try(p_nn   <- neuralnet(logPro~lon+lat+Depth+DOY+logChl0+T0+I0,
                            data=Train,hidden=LAYERS[[k]],
                            stepmax = STEPMAX,linear.output=FALSE))
    if (!is.null(p_nn)){
      p.p    <- compute(p_nn,Test[,v])$net.result
      p.p    <- p.p*(max(np$logPro)-min(np$logPro))+min(np$logPro)
      Test$logPro1  <- Test$logPro*(max(np$logPro)-min(np$logPro))+min(np$logPro)
      R2[2]    <- cor(p.p,Test$logPro1)**2
      RMSEv[2] <- rmse(p.p,Test$logPro1)
      MB[2]    <- mean(p.p-Test$logPro1)
    }else{
      R2[2]    <- NA
      RMSEv[2] <- NA
      MB[2]    <- NA
    }
    
    try(s_nn   <- neuralnet(logSyn~lon+lat+Depth+DOY+logChl0+T0+I0,
                            data=Train,hidden=LAYERS[[k]],
                            stepmax = STEPMAX,linear.output=FALSE))
    if (!is.null(p_nn)){
      s.p    <- compute(s_nn,Test[,v])$net.result
      s.p    <- s.p*(max(np$logSyn)-min(np$logSyn))+min(np$logSyn)
      Test$logSyn1  <- Test$logSyn*(max(np$logSyn)-min(np$logSyn))+min(np$logSyn)
      R2[3]  <- cor(s.p,Test$logSyn1)**2
      RMSEv[3] <- rmse(s.p,Test$logSyn1)
      MB[3]    <- mean(s.p-Test$logSyn1)
    }else{
      R2[3]    <- NA
      RMSEv[3] <- NA
      MB[3]    <- NA
    }
    
    try(e_nn   <- neuralnet(logPeuk~lon+lat+Depth+DOY+logChl0+T0+I0,
                            data=Train,hidden=LAYERS[[k]],
                            stepmax = STEPMAX,linear.output=FALSE))
    if (!is.null(e_nn)){
      e.p    <- compute(e_nn,Test[,v])$net.result
      e.p    <- e.p*(max(np$logPeuk)-min(np$logPeuk))+min(np$logPeuk)
      Test$logPeuk1 <- Test$logPeuk*(max(np$logPeuk)-min(np$logPeuk))+min(np$logPeuk)
      R2[4]    <- cor(e.p,Test$logPeuk1)**2
      RMSEv[4] <- rmse(e.p,Test$logPeuk1)
      MB[4]    <- mean(e.p - Test$logPeuk1)
    }else{
      R2[4]    <- NA
      RMSEv[4] <- NA
      MB[4]    <- NA
    }
    c(R2, RMSEv, MB)
  }
  MEANR2[k,] <- apply(nn.result,2,function(x)mean(x,na.rm=T))
  SDR2[k,]   <- apply(nn.result,2,function(x)sd(x,na.rm=T))
  save(LAYERS, MEANR2, SDR2, file = 'neuralnet_fit.Rdata')
  print(paste('Neuralnet completed for layer',LAYERS[[k]]))
}

#plot figures
pdf("FigS2_neuralnet_optim.pdf", width=8, height=8)
op   <- par(font.lab=1,
            family ="serif",
            mar    = c(4,4,2,2),
            mgp    = c(2,1,0),
            mfrow  = c(4,1),
            cex.lab= 1.2,
            cex.axis=1.2)

txts <- c('Chl', 'Pro', 'Syn', 'Peuk')
txts <- paste0(letters[1:4],') ',txts)

Xaxis <- 1:length(LAYERS)

for (j in 1:4){
  #Chl fitting
  plot(Xaxis, MEANR2[1:length(LAYERS), j], 
       type = 'b',
       ylim = c(0,1),
       xaxt = 'n',
       xlab = 'Number of neurons',
       ylab = bquote(R^2))
  
  #Add x axis labels
  y2 <- Xaxis
  Y1 <- c(ONE_LAYER, '2*2', '5*5','5*10','10*5','10*10')
  axis(1, at = y2, label = Y1)
  #Add CI
  segments(Xaxis, MEANR2[1:length(LAYERS), j] - 2*SDR2[1:length(LAYERS), j],
           Xaxis, MEANR2[1:length(LAYERS), j] + 2*SDR2[1:length(LAYERS), j])
  
  mtext(txts[j], side=3, adj = 0)
}
dev.off()
