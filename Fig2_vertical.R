#Fig. 2 
setwd("~/OneDrive/SCS_pico_niche")
load('np.Rdata')
library(mgcv)
today   = Sys.Date()
Fig2file = paste0('Fig2_Vertial_',today,'.pdf')
pdf(Fig2file,width = 7, height = 7)
op     <- par(font.lab=1,
       family = "serif",
       mar    = c(3,2,1,0.2),
       mgp    = c(2.2,1,0),
       mfrow  = c(2,2),
       oma    = c(0,4,0,0),
       cex.lab= 1.2,
       cex.axis=1.2)
NP1$Dep_ <- -NP1$Depth
ii       <- 0
for (j in c('logChl','logPro','logSyn','logPeuk')){
   if (j == 'logChl'){
      XLab <-  expression(paste('Chl (mg '*m^-3*')'))
   } else if (j == 'logPro'){
      XLab <-  expression(paste('Pro (cells '*mL^-1*')'))
     
   } else if (j == 'logSyn'){
      XLab <-  expression(paste('Syn (cells '*mL^-1*')'))
     
   } else if (j == 'logPeuk'){
      XLab <-  expression(paste('Peuk (cells '*mL^-1*')'))
   }
   
   if (j == 'logChl'){
      x1   <-  c(0.01, 0.1, 1, 10)
      X1   <- x1
   }else{
      x1   <-  c(10, 10^2, 10^3, 10^4, 10^5, 10^6)
      X1   <-  c(10, expression(paste(10^2)), expression(paste(10^3)), expression(paste(10^4)),
                 expression(paste(10^5)), expression(paste(10^6)))
   }
   ii    <- ii + 1
   plot(NP1[,j], NP1$Dep_, type = 'n',
        xaxt = 'n',
        xlab = XLab, 
        ylab = '')
   
   y1 = log(x1)
   axis(1, at= y1, label = X1)
   
   seasons  <- c('Winter','Spring','Summer','Fall')
   for (i in 1:4){
     k        <- 1 + (i - 1) * 3
     dat      <- NP1[NP1$Month %in% k:(k + 2), ]
     dat$x    <- dat[,j]
     if (j=='logPro') dat <- dat[dat$x > 0,]  #Remove zero Pro
     points(dat[,j], dat$Dep_, cex = .2, pch = 16, col = i)
     Gam      <- gam(x ~ s(Dep_, bs = 'cr',k = 4), data = dat, gamma = 1.4)
     newx     <- data.frame(Dep_ = seq(-150,0,1))
     newy     <- predict(Gam,newx)
     lines(newy,newx$Dep_,col = i, lwd = 2)
   }
   if (ii == 4){
      legend('topleft', legend = seasons,
              col = 1:4, lty = 1, lwd = 1.5, cex = 1)
   }
   mtext(LETTERS[ii],adj = 0)
}

mtext('Depth (m)', side = 2, adj = .5, outer = T)
dev.off()
