setwd("~/OneDrive/SCS_pico_niche")
library(foreach)
library(plot3D)
load('np.Rdata')

rmse   <- function(x,y) sd(x-y)

set.seed(1)

#GAM
library(mgcv)

np <- np[,c("lon","lat","Depth","DOY",
            "logChl0","T0","I0","logChl",
                                        "logPro","logSyn","logPeuk")]

#Optimize k for k for s(Depth), and gamma
K     <- c(5, 10, 20, 40)
GAMMA <- c(1.2,1.4,1.6)

MEANR2 <- array(NA, dim=c(12,  length(K), length(GAMMA)))
SDR2   <- MEANR2

#Need to test the optimal learning rate (lr), tree complexity (tc), and total number of trees
for (n in 1:length(K)){
  for (p in 1:length(GAMMA)){
    gam.result <- foreach(i=1:10,.combine='rbind') %do% {
      x      <- sample(rownames(np), 0.5*nrow(np))
      y      <- setdiff(rownames(np),x)         #Find the elements of rownames(NP1) which did not belong to x
      Train  <- np[x,]   #Data for training
      Test   <- np[y,]
      c_gam  <- gam(logChl~te(lon,lat) + 
                      s(Depth, k = K[n]) + 
                      s(DOY, bs="cc", k=K[n]) + 
                      s(logChl0) + 
                      s(T0) + s(I0), data = Train, gamma = GAMMA[p])
      p_gam  <- gam(logPro~te(lon,lat)+s(Depth, k = K[n]) + 
                      s(DOY,bs="cc", k=K[n])+s(logChl0)+s(T0)+s(I0),
                    data = Train, gamma = GAMMA[p])
      
      s_gam  <- gam(logSyn ~ te(lon,lat)+s(Depth, k = K[n]) + 
                      s(DOY,bs="cc", k = K[n])+s(logChl0)+s(T0)+s(I0),
                    data = Train, gamma = GAMMA[p])
      
      e_gam  <- gam(logPeuk ~ te(lon,lat)+s(Depth, k = K[n]) + 
                      s(DOY,bs="cc", k = K[n])+s(logChl0)+s(T0)+s(I0),
                    data = Train, gamma = GAMMA[p])
      
      c.p    <- predict(c_gam,Test)
      p.p    <- predict(p_gam,Test)
      s.p    <- predict(s_gam,Test)
      e.p    <- predict(e_gam,Test)
      c(cor(c.p,Test$logChl)**2,
        cor(p.p,Test$logPro)**2,
        cor(s.p,Test$logSyn)**2,
        cor(e.p,Test$logPeuk)**2,
        rmse(c.p,Test$logChl),
        rmse(p.p,Test$logPro),
        rmse(s.p,Test$logSyn),
        rmse(e.p,Test$logPeuk),
        mean(c.p-Test$logChl),
        mean(p.p-Test$logPro),
        mean(s.p-Test$logSyn),
        mean(e.p-Test$logPeuk))
    }
    MEANR2[,n,p] <- apply(gam.result,2,mean)
    SDR2[,n,p]   <- apply(gam.result,2,sd)	
    save(MEANR2, SDR2, file = 'GAM_optimsetting2.Rdata')
  }
}

#Plot out
load('GAM_optimsetting.Rdata')
pdf("FigS3_GAM_optim.pdf", width=8, height=8)
op   <- par(font.lab=1,
            family ="serif",
            mar    = c(4,4,2,2),
            mgp    = c(2,1,0),
            mfcol  = c(4,3),
            cex.lab= 1.2,
            cex.axis=1.2)

txts <- c('Chl', 'Pro', 'Syn', 'Peuk')
txts <- paste0(letters[1:12],') ',rep(txts,3))

nt <- 0
for (i in 1:length(GAMMA)){
  for (j in 1:4){
    
    nt <- nt + 1
    #Chl fitting
    plot(K, MEANR2[j,,i], 
         type = 'b',
         ylim = c(0.6,0.9),
         xlab = 'k',
         ylab = bquote(R^2))
    
    #Add CI
    segments(K, MEANR2[j,,i] - 2*SDR2[j,,i],
             K, MEANR2[j,,i] + 2*SDR2[j,,i])
    
    mtext(paste0(txts[nt], ', gamma=', 
                 GAMMA[i]), side=3, adj = 0)
  }
}
dev.off()


