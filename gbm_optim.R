setwd("~/OneDrive/SCS_pico_niche")
library(foreach)
library(plot3D)
#boosted regression trees
library(dismo)
library(gbm)
load('np.Rdata')

rmse <- function(x,y) sd(x-y)
set.seed(1)

TREES  <- c(2,5,10,15)
LR     <- (1:10)/1e3

MEANR2 <- array(NA, dim=c(12, length(TREES), length(LR)))
SDR2   <- MEANR2

#Need to test the optimal learning rate (lr), tree complexity (tc), and total number of trees
for (k in 1:length(TREES)){
  for (l in 1:length(LR)){ 
    gbm.result <- foreach(i=1:10, .combine='rbind') %do% {
      x      <- sample(rownames(np), 0.5*nrow(np))
      y      <- setdiff(rownames(np),x)         #Find the elements of rownames(NP1) which did not belong to x
      Train  <- np[x,]   #Data for training
      Test   <- np[y,]
      
      #All predictors included
      c_brt  <- gbm.step(data=Train, gbm.x = 1:7, gbm.y = 8,  
                         tree.complexity = TREES[k],
                         learning.rate   = LR[l],
                         max.trees = 20000,
                         family = "gaussian", silent = T)
      
      c.p    <- predict.gbm(c_brt, Test,
                            n.trees=c_brt$gbm.call$best.trees, 
                            type="response")
      
      p_brt  <- gbm.step(data=Train, gbm.x = 1:7, gbm.y = 9,  
                         tree.complexity = TREES[k],
                         learning.rate   = LR[l],
                         max.trees = 20000,
                         family = "gaussian", silent = T)
      
      p.p    <- predict.gbm(p_brt, Test,
                            n.trees=p_brt$gbm.call$best.trees, 
                            type="response")
      
      s_brt  <- gbm.step(data=Train, gbm.x = 1:7, gbm.y = 10,  
                         tree.complexity = TREES[k],
                         learning.rate   = LR[l],
                         max.trees = 20000,
                         family = "gaussian", silent = T)
      s.p    <- predict.gbm(s_brt, Test,
                            n.trees=s_brt$gbm.call$best.trees, 
                            type="response")
      
      e_brt  <- gbm.step(data=Train, gbm.x = 1:7, gbm.y = 11,  
                         tree.complexity = TREES[k],
                         learning.rate   = LR[l],
                         max.trees = 20000,
                         family = "gaussian", silent = T)
      
      e.p    <- predict.gbm(e_brt, Test,
                            n.trees=e_brt$gbm.call$best.trees, 
                            type="response")
      c(cor(c.p,Test$logChl)**2,
        cor(p.p,Test$logPro)**2,
        cor(s.p,Test$logSyn)**2,
        cor(e.p,Test$logPeuk)**2,
        rmse(c.p,Test$logChl),
        rmse(p.p,Test$logPro),
        rmse(s.p,Test$logSyn),
        rmse(e.p,Test$logPeuk),
        mean(c.p - Test$logChl),
        mean(p.p - Test$logPro),
        mean(s.p - Test$logSyn),
        mean(e.p - Test$logPeuk))
    }
    
    MEANR2[,k,l] <- apply(gbm.result,2,mean)
    SDR2[,k,l]   <- apply(gbm.result,2,sd)
    save(np, TREES, LR, MEANR2, SDR2, file = 'BRTs_optimsetting2.Rdata')
    print(paste('Completed for BRT Tree',TREES[k],'with learning rate',LR[l]))
  }
}

#load('BRTs_optimsetting2.Rdata')
#plot out
pdf("FigS5_gbm_optim2.pdf", width=5, height=8)
op   <- par(font.lab=1,
            family ="serif",
            mar    = c(3.5,3.5,1,.5),
            mgp    = c(2,1,0),
            mfcol  = c(4,1),
            cex.lab= 1.2,
            cex.axis=1.2)

txts <- c('Chl', 'Pro', 'Syn', 'Peuk')
txts <- paste0(letters[1:4],') ', txts)

for (j in 1:4){
  for (i in 1:length(TREES)){
    #Chl fitting
    if (i == 1){
      plot(LR, MEANR2[j,i,], 
           ylim = c(0.6,0.9),
           cex.lab = 1.4, 
           cex.axis= 1.4, cex.main = 1.4,
           type = 'b',
           xlab = "Learning rate", 
           ylab = bquote(R^2))
      
      if (j == 1){
        legtxt <- paste0('Tree complexity =', TREES)
        legend('bottomright', legend = legtxt,
               col = 1:4, lty = 1)
      }
    }else{
      points(LR, MEANR2[j,i,], type = 'b', col = i)
    }
    
    #Add CI
    segments(LR, MEANR2[j,i,] - 2*SDR2[j,i,],
             LR, MEANR2[j,i,] + 2*SDR2[j,i,],
             col = i)
    
    
  }
  mtext(txts[j], side=3, adj = 0)
  
}
dev.off()

#Find best parameter combinations
for (j in 1:4){
  which(MEANR2[j,,] == max(MEANR2[j,,]), arr.ind = TRUE)
  which(MEANR2[j+4,,] == min(MEANR2[j+4,,]), arr.ind = TRUE)
}



# Simplify
#c_brt_sim <- gbm.simplify(c_brt)
#c_brt_sim$pred.list[[4]]
