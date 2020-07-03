setwd("~/OneDrive/SCS_pico_niche")
#Calculate the relative importance of each envr. factor
library(dismo)
library(gbm)
load('np.Rdata')
set.seed(1)

#Use tree complexity = 15 and learning rate = 0.01 to examine individual effects
NVAR <- 7
NPop <- 4
varnames <- names(np[,1:7])
Varnames <- c('Lon', 'Lat', 'Depth', 'DOY', 'SSChl', 'SST', 
              expression(paste(PAR[sat])))

Popnames <- c('Chl', 'Pro', 'Syn', 'Peuk')
Nrep     <- 10
imp      <- array(NA, dim = c(NVAR, NPop, Nrep))

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
  
  r1 <- c_brt$contributions
  
  for (k in 1:NVAR){
    imp[k,1,i] <- r1[r1$var == varnames[k], 2]
  }
  
  p_brt  <- gbm.step(data=Train, gbm.x = 1:7, gbm.y = 9,  
                     tree.complexity = 15,
                     learning.rate   = .01,
                     family = "gaussian", silent = T)
  
  r1 <- p_brt$contributions
  
  for (k in 1:NVAR){
    imp[k,2,i] <- r1[r1$var == varnames[k], 2]
  }
  
  s_brt  <- gbm.step(data=Train, gbm.x = 1:7, gbm.y = 10,  
                     tree.complexity = 15,
                     learning.rate   = .01,
                     family = "gaussian", silent = T)
  
  r1 <- s_brt$contributions
  for (k in 1:NVAR){
    imp[k,3,i] <- r1[r1$var == varnames[k], 2]
  }
  
  e_brt  <- gbm.step(data=Train, gbm.x = 1:7, gbm.y = 11,  
                     tree.complexity = 15,
                     learning.rate   = .01,
                     family = "gaussian", silent = T)
  r1 <- e_brt$contributions
  for (k in 1:NVAR){
    imp[k,4,i] <- r1[r1$var == varnames[k], 2]
  }
}

#Taking average
Medianimp <- apply(imp, c(1,2), function(x)quantile(x,probs=.5))
   Lowimp <- apply(imp, c(1,2), function(x)quantile(x,probs=.025))
    Upimp <- apply(imp, c(1,2), function(x)quantile(x,probs=.975))
  
#plot out
pdf("Rel_imp.pdf", width=7, height=6)
  op   <- par(font.lab=1,
              family ="serif",
              mar    = c(2,2,2,1),
              mgp    = c(2.2,1,0),
              oma    = c(3.5,3.5,0,0),
              mfrow  = c(2,2),
              cex.axis=1.2, cex.lab=1.1)
  for (i in 1:NPop){
    plot(Medianimp[,i], pch=16,
         ylim=c(0,80),
         xaxt='n',
         xlab='', 
         ylab='')
    axis(1, at=1:NVAR, labels=Varnames)
    
    #add 95% CI
    segments(1:NVAR, Lowimp[,i],
             1:NVAR,  Upimp[,i])
    
    mtext(paste0(LETTERS[i],') ',Popnames[i]), side=3, adj=0)
  }
  mtext('Environmental factor', 
        side=1, adj=.5, cex=1.2, outer=T, line=.3)
  mtext('Percentage of contribution (%)', 
        side=2, adj=.5, cex=1.2, outer=T, line=.4)
dev.off()
  