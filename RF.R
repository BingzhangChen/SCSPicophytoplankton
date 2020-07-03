setwd("~/OneDrive/SCS_pico_niche")
library(foreach)
library(randomForest)
load('np.Rdata')

rmse   <- function(x,y) sd(x-y)
set.seed(1)
np <- np[,c("lon","lat","Depth","DOY",
            "logChl0","T0","I0","logChl",
            "logPro","logSyn","logPeuk")]
#Random forests
#parameters to be optimized
#ntree: number of trees to grow
#mtry: number of variables randomly selected
#nperm: number of permutations
NTs <- c(500, 1000, 2000)
MTs <- c(3,6)

MEANR2 <- array(NA, dim=c(12, length(NTs), 
                          length(MTs)))
SDR2   <- MEANR2

for (m in 1:length(NTs)){
  for (n in 1:length(MTs)){
      rf.result <- foreach(i=1:10,.combine='rbind') %do% {
        x      <- sample(rownames(np), 0.5*nrow(np))
        y      <- setdiff(rownames(np),x)         #Find the elements of rownames(NP1) which did not belong to x
        Train  <- np[x,]   #Data for training
        Test   <- np[y,]
        c_rf   <- randomForest(logChl ~ lon + lat + 
                                 Depth + DOY + 
                                 logChl0 + T0 + I0,
                               ntree=NTs[m],
                               mtry =MTs[n],
                               data =Train)
        p_rf   <- randomForest(logPro ~ lon + lat + 
                                 Depth + DOY + 
                                 logChl0 + T0 + I0,
                               ntree=NTs[m],
                               mtry =MTs[n],
                               data =Train)
        s_rf   <- randomForest(logSyn ~ lon + lat + 
                                 Depth + DOY + 
                                 logChl0 + T0 + I0,
                               ntree=NTs[m],
                               mtry =MTs[n],
                               data =Train)
        e_rf   <- randomForest(logPeuk ~ lon + lat + 
                                 Depth + DOY + 
                                 logChl0 + T0 + I0,
                               ntree=NTs[m],
                               mtry =MTs[n],
                               data =Train)
        
        c.p    <- predict(c_rf,Test)
        p.p    <- predict(p_rf,Test)
        s.p    <- predict(s_rf,Test)
        e.p    <- predict(e_rf,Test)
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
      MEANR2[,m,n]=apply(rf.result,2,mean)
      SDR2[,m,n]  =apply(rf.result,2,sd)	
      save(MEANR2, SDR2, file = 'RF_optimsetting3.Rdata')
  }
}

load('RF_optimsetting3.Rdata')
#plot out
pdf("FigS4_RF_optim.pdf", width=6, height=8)
op   <- par(font.lab=1,
            family ="serif",
            mar    = c(4,4,2,2),
            mgp    = c(2,1,0),
            mfcol  = c(4,2),
            cex.lab= 1.2,
            cex.axis=1.2)

txts <- c('Chl', 'Pro', 'Syn', 'Peuk')
txts <- paste0(letters[1:8],') ',rep(txts,2))
nt   <- 0
for (i in 1:length(MTs)){
  for (j in 1:4){
    
    nt <- nt + 1
    #Chl fitting
    plot(NTs, MEANR2[j,,i], 
         type = 'b',
         ylim = c(0.6,0.9),
         xlab = 'Number of trees',
         ylab = bquote(R^2))
    
    #Add CI
    segments(NTs, MEANR2[j,,i] - 2*SDR2[j,,i],
             NTs, MEANR2[j,,i] + 2*SDR2[j,,i])
    
    mtext(paste0(txts[nt], ', mtry=', MTs[i]), side=3, adj = 0)
  }
}
dev.off()


#Final randomForest model
# rf.chl <-
#   randomForest(logChl ~ lon + lat + Depth + DOY + logChl0 + T0 + I0, data =
#                  np)
# rf.pro <-
#   randomForest(logPro ~ lon + lat + Depth + DOY + logChl0 + T0 + I0, data =
#                  np)
# rf.syn <-
#   randomForest(logSyn ~ lon + lat + Depth + DOY + logChl0 + T0 + I0, data =
#                  np)
# rf.peu <-
#   randomForest(logPeuk ~ lon + lat + Depth + DOY + logChl0 + T0 + I0, data =
#                  np)
# 
# print(rf.chl)
# varImpPlot(rf.chl)
# 
# print(rf.pro)
# imp    <- varImpPlot(rf.pro)
# impvar <- rownames(imp)[order(imp[, 1], decreasing = TRUE)]
# 
# print(rf.syn)
# imp    <- varImpPlot(rf.syn)
# impvar <- rownames(imp)[order(imp[, 1], decreasing = TRUE)]
# 
# print(rf.peu)
# imp    <- varImpPlot(rf.peuk)
# impvar <- rownames(imp)[order(imp[, 1], decreasing = TRUE)]
# 
# pdf("partialPlot_Chl_RF.pdf",
#     width = 12,
#     height = 12)
# op     <- par(
#   font.lab = 1,
#   family = "serif",
#   mar    = c(3, 3, 2, .2),
#   mgp    = c(2, 1, 0),
#   mfrow  = c(3, 3)
# )
# 
# imp    <- varImpPlot(rf.chl)
# for (i in 1:nrow(imp)) {
#   partialPlot(rf.chl,
#               np,
#               rownames(imp)[i],
#               xlab = rownames(imp)[i],
#               main = "")
# }
# dev.off()
# 
# pdf("partialPlot_Pro_RF.pdf",
#     width = 12,
#     height = 12)
# op     <- par(
#   font.lab = 1,
#   family = "serif",
#   mar    = c(3, 3, 2, .2),
#   mgp    = c(2, 1, 0),
#   mfrow  = c(3, 3)
# )
# 
# imp    <- varImpPlot(rf.pro)
# for (i in 1:nrow(imp)) {
#   partialPlot(rf.pro,
#               np,
#               rownames(imp)[i],
#               xlab = rownames(imp)[i],
#               main = "")
# }
# dev.off()
# 
# pdf("partialPlot_Syn_RF.pdf",
#     width = 12,
#     height = 12)
# op     <- par(
#   font.lab = 1,
#   family = "serif",
#   mar    = c(3, 3, 2, .2),
#   mgp    = c(2, 1, 0),
#   mfrow  = c(3, 3)
# )
# 
# imp    <- varImpPlot(rf.syn)
# for (i in 1:nrow(imp)) {
#   partialPlot(rf.syn,
#               np,
#               rownames(imp)[i],
#               xlab = rownames(imp)[i],
#               main = "")
# }
# dev.off()
# 
# pdf("partialPlot_Peuk_RF.pdf",
#     width = 12,
#     height = 12)
# op     <- par(
#   font.lab = 1,
#   family = "serif",
#   mar    = c(3, 3, 2, .2),
#   mgp    = c(2, 1, 0),
#   mfrow  = c(3, 3)
# )
# 
# imp    <- varImpPlot(rf.peu)
# for (i in 1:nrow(imp)) {
#   partialPlot(rf.peu,
#               np,
#               rownames(imp)[i],
#               xlab = rownames(imp)[i],
#               main = "")
# }
# dev.off()
