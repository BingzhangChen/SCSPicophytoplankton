setwd("~/OneDrive/SCS_pico_niche")
#First, load disco package
library(dismo)
library(gbm)

#Load fitted models
load('Full_BRT_Model.Rdata')

tDOY <- function(x) cos(x/365 * 2*pi) #Transform DOY

#Create the dataframe for prediction
newdat <- data.frame(
  lon   = 116,
  lat   = 18,
  Depth = 5,
  DOY   = tDOY(15),
  logChl0 = log(0.1),
  T0 = 30,
  I0 = 50
)

c.p    <- predict.gbm(c_brt_full, newdat,
                      n.trees=c_brt_full$gbm.call$best.trees, 
                      type="response")

(exp(c.p)) #Predict Chl a concentration
p.p    <- predict.gbm(p_brt_full, newdat,
                      n.trees=p_brt_full$gbm.call$best.trees, 
                      type="response")
(exp(p.p)) #Predict Prochlorococcus abundance

s.p    <- predict.gbm(s_brt_full, newdat,
                      n.trees=s_brt_full$gbm.call$best.trees, 
                      type="response")
(exp(s.p)) #Predict Synechococccus abundance

e.p    <- predict.gbm(e_brt_full, newdat,
                      n.trees=e_brt_full$gbm.call$best.trees, 
                      type="response")
(exp(e.p)) #Predict Picoeukaryote abundance
