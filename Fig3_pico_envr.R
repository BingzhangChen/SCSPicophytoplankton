#Plot picophytoplankton vs. Chl
setwd("~/OneDrive/SCS_pico_niche")
source('Make_newdata.R')

#Use normal plotting functions
env  <- c('Temp', 'logpar', 'logno3', 'logChl')
phy  <- c('logPro', 'logSyn', 'logPeuk')
vars <- expand.grid(env, phy)
vars <- apply(vars, c(1,2), function(x)as.character(x))

env.xlab <- c(expression(paste(Temp*' (ºC)')),
              expression(paste(PAR[z]*' (E '*m^-2*' '*d^-1*')')),
              'Nitrate (µM)',
              expression(paste(Chl*' (mg '*m^-3*')')))

phy.ylab <- c(expression(paste('Pro (cells '*mL^-1*')')), 
              expression(paste('Syn (cells '*mL^-1*')')),
              expression(paste('Peuk (cells '*mL^-1*')')))

env.xaxt <- c('s', rep('n',3))

#Adjust of xlab
XADJ <- c(.08, 0.35, 0.65, 0.93)

#Axis ticks for Y axis
y1   <-  c(1, 10^2, 10^4, 10^5)
Y1   <-  c(0, expression(paste(10^2)), expression(paste(10^4)),
           expression(paste(10^5)))
y2   <-  log(y1)
ii   <- 0
today    <- Sys.Date()
Fig3file <- paste0('Fig3_Env_',today,'.pdf')
pdf(Fig3file, width = 9, height = 7, paper = 'a4')

op     <- par(font.lab=1,
              family = "serif",
              mar    = c(2,2,1.2,1.2),
              mgp    = c(2.2,1,0),
              mfrow  = c(3,4),
              oma    = c(4,4,0,0),
              cex.lab= 1.2,
              cex.axis=1.2)


for (j in 1:length(phy)){
  for (i in 1:length(env)){
    ii <- ii + 1
    plot(NP1[,env[i]], NP1[,phy[j]], pch=16, cex=.5,
         xlab = '',
         ylab = '',
         xaxt = env.xaxt[i],
         yaxt = 'n')
    
    #Add surface data for comparison
    points(NP1[NP1$Depth <= 5, env[i]], NP1[NP1$Depth <= 5, phy[j]], 
           pch=16, cex=.5, col=2)
    
    axis(2, at = y2, label = Y1)
    
    if (i != 1){
      #Axis ticks for X axis
      x1   <-  c(0.01, 0.1, 1, 10)
      x2   <-  log(x1)
      axis(1, at = x2, label = x1)
    }
    mtext(LETTERS[ii], adj = 0)
    
    #Add xlab
    if (j == length(phy)){
      mtext(env.xlab[i], side = 1, line=.5,
            adj = XADJ[i], 
            outer = T)
    }
  }
  mtext(phy.ylab[j], side = 2, line=.3,
        adj = 1-(j-1)/4-j/9, outer = T)
}
legend('bottomright', legend = c('All data', 'Surface'),
       col = 1:2, pch=16, cex = 1.2)
dev.off()


# plots <- function(fname=NP1, x='Temp', y='logPro', id = 1){
#   fname$x  <- fname[,x]
#   fname$y  <- fname[,y]
# 
#   lab <- expand.grid(env.xlab, phy.ylab)
#   lab <- apply(lab, c(1,2), function(x)as.character(x))
#   
#   LO_theme <- theme_bw(base_family = "serif") + theme(
#     panel.grid.major = element_line(size = 0, color="white"),
#     panel.grid.minor = element_line(size = 0, color="white"),
#     panel.border     = element_rect(size = .8, color="black", fill = NA),
#     axis.line        = element_line(size = .8, color = "black"),
#     axis.text        = element_text(size = 5),
#     axis.title       = element_text(size = 8),
#     legend.position  = c(1.5,.25),
#     text             = element_text(size = 6),
#     plot.background  = element_rect(color = "white"),
#     panel.background = element_rect(size=.2, color = "black"),
#     panel.margin     = unit(c(0,0,0,0), "lines")
# )
# 
# return(graph  <-   ggplot(fname, aes(x=x, y=y)) +
#   stat_density2d(aes(fill = ..level..), geom="polygon")+
#   scale_fill_gradientn(colours=rev(rainbow(100, start = 0.1, end=.8)))+  
#   xlab(lab[id,1]) +
#   ylab(lab[id,2]) +
#   geom_smooth(colour="red",size=.6, method='loess')    +
#     
#   geom_point(data = fname[fname$Depth <= 5, ], 
#              aes(x,y), shape = 16, size=.2)            +
#   scale_x_discrete(breaks=c("0","4.605","9.21", "13.82"),
#                    labels=c("0", expression(paste(10^2)), 
#                             expression(paste(10^4)), 
#                             expression(paste(10^6))) ) +
#     
#   annotate("text",  x = min(fname$x,na.rm=T) + .1, 
#                     y = max(fname$y,na.rm=T), 
#            label  = LETTERS[id], 
#            family = "serif", size=3)                   +
#   LO_theme)
# }

#N = 100
#z = array(NA, dim = c(N,2,3))
#for (i in 1:3){
#    x   = rnorm(N)
#    y   = rnorm(N)
#    z[,1,i] = x
#    z[,2,i] = y
#}
#
#lm_eqn <- function(df){
#    m <- lm(y ~ x, df);
#    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
#         list(a = format(coef(m)[1], digits = 2), 
#              b = format(coef(m)[2], digits = 2), 
#             r2 = format(summary(m)$r.squared, digits = 3)))
#    as.character(as.expression(eq));                 
#}
#
#par(mfrow=c(3,1))
#for (i in 1:3){
#   df = data.frame(x = z[,1,i], y = z[,2,i])
#p  <-   ggplot(df, aes(x=x, y=y)) +
#  annotate("text",  x = 0+.1*i, 
#                    y = 0+.1*i, 
#           label  = lm_eqn(df), 
#           family = "serif", size=8,parse=T)
#p
#}
#p1 <- p + 
#annotate('text',x = 25, y = 300, label = lm_eqn(df), parse = TRUE)+
#annotate('text',x = 30, y = 350, label = lm_eqn(df), parse = TRUE)


#allplots <- lapply(1:nrow(vars), function(i) plots(x=vars[i,1], y=vars[i,2],id=i))


 
#    do.call(grid.arrange,
#            c(allplots, nrow = 3, 
#              padding = unit(0, "line")))
# dev.off()


