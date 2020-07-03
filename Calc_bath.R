#Calculate bathymetry for a given grid
#Find the bottom depth of each new grid
library(ncdf4)
library(class)
ncread  <- function(file, VAR, start = NA, count = NA){
  if(file.exists(file)){
    nc    <- nc_open(file)
  }else{
    stop(paste(file,'does not exist!'))
  }
  data    <- ncvar_get(nc, VAR, start = start, count = count)
  nc_close(nc)
  return(data)
}

bathfile <- "etopo05.nc"
depth    <- ncread(bathfile, 'ROSE')
lon      <- ncread(bathfile, 'ETOPO05_X')
lat      <- ncread(bathfile, 'ETOPO05_Y')
i        <- which(lon > 109 & lon < 122)
j        <- which(lat > 11 & lat < 24)
depth.scs <- depth[i,j]
Depth.scs <- matrix(NA, nr = length(i)*length(j), nc = 3)
for (u in 1:length(i)){
  for (v in 1:length(j)){
    Depth.scs[(u-1)*length(j)+v,1] = lon[i[u]]
    Depth.scs[(u-1)*length(j)+v,2] = lat[j[v]]
    Depth.scs[(u-1)*length(j)+v,3] = depth.scs[u,v]    
  }
}
colnames(Depth.scs) = c('lon','lat','Depth')
Depth.scs           = as.data.frame(Depth.scs)

#Estimate bathymetry using the nearest neighbourhood method
latlon.train <- Depth.scs[,1:2]
Depth.train  <- Depth.scs[,3]


get_bath <- function(latlon.sim){
  return(as.numeric(as.character(knn(latlon.train, latlon.sim[,c('lon','lat')], Depth.train))))
}
  

