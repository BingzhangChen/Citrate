#Estimate the temperature variation of a site
source('~/Working/FlexEFT1D/Rscripts/Interpolate_WOA.R')
require(class)
load('~/ROMS/data/WOA13/WOA13_Temp.Rdata') #Contains the object woa13temp
get_t <- function(latlon){
   #Lonlat needs to be 2-column dataframe
   #And the first column needs to be longitude
latlon[,1]   <- replace(latlon[,1], latlon[,1]<0,
                        latlon[latlon[,1]<0,1] + 360)
   lon       <- woa13temp$lon
   lat       <- woa13temp$lat
   ssttop    <- woa13temp$data[ , , 1 , ]
   sstannual <- apply(ssttop, c(1,2), mean, na.rm=TRUE)              #Get annual mean value
   pts       <- expand.grid(lon, lat)
   Y <- as.numeric(as.character(knn(pts[!is.na(sstannual),], latlon, sstannual[!is.na(sstannual)])))               #Average temperature
   return(Y)
}
