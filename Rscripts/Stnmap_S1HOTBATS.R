#a) Plot station map
library(plot3D)
source('~/Working/FlexEFT1D/Rscripts/Interpolate_WOA.R')
source('~/Working/FlexEFT1D/Rscripts/interp.R')
#source('~/Working/FlexEFT1D/Rscripts/MLD_obs.R')
Chl = readnc('chl', sourcefile = ChlFile, ROMS=F)
lon = Chl$lon
lat = Chl$lat
CHL = Chl$data
fLon= lon >= 120 & lon <= 360-1
fLat= lat >= 20  & lat <= 60
lon = lon[fLon]
lat = lat[fLat]

#BATS: 31.6667 N (31°40′N), -64.16667 (64°10′W) 
cff  <- c(which(lon <= 180),  which(lon > 180))
Nlon <- lon[cff]

CHL <- CHL[fLon,fLat,]
CHL <- apply(CHL,c(1,2),mean,na.rm=T)
CHL <- CHL[cff,]

LON= c(145, -158 + 360, -64.16667 + 360) #S1, ALOHA, and BATS
LAT= c(30,  22.75,      31.6667)

#Calculate annual mean Chl:
jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))  #A good color

Stn_names = c('S1', 'ALOHA','BATS')

pdf('~/Working/FlexEFT1D/DRAM/S1HOTBATS_map.pdf',
                          width=8, height=4,paper='a4')
op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(3.5,4,2,2),
             mgp     = c(2.3,1,0),
             cex.lab = 1.2,
             lwd     = 1.5,
             pch     = 16,
             mfrow   =c(1,1),
             cex.axis=1) 


image2D(x=Nlon, y=lat, z=CHL,zlim=c(0,2), col = jet.colors(18),xaxt='n',
         xlab = "Longitude (ºE)", ylab = "Latitude (ºN)")
lon1 = seq(140,360,by=20)
lon2 = lon1
lon2[lon2>180]=lon2[lon2>180]-360
axis(1, at = lon1, labels=lon2)

title = expression(paste('Annual mean Chl (mg '*m^-3*')'))
text(x=LON+3.5,y=LAT+2.5,labels=Stn_names,
     cex=1.5,col='yellow')
mtext(title,line=.3)
points(LON,LAT,pch=17,col=2)
dev.off()
