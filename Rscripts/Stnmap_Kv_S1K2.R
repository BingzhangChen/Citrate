#a) Plot station map
library(plot3D)
source('~/Working/FlexEFT1D/Rscripts/Interpolate_WOA.R')
source('~/Working/FlexEFT1D/Rscripts/interp.R')
source('~/Working/FlexEFT1D/Rscripts/MLD_obs.R')
Chl = readnc('chl', sourcefile = ChlFile, ROMS=F)
lon = Chl$lon
lat = Chl$lat
CHL = Chl$data
fLon= lon >= 120 & lon <= 360-90
fLat= lat >= 20  & lat <= 60
lon = lon[fLon]
lat = lat[fLat]

cff  <- c(which(lon <= 180),  which(lon > 180))
Nlon <- lon[cff]

CHL <- CHL[fLon,fLat,]
CHL <- apply(CHL,c(1,2),mean,na.rm=T)
CHL <- CHL[cff,]

LON= c(145, -158 + 360, 160) #S1, ALOHA, and K2
LAT= c(30,  22.75,      47)

#Calculate annual mean Chl:
jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))  #A good color

Stn_names = c('S1', 'ALOHA','K2')
var       = 'Aks'

plot_forcing = function(Stn_name, var, Label){
   
    if (Stn_name == 'K2' && var == 'Aks'){
       #Use the Aks data from Taketo
       outfile  =  '~/Working/FlexEFT1D/K2/K2_AksT.dat'
       timefile =  '~/Working/FlexEFT1D/K2/K2_Aks_timeT.dat'
    } else {
       #Write out timefile:
       timefile =  paste('~/Working/FlexEFT1D/',
                         Stn_name,'/',Stn_name,'_',var,'_time.dat',sep='')
       outfile  =  paste('~/Working/FlexEFT1D/',
                         Stn_name,'/',Stn_name,'_',var,'.dat',sep='')

    }
    
    time     =  read.table(timefile,header=T)
    data     =  read.table(outfile, header=T)
    time     =  as.numeric(time)
    M        =  length(time)
    depth    =  data[,1]
    data1    =  as.matrix(t(data[,2:ncol(data)]))
    
    if (var == 'Aks'){
       #Calculate MLD based on threshold of Kv (1E-4 m2/s):
       MLD <- numeric(M)
       for (i in 1:M){
           MLD[i] <- interp(x=depth,y=data[,i+1],Xth=1E-4)
       } 
       #Calculate MLD based on observation profiles:
       MLD_ob <- MLD_obs(Stn_name)
       ZLIM   <- c(-5.1,0)
    }else{
       if (Stn_name == 'K2'){
          ZLIM   <- c(0,13) #Temperature
       } else{
          ZLIM   <- c(14,28) #Temperature
       }
    }
    x      = which(depth>=-200)
    depth  = depth[x]
    if (var == 'Aks'){
      data1 = log10(data1[,x])
     title1 = bquote(.(Label) ~ ' Log'[10] 
              ~ italic(' K')[v])
    } else if (!(var %in% c('par','solfe','dust'))){
        data1  <- data1[,x]
        title1 <- paste0(Label,var)
       if (var=='temp'){
           title1 <- bquote(.(Label) ~ ' Temperature (ºC) ')
       }
    }

    if (var == 'par'){
       plot(time,data1[,1]/.4,pch=16,type='b', ylim=c(20,150),
            xlab = 'Month',
            ylab = expression(paste('PAR (W '*m^-2*')'))
       )
       title1= bquote(.(Label) ~ ' PAR ')
    } else if(var == 'solfe'){
       plot(time,data1[,1]*1E12,pch=16,type='b', ylim=c(0,.15),
            xlab = 'Month',
            ylab = expression(paste('Deposited iron (ng '*m^-2*' '*s^-1*')'))
       )
       title1= bquote(.(Label) ~ ' Iron deposition')

    } else{
       image2D(data1, x=as.numeric(time), y=depth,lwd=2, zlim=ZLIM,
                     xlab="Month", ylab="Depth (m)")
       
       if (var == 'Aks'){
          lines(time,MLD,col='tan',lwd=5)
          points(MLD_ob$DOY/30,MLD_ob$MLD,pch=16,col='white')
          axis(1, at = seq(1,11,by=2))
       }
    }
    mtext(title1,adj=0,cex=.8)
}

pdf('~/Working/FlexEFT1D/DRAM/Stn_map_S1_K2_HOT.pdf',
                          width=9, height=9,paper='a4')
op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(3.5,4,2,2),
             mgp     = c(2.3,1,0),
             cex.lab = 1.2,
             lwd     = 1.5,
             pch     = 16,
             cex.axis=1) 

nf <- layout(matrix(c(1,1,1,1:13), 4, 4, byrow = TRUE), respect = TRUE)

image2D(x=Nlon, y=lat, z=CHL,zlim=c(0,2), col = jet.colors(18),xaxt='n',
         xlab = "Longitude (ºE)", ylab = "Latitude (ºN)")
lon1 = seq(140,260,by=20)
lon2 = lon1
lon2[lon2>180]=lon2[lon2>180]-360
axis(1, at = lon1, labels=lon2)

title = expression(paste(' Annual mean Chl (mg '*m^-3*')'))
text(x=LON+2.5,y=LAT+2.5,labels=Stn_names,
     cex=1.5,col='yellow')
mtext(title,line=.3)
points(LON,LAT,pch=17,col=2)

plot_forcing('S1', 'Aks', 'b) S1')
plot_forcing('S1', 'temp','c)')
plot_forcing('S1', 'par', 'd)')
plot_forcing('S1', 'solfe','e)')
plot_forcing('K2', 'Aks', 'f) K2')
plot_forcing('K2', 'temp','g)')
plot_forcing('K2', 'par', 'h)')
plot_forcing('K2', 'solfe','i)')
plot_forcing('HOT','Aks','j) ALOHA')
plot_forcing('HOT','temp','k)')
plot_forcing('HOT','par','l)')
plot_forcing('HOT','solfe','m)')
dev.off()
