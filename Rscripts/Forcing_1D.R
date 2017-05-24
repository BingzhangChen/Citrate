#Read in external parameters
#and compute biological parameters only
setwd('~/Working/FlexEFT1D')
library(plot3D)
library(oce)
source('Rscripts/getdata_station.R')
source('Rscripts/interp.R')
source('Rscripts/Hz.R')
source('Rscripts/MLD_obs.R')

WOA13NFile <- '~/ROMS/Data/WOA13/WOA13_NO3.Rdata'
load(WOA13NFile)

WOA13PFile <- '~/ROMS/Data/WOA13/WOA13_PO4.Rdata'
load(WOA13PFile)

WOA13TFile <- '~/ROMS/Data/WOA13/WOA13_Temp.Rdata'
load(WOA13TFile)

#Compile external fortran file:
#system('ifort -c FlexEFT.f90')
#system('ifort -shared -o FlexEFT.so FlexEFT.o')

Stn_name <- 'HOT'
depth    <- 250 #Total depth (m) of the station

if (Stn_name == 'S1'){
    stn_lon  = 145
    stn_lat  = 30
} else if(Stn_name == 'K2'){
    stn_lon  = 160 
    stn_lat  = 47
} else if(Stn_name == 'HOT'){
    stn_lon  = -158
    stn_lat  = 22.75
}

#For each depth, get the profile at the targeted coordinates
 Temp_data <- woa13temp
 Par_data  <- readnc('par')
 NO3_data  <- woa13no3
 PO4_data  <- woa13po4

 Aks_data  <- readnc('Aks')       #From my ROMS output, unit: m2/s
wSODA_data <- readnc('w_SODA')    #unit: m/s
wROMS_data <- readnc('w_ROMS')    #unit: m/s

#Correct timing for Aks and wROMS
m_per_s   <- 1/(3600*24*30)
time      <- Aks_data$time*m_per_s
time      <- time%%12
  Aks_data$time  <- time
wROMS_data$time  <- time


#Write into data files:
for (var in c('temp','par','Aks','NO3','PO4','wROMS','wSODA','wstr')){

    if (var == 'wstr'){
       taux = getdata_station('taux3',stn_lon,stn_lat)
       tauy = getdata_station('tauy3',stn_lon,stn_lat)
       time = taux$time
       data = sqrt(taux$data**2 + tauy$data**2)
    }else{
       cff  = getdata_station(var,stn_lon,stn_lat)
       time = cff$time 
       data = cff$data 
    }

    if (!(var %in% c('Aks','wROMS', 'par', 'wstr')) ){
       data = data[nrow(data):1,]          #Bottom layer first
    }

    outfile  <-  paste('~/Working/FlexEFT1D/',
                        Stn_name,'/',Stn_name,'_',var,'.dat',sep='')
    write.table(data,file=outfile,  
                  append = F,row.names=FALSE,col.names=TRUE) 
    data      <-  read.table(outfile,header=T)
    #Write out timefile:
    timefile  <-  paste('~/Working/FlexEFT1D/',
                         Stn_name,'/',Stn_name,'_',var,'_time.dat',sep='')
    write.table(time,file=timefile,  
                  append = F,row.names=FALSE,col.names=TRUE)  
    time      <-  read.table(timefile,header=T)
}

plot_forcing <- function(Stn_name,var,useTaketo=FALSE){
    outfile  =  paste('~/Working/FlexEFT1D/',
                      Stn_name,'/',Stn_name,'_',var,'.dat',sep='')
    #Write out timefile:
    timefile =  paste('~/Working/FlexEFT1D/',
                      Stn_name,'/',Stn_name,'_',var,'_time.dat',sep='')

    pdffile  = paste('~/Working/FlexEFT1D/',
                     Stn_name,'/',Stn_name,'_',var,'.pdf',sep='')

    if (useTaketo && var == 'Aks'){
     timefile= paste('~/Working/FlexEFT1D/',Stn_name,'/',
                          Stn_name,'_Aks_timeT.dat',sep='')

     outfile = paste('~/Working/FlexEFT1D/',Stn_name,'/',Stn_name,'_AksT.dat',sep='')
     pdffile = paste('~/Working/FlexEFT1D/',
                     Stn_name,'/',Stn_name,'_',var,'T.pdf',sep='')

    }
    time     =  read.table(timefile,header=T)
    data     =  read.table(outfile, header=T)

       pdf(pdffile,width=5,height=5,paper='a4')
          par(font.lab  = 1,
              family    = "serif",
              mar       = c(4,4,4,4),
              mgp       = c(2.2,1,0))
          time  = as.numeric(time)
          M     = length(time)
          depth = data[,1]
          data1 = as.matrix(t(data[,2:ncol(data)]))

          if (var == 'Aks'){
             #Calculate MLD based on threshold of Kv (1E-4 m2/s):
             MLD <- numeric(M)
             for (i in 1:M){
                 MLD[i] <- interp(x=depth,y=data[,i+1],Xth=1E-4)
             } 
             #Calculate MLD based on observation profiles:
             MLD_ob <- MLD_obs(Stn_name)
          }
          x     = which(depth>=-250)
          depth = depth[x]

          if (var == 'Aks'){
            data1 = log10(data1[,x])
           title1 = bquote('Log'[10] 
                    ~ italic(' K')[v] ~ ' at ' ~
                    .(Stn_name))
        
          } else{
              data1  <- data1[,x]
              title1 <- var
             if (var=='temp'){
                 title1 <- bquote('Temperature (ºC) at ' ~ .(Stn_name))
             }
          }

          image2D(data1, x=as.numeric(time), y=depth,lwd=2,
                      xlab="Month", ylab="Depth (m)")
          mtext(title1,line = 0.3)
          if (var=='Aks') {
              lines(time,MLD,col='tan',lwd=5)
              points(MLD_ob$DOY/30,MLD_ob$MLD,pch=0,cex=1.5,col='white')
          }
          axis(1, at = seq(1,11,by=2))
       dev.off()
}


#Get data from 3D ROMS output for back calculation:
setwd('~/Roms_tools/Run/NPacS') 
source('~/Working/FlexEFT1D/Rscripts/getdata_station.R')
#Get surface PAR for the given station
sPAR <- getdata_station('radsw',stn_lon,stn_lat)
sPAR$data <- sPAR$data*.43

#Calculate vertical profiles of PAR:
#Get mixed layer depth:
HBL <- getdata_station('hbl',stn_lon,stn_lat)

#Get seasonal CHL data:
CHL <- getdata_station('Chl_roms',stn_lon,stn_lat) #Chl largely OK

#Get average size:
LNV <- getdata_station('LNV',stn_lon,stn_lat)

#Get temperature:
Temp <- getdata_station('temp_roms',stn_lon,stn_lat)

#Get iron:
DFe <- getdata_station('DFE',stn_lon,stn_lat)

#Get no3:
NO3 <- getdata_station('NO3_roms',stn_lon,stn_lat)

#Get Hz:
depth <- CHL$data$Depth
depth <- c(depth,0)
#Calculate Z_w:
N     <- length(depth)
Z_w   <- as.double(N)
Z_w[N]   <- 0
Z_w[N-1] <- depth[N-1]*2
Hz       <- as.double(N-1)
Hz[N-1]  <- Z_w[N]-Z_w[N-1]       
for (i in (N-2):1){
    Z_w[i] <- depth[i] - (Z_w[i+1]-depth[i]) 
     Hz[i] <- Z_w[i+1] - Z_w[i]
}

#Get vertical PAR:
source('~/Working/FlexEFT1D/Rscripts/Calculate_PAR.R')
source('~/Working/FlexEFT1D/Rscripts/NPZDcont.R')
M   <- length(sPAR$data)
I0_ <- as.double(sPAR$data[M])
chl <- CHL$data[M]
chl <- chl[,1]
Par <- Cal_PAR(I0_,N-1,Hz,chl)
hbl <- HBL$data[M]
hbl <- hbl[,1]
PARavg <- 0
for (i in (N-1):1){
   if (depth[i] > -hbl){
      PARavg <- PARavg+Par[i]*Hz[i]
   }
}
PARavg = PARavg/abs(hbl)

Pmu <- LNV$data[M]
Pmu <- Pmu[,1]

no3 <- NO3$data[M]
no3 <- no3[,1]

fer <- DFe$data[M]
fer <- fer[,1]

temp <- Temp$data[M]
temp <- temp[,1]

#Calculate growth rate:
N <- 40
cff <- NPZDCONT(Temp_=temp[N], PAR_=PARavg, NO3=no3[N],
         PMU=Pmu[N],
          Fe=fer[N] )


#Get Aks:

Aks <- getdata_station('Aks',stn_lon,stn_lat)
