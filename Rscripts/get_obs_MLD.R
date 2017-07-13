source('~/Working/FlexEFT1D/Rscripts/get_Deu.R')
source('~/Working/FlexEFT1D/Rscripts/interp.R')
source('~/Working/FlexEFT1D/Rscripts/get_MLD.R')
#setwd('~/Working/FlexEFT1D/Forcing')
#Get HOT N:C data
#dat = read.table('HOT_N2C.dat',header=T)
#dat[,3] = 1/dat[,3]
#colnames(dat)[3] = 'C2N'
#write.table(dat, file = 'HOT_QN.dat', row.names=F)

get_obs_MLD <- function(Stn,Var){
#This function calculates avg. value of one variable within MLD
    #Obtain MLD data:
    MLD  <- get_MLD(Stn)
    #Read observational data
    dat  <- paste0('~/Working/FlexEFT1D/Forcing/',Stn,'_',Var,'.dat')
    dat  <- read.table(dat,header=T)
    if(Var=='QN') dat[,3]=1/dat[,3]
   #Calculate the average data within the MLD:
    DOY_uni = unique(sort(dat$DOY))
    N       = length(DOY_uni)
    dat2    = numeric(N)
    #if (Var == 'NPP') Deu = abs(get_Deu(Stn,'EFTsimple'))
    for (i in 1:N){
       DOY_ = DOY_uni[i]
       MLD_ = MLD[DOY_]
       dat_ = dat[dat$DOY==DOY_,]
       xdep = dat_$Depth
    #   if (Var == 'NPP'){
    #   #Calculate integrated NPP (from surface to 0.1% light depth):
    #      #Get the depth of 0.1% of I0
    #      xff = 0
    #    # Deu_ = Deu[DOY_]
    #      Y   = dat_[,3]
    #    #  if (max(xdep) < Deu_) {
    #    #     xdep = c(xdep,Deu_)
    #    #      Y   = c(Y,0)
    #    #  }
    #      for (k in 1:(nrow(dat_)-1)){
    #          xff=xff+(xdep[k+1]-xdep[k])*(Y[k]+Y[k+1])/2
    #      }
    #      dat2[i]=xff

    #   } else{
       #If not NPP:
            dat2_ = dat_[dat_$Depth <= abs(MLD_),]
          if (nrow(dat2_) < 1) {
             dat2[i]= dat_[which.min(xdep),3]
          } else{
             dat2[i]= mean(dat2_[,3])
          }
    #   }
    }
    return(list(DOY=DOY_uni,dat=dat2))
}

source('~/Working/FlexEFT1D/Rscripts/getdata_station.R')
get_latlon = function(Stn_name){
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
  return(list(Lon=stn_lon,Lat=stn_lat))
}
get_theta <- function(Stn_name){
  stn_lon  =  get_latlon(Stn_name)$Lon
  stn_lat  =  get_latlon(Stn_name)$Lat
  cff      = getdata_station('Chl_C',stn_lon,stn_lat)
  cff$data = cff$data[,-1]
  cff$time = 360/length(cff$time)*(cff$time-.5)
  return(cff)
}

