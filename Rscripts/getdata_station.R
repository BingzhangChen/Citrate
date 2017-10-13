source('~/Working/FlexEFT1D/Rscripts/Interpolate_WOA.R')
getdata_station <- function(varname,lon,lat){
    #First, get the total dataset
    ROMS <- FALSE
    if (varname == 'temp'){
      data_s <- Temp_data
    } else if (varname == 'solfe'){
      data_s=readnc(varname,sourcefile=Dustfile,ROMS=F)
    } else if (varname == 'fer'){
      data_s=readnc(varname,sourcefile=DFeFile,ROMS=F)
    } else if (varname == 'par'){
      data_s <- Par_data
    } else if (varname %in% c('taux3','tauy3')){

      coads_dir='/Users/apple/ROMS/Data/Roms_tools/COADS05/';
      #  Wind stress
      #
      file=paste0(coads_dir,varname,'.nc')
      data_s=readnc(varname,sourcefile=file,ROMS=F)

    } else if (varname %in% c('radsw','hbl','Aks','Chl_roms','LNV','temp_roms','DFE','NO3_roms')){
      avgfile <- '~/Roms_tools/run/NpacS/npacS_avg11.nc'
      data_s  <- readnc(varname, sourcefile = avgfile)
      ROMS    <- TRUE
      #readnc <- function(variable, sourcefile = ROMSFile, ROMS=T, nameit='npacS'){
    } else if (varname == 'NO3'){
      WOA13NFile <- '~/ROMS/Data/WOA13/WOA13_NO3.Rdata'
      load(WOA13NFile)
      data_s <- woa13no3
    } else if (varname == 'PO4'){
      WOA13PFile <- '~/ROMS/Data/WOA13/WOA13_PO4.Rdata'
      load(WOA13PFile)
      data_s <- woa13po4
    } else if (varname == 'wSODA'){
      data_s <- wSODA_data
    } else if (varname == 'wROMS'){
      data_s <- wROMS_data
    } else if (varname == 'ssh'){
      data_s <- SSH_data
    } else if (varname == 'Chl_C'){
      data_s <- Chl_C
    } else if (varname == 'CESM.NPP'){
      data_s <- NPP.CESM
    } else{
      stop('Variable name incorrect!')
    }
    Lon_s  <- data_s$lon
    Lat_s  <- data_s$lat
    Lon_s[Lon_s>0] <- Lon_s[Lon_s>0]-360
    # Obtain the data:
    Z      <- data_s$data

    # Determine the time dimension of the data
    time            <- data_s$time
    lon[lon>0]      <- lon[lon>0]-360
    gridlist        <- list(x=lon,y=lat)
    obs_time        <- data.frame(matrix(NA, nr=1, nc=length(time) ) )
    names(obs_time) <- c( paste('M',1:length(time),sep=''))

    if (ROMS){
      obs_time[1,]    <- time/30 
    } else{
      obs_time[1,]    <- time 
    }

    if(length(dim(Z)) == 3) {
       dat  <- data.frame(matrix(NA, nr=1, nc=length(time)+1 ) )
      names(dat) <- c('Depth', paste('M',1:length(time),sep=''))
       dat[,1]   <- 0

      for (i in 1:length(time) ) {
          d    <- list(x=Lon_s, y=Lat_s, z=Z[,,i])
        dat[1,i+1] <- interp.surface.grid(d,gridlist)$z
      }

     } else if (length(dim(Z)) == 4){
      depth <- data_s$depth

      if (ROMS){
        N = dim(depth)[3]
       dat= numeric(N)
        for (k in 1:N){
            d      <- list(x=Lon_s, y=Lat_s, z=depth[,,k] )    
            dat[k] <- interp.surface.grid(d,gridlist)$z 
        }
       depth <- abs(dat)
      } 

       dat   <- data.frame(matrix(NA, nr=length(depth), nc=length(time)+1 ) )

      names(dat) <- c('Depth', paste('M',1:length(time),sep=''))
       dat[,1]   <- -depth
      for (i in 1:length(depth)) {
          for (j in 1:length(time) ){
              d <- list(x=Lon_s, y=Lat_s, z=Z[,,i,j])    
            dat[i,j+1] <- interp.surface.grid(d,gridlist)$z
          }
      }

     } else{
     stop('Data dimension incorrect!')
     }
     return(list(time=obs_time,data=dat))
}


