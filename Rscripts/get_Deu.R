get_Deu <- function(Stn,model,th=1E-2){
    #Par data: 
    DIR     <- paste('~/Working/FlexEFT1D/DRAM_0.8/',model,'/BOTH/',sep='')
    dat     <- getData(DIR,Stn,'PAR')
    days    <- dat$days
      depth <- dat$depth
       dat  <- dat$data

    d_per_y = 360
    #Get the data of the final year
    w       = (nrow(dat)-d_per_y+1):nrow(dat)
    dat     = dat[w,]
    Deu     = numeric(d_per_y)
    for (i in 1:d_per_y){
        cff = as.numeric(dat[i,])
        cff = cff/cff[length(cff)]
        Deu[i] = depth[which.max(cff >= th) ]
    }
    return(Deu)
}

