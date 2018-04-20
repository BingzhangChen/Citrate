source('~/Working/FlexEFT1D/Rscripts/Hz.R')
getData <- function(fileDIR,Stn,VAR){
   #Get all the data:
   dir    <- paste(fileDIR,Stn,sep='')
   dir    <- as.character(dir)
   file   <- paste(dir,'.out',sep='')
   data   <- read.table(file,header=F)
   data[,1] <- as.character(data[,1])
   depth    <- as.numeric(data[1,4:ncol(data)])

   if (VAR == 'TN'){
     Hz     <- Hz()
      #Read PHY data:
     PHY    <- data[data[,1] == 'PHY1', ]
     days   <- as.numeric(as.character(PHY[, 3]))
     PHY    <- PHY[,  4:ncol(PHY)]
     #Read NO3 data:
     NO3    <- data[data[,1] == 'NO3', ]
     NO3    <- NO3[,  4:ncol(NO3)]
     #Read ZOO data:
     ZOO    <- data[data[,1] == 'ZOO', ]
     ZOO    <- ZOO[,  4:ncol(ZOO)]
     #Read DET data:
     DET    <- data[data[,1] == 'DET', ]
     DET    <- DET[,  4:ncol(DET)]
     TN     <- PHY+NO3+ZOO+DET
     for (i in 1:nrow(TN)) TN[i,] <- TN[i,]*Hz
     data   <- apply(TN, 1, sum)

   }else if (VAR == 'R_PMU' || VAR == 'R_VAR'){
     #Read PHY data:
     PHY    <- data[data[,1] == 'PHY1', ]
     PHY    <- PHY[,  4:ncol(PHY)]
     #Read PMU or VAR data:
     NEWVAR <- substr(VAR,3,5)
     dat    <- data[data[,1] == NEWVAR, ] 
     days   <- as.numeric(as.character(dat[, 3]))
     dat    <- dat[, 4:ncol(dat)]
     if (VAR == 'R_PMU'){
     #Get real PMU
        data <- dat/PHY
        data <- (exp(data-log(10)) * 6/pi)**0.3333333333333
     }else if (VAR == 'R_VAR'){
        #Get PMU
        PMU  <- data[data[,1] == 'PMU', ]
        PMU  <- PMU[, 4:ncol(PMU)]
        PMU  <- PMU/PHY
        data <- dat/PHY - PMU^2
     }
   } else if (VAR == 'TD_VAR'){
     #TD_VAR is the contribution of "trait diffusion" to changes of size variance
    # (VAR*(VAR*(d2muNetdl2-d2gdl2bar+VTR*d4mudl4)-5d0*VTR*d2muNetdl2)+2d0*VTR*muNet)

    #Read d4mudl4
     d4mu   <- data[data[,1] == 'd4mu', ]
     days   <- as.numeric(as.character(d4mu[, 3]))
     d4mu   <- d4mu[, 4:ncol(d4mu)]
    #Read d2mudl2
     d2mu   <- data[data[,1] == 'd2mu', ]
     d2mu   <- d2mu[, 4:ncol(d2mu)]
    #Read muNet
     muNet  <- data[data[,1] == 'muN1', ]
     muNet  <- muNet[, 4:ncol(muNet)]
    #Read real VAR data:
     PHY    <- data[data[,1] == 'PHY1', ]
     PHY    <- PHY[,  4:ncol(PHY)]
     VAR    <- data[data[,1] == 'VAR', ] 
     VAR    <- VAR[, 4:ncol(VAR)]
     #Get real VAR
     VAR    <- VAR/PHY**2

    #Calculate Trait diffusion values:
     #Read VTR:
     enspar <- paste0(fileDIR,'enspar')
     if (file.exists(enspar)){
        pars   <- read.table(enspar, header=T)
        VTR    <- pars[nrow(pars),'VTR']
     }else {
        VTR    <- 0.08  # single run
     }
    #     Total  <- VAR*(VAR*(d2mu+VTR*d4mu)-5*VTR*d2mu)+2*VTR*muNet
     TD     <- VAR*(VAR*(VTR*d4mu)-5*VTR*d2mu)+2*VTR*muNet 
     data   <- TD

   } else{
     stopifnot(VAR %in% data[,1])
     data <- data[data[,1] == VAR, ]
     days <- as.numeric(as.character(data[, 3]))
     data <- data[, 4:ncol(data)]
   }
   return(list(days=days,depth=depth,data=data))
}

