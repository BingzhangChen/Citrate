source('~/Working/FlexEFT1D/Rscripts/Hz.R')
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
getData <- function(fileDIR,Stn,VAR){
   #Get all the data:
   dir    <- paste(fileDIR,Stn,sep='')
   dir    <- as.character(dir)
   file   <- paste(dir,'.out',sep='')
   data   <- read.table(file,header=F)
   NC     <- ncol(data)
   NR     <- nrow(data)
   data[2:NR,1] <- as.character(data[2:NR,1])
   depth        <- unlist(data[1,4:ncol(data)])

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

#Stn <- 'S1'
#
#simP <- function(Stn='S1'){
#  setwd(paste0('~/Working/FlexEFT1D/DRAM/NPZclosure_sRun/',Stn))
#  file     <- paste0(Stn,'.out')
#  data     <- read.table(file,header=F)
#  data[,1] <- as.character(data[,1])
#  
#  depth <- as.numeric(data[1,4:ncol(data)])
#  
#  NO3   = data[data[,1] == 'NO3', ]
#  days  = as.numeric(as.character(NO3[, 3]))
#  PHY   = data[data[,1] == 'PHY', ]
#  ZOO   = data[data[,1] == 'ZOO', ]
#  VPHY  = data[data[,1] == 'VPHY', ]
#  VNO3  = data[data[,1] == 'VNO3', ]
#  VZOO  = data[data[,1] == 'VZOO', ]
#  COVNP = data[data[,1] == 'COVNP', ]
#  COVNZ = data[data[,1] == 'COVNZ', ]
#  COVPZ = data[data[,1] == 'COVPZ', ]
#  #for (i in 4:ncol(data)) {
#  #    data[2:nrow(data), i] = as.numeric(as.character(data[2:nrow(data), i]))
#  #}
#  pdf('trial.pdf', width=5, height = 6)
#  par(mfrow=c(2,1))
#  plot(days, PHY[, ncol(PHY)], type = 'l', ylim = c(0,5))
#  points(days, NO3[, ncol(NO3)], type = 'l', col=2)
#  points(days, ZOO[, ncol(ZOO)], type = 'l', col=3)
#  legend('topright',c('PHY','NO3','ZOO'), lty=1, col=1:3)
#  
#  plot(days, VPHY[, ncol(VPHY)], type = 'l', ylim = c(-5, 10))
#  points(days, VNO3[, ncol(VNO3)], type = 'l', col=2)
#  points(days, VZOO[, ncol(VZOO)], type = 'l', col=3)
#  points(days, COVNP[, ncol(COVNP)], type = 'l', col=2, lty=2)
#  points(days, COVNZ[, ncol(COVNZ)], type = 'l', col=1, lty=2)
#  points(days, COVPZ[, ncol(COVPZ)], type = 'l', col=3, lty=2)
#  legend('topright',
#       c('VPHY','VNO3','VZOO','COVNP','COVNZ','COVPZ'), lty=c(1,1,1,2,2,2), col=c(1:3,2,1,3))
#  dev.off()
#}
#
#
##Check why PHY goes to extinction:
#P_N = data[data[,1] == 'P_N', ]
#Z_P = data[data[,1] == 'Z_P', ]
#N_P = data[data[,1] == 'N_P', ]
#
#plot(days, P_N[, ncol(P_N)], type = 'l', ylim = c(-50, 100))
#points(days, Z_P[, ncol(Z_P)], type = 'l', col=2)
#points(days, N_P[, ncol(N_P)], type = 'l', col=3)
#
#
##Check why COVNP ==> -inf
#COVNP = COVNP +    (PP_NP_PZ + PP_NP_PP + PP_NP_NN - PP_PP_NP - PP_NN_NP - PP_NZ_NP)*dtdays
#
#PP_NP = data[data[,1] == 'PP_NP', ]
#NZ_NP = data[data[,1] == 'NZ_NP', ]
#NP_PZ = data[data[,1] == 'NP_PZ', ]
#NP_NN = data[data[,1] == 'NP_NN', ]
#NN_NP = data[data[,1] == 'NN_NP', ]
#NP_PP = data[data[,1] == 'NP_PP', ]
#plot(days, PP_NP[, ncol(PP_NP)], type = 'l', lty = 2, ylim = c(-100,100))
#points(days, NN_NP[, ncol(NN_NP)], type = 'l', lty=2, col=2)
#points(days, NZ_NP[, ncol(NZ_NP)], type = 'l', lty=2, col=3)
#
#points(days, NP_PZ[, ncol(NP_PZ)], type = 'l', lty=1, col=1)
#points(days, NP_PP[, ncol(NP_PP)], type = 'l', lty=1, col=2)
#points(days, NP_NN[, ncol(NP_NN)], type = 'l', lty=1, col=3)
#
#plot(days, PP_NP[, ncol(PP_NP)], type = 'l')
#points(days, NP_PP[, ncol(NP_PP)], type = 'l', col=2)
