plot_TN_TVAR <- function(model_ID = 'NPclosure', Stns=c('HOT','S1')){
   source('~/Working/FlexEFT1D/Rscripts/Hz.R')
   source('~/Working/FlexEFT1D/Rscripts/getData.R')
   Hz        <- Hz(hmax=500, thetaS=2, nlev=40)
   BOTH      <- FALSE
   for (Stn in Stns){
     if (BOTH){
      filedir  <- paste0('~/working/FlexEFT1D/DRAM/',model_ID,'/BOTH_TD/')
     } else{
      filedir  <- paste0('~/working/FlexEFT1D/DRAM/',model_ID,'/',Stn,'/')
     }
     setwd(filedir)
     NO3   <- getData(filedir,Stn,'NO3')
     days  <- NO3$days
     depth <- NO3$depth
     NO3   <- as.matrix(NO3$data)
     PHY   <- as.matrix(getData(filedir,Stn,'PHY')$data)
     #Calculate total N for each date:
     TNO3 <- NO3 %*% Hz
     TPHY <- PHY %*% Hz
     
     #Obtain variances and covariances:
     VPHY    <- as.matrix(getData(filedir,Stn,'VPHY')$data)
     VNO3    <- as.matrix(getData(filedir,Stn,'VNO3')$data)
     COVNP   <- as.matrix(getData(filedir,Stn,'COVNP')$data)
     TVNO3   <- VNO3  %*% Hz
     TVPHY   <- VPHY  %*% Hz
     TCOVNP  <- COVNP %*% Hz

     if (model_ID == 'NPZclosure'){
         ZOO   <- as.matrix(getData(filedir,Stn,'ZOO')$data)
        VZOO   <- as.matrix(getData(filedir,Stn,'VZOO')$data)
        COVNZ  <- as.matrix(getData(filedir,Stn,'COVNZ')$data)
        COVPZ  <- as.matrix(getData(filedir,Stn,'COVPZ')$data)
        TZOO   <-  ZOO  %*% Hz
        TVZOO  <- VZOO  %*% Hz
        TCOVNZ <- COVNZ %*% Hz
        TCOVPZ <- COVPZ %*% Hz
          TN   <- TNO3  + TPHY  + TZOO
        TVAR   <- TVPHY + TVNO3 + TVZOO + 2.*(TCOVNP + TCOVPZ + TCOVNZ) 
     }else if (model_ID == 'NPclosure'){
          TN   <- TNO3  + TPHY
        TVAR   <- TVPHY + TVNO3 + 2.*TCOVNP 
     }else{
        stop("Model name incorrect!")
     }
     pdffile <- paste(Stn,model_ID,'_TN_TVAR.pdf',sep='')
     pdf(pdffile, width=5,height=3*length(Stns),paper='a4')
       op <- par(font.lab = 1,
                  family ="serif",
                  mar    = c(4,4,1,3),
                  oma    = c(2,2,0,0),
                  mgp    = c(2.3,1,0),
                  lwd    = 1.5,
                  mfrow  = c(length(Stns),1)) 
       
       plot(days,TN,type='l',lwd=2,ylim=c(0,max(TN)))
       lines(days,TNO3,col=2,lty=2)
       lines(days,TPHY,col=3,lty=2)

       if (model_ID == 'NPZclosure'){
          lines(days,TZOO,col=4,lty=2)
          legend('topright',legend=c('TN','NO3','PHY','ZOO'),col=1:4,lty=c(1,2,2,2))
       }else{
          legend('topright',legend=c('TN','NO3','PHY'),col=1:3,lty=c(1,2,2))
       }
       
       plot(days,TVAR,type='l',lwd=2,ylim=c(0,max(TVAR)))
       lines(days,TVNO3,  col=2,lty=2)
       lines(days,TVPHY,  col=3,lty=2)
       lines(days,TCOVNP, col=4,lty=2)
       if (model_ID == 'NPZclosure'){
          lines(days,TVZOO, col=5,lty=2)
          lines(days,TCOVPZ,col=6,lty=2)
          lines(days,TCOVNZ,col=7,lty=2)
          legend('topright',legend=c('TVAR','VNO3','VPHY','COVNP','VZOO','COVPZ','COVNZ'),col=1:7,lty=c(1,rep(2,6)))
       }else{
          legend('topright',legend=c('TVAR','VNO3','VPHY','COVNP'),col=1:4,lty=c(1,2,2,2))
       }

     dev.off()
   }
}
