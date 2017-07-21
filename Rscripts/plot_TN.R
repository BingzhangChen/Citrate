source('~/Working/FlexEFT1D/Rscripts/Hz.R')
Hz      <- Hz(hmax=250, thetaS=2, nlev=30)
plot_TN <- function(Stn, model_ID,BOTH = T){
   if (BOTH){
    filedir  <- paste0('~/working/FlexEFT1D/DRAM_0.9/',model_ID,'/BOTH/')
   } else{
    filedir  <- paste0('~/working/FlexEFT1D/DRAM_0.9/',model_ID,'/',Stn,'/')
   }

   setwd(filedir)
   NO3   <- getData(filedir,Stn,'NO3')
   ZOO   <- getData(filedir,Stn,'ZOO')$data
   DET   <- getData(filedir,Stn,'DET')$data
   days  <- NO3$days
   depth <- NO3$depth
   NO3   <- NO3$data
   PHY   <- getData(filedir,Stn,'PHY_T')$data
   #Calculate total N for each date:
   TNO3 <- numeric(length(days)) 
   TZOO <- numeric(length(days)) 
   TDET <- numeric(length(days)) 
   TPHY <- numeric(length(days)) 
   
   for (i in 1:length(days)){
   
       TNO3[i] <- sum(NO3[i,]*Hz)   #Unit: mmol/m2
       TDET[i] <- sum(DET[i,]*Hz)   #Unit: mmol/m2
       TZOO[i] <- sum(ZOO[i,]*Hz)   #Unit: mmol/m2
       TPHY[i] <- sum(PHY[i,]*Hz)   #Unit: mmol/m2
   }
 
   TN      <- TNO3+TZOO+TDET+TPHY
   pdffile <- paste(Stn,model_ID,'_TN.pdf',sep='')
   pdf(pdffile, width=5,height=6,paper='a4')
   plot(days,TN,type='l',lwd=2,ylim=c(0,max(TN)))
   lines(days,TNO3,col='red')
   lines(days,TDET,col='blue')
   lines(days,TZOO,col='green')
   lines(days,TPHY,col='cyan',lty=2)
   legend('topright',legend=c('TN','NO3','DET','ZOO','PHY'),
          col=c(1,'red','blue','green','cyan'),
          lty=c(1,1,1,1,3))
   dev.off()
}

