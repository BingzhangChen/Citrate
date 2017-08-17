source('~/Working/FlexEFT1D/Rscripts/Hz.R')
source('~/Working/FlexEFT1D/Rscripts/getData.R')
Hz        <- Hz(hmax=500, thetaS=2, nlev=30)
model_ID  <- 'NPZDcont'
BOTH      <- T
Stns      <- c('S1', 'K2')
for (Stn in Stns){
  if (BOTH){
   filedir  <- paste0('~/working/FlexEFT1D/DRAM/',model_ID,'/BOTH_TD/')
  } else{
   filedir  <- paste0('~/working/FlexEFT1D/DRAM/',model_ID,'/',Stn,'/')
  }
  
  setwd(filedir)
  NO3   <- getData(filedir,Stn,'NO3')
  MIC   <- as.matrix(getData(filedir,Stn,'MIC')$data)
  MES   <- as.matrix(getData(filedir,Stn,'MES')$data)
  DET   <- as.matrix(getData(filedir,Stn,'DET')$data)
  days  <- NO3$days
  depth <- NO3$depth
  NO3   <- as.matrix(NO3$data)
  PHY   <- as.matrix(getData(filedir,Stn,'PHY_T')$data)
  #Calculate total N for each date:
  TNO3 <- NO3 %*% Hz
  TMIC <- MIC %*% Hz
  TMES <- MES %*% Hz
  TDET <- DET %*% Hz
  TPHY <- PHY %*% Hz
  TN   <- TNO3+TMIC+TDET+TPHY+TMES
  
  pdffile <- paste(Stn,model_ID,'_TN_TFe.pdf',sep='')
  
  pdf(pdffile, width=5,height=3*length(Stns),paper='a4')
    op <- par(font.lab = 1,
               family ="serif",
               mar    = c(4,4,1,3),
               oma    = c(2,2,0,0),
               mgp    = c(2.3,1,0),
               lwd    = 1.5,
               mfrow  = c(length(Stns),1)) 
    
    plot(days,TN,type='l',lwd=2,ylim=c(0,max(TN)))
    lines(days,TNO3,col='red')
    lines(days,TDET,col='blue')
    lines(days,TMIC,col='green')
    lines(days,TMES,col='tan')
    lines(days,TPHY,col='cyan',lty=2)
    legend('topright',legend=c('TN','NO3','DET','MIC','MES','PHY'),
           col=c(1,'red','blue','green','tan','cyan'),
           lty=c(1,1,1,1,1,2))
    
    #Calculate total Iron:
    fer   <- as.matrix(getData(filedir,Stn,'Fer'  )$data)
    DETfe <- as.matrix(getData(filedir,Stn,'DETFe')$data)
    Fe_N  <- 0.0265
    Tfer  <- fer   %*% Hz 
    TDETfe<- DETfe %*% Hz
    TFE   <- Tfer + TDETfe + (TMIC + TPHY)*Fe_N
    plot(days,TFE,type='l',lwd=2,ylim=c(0,max(TFE)))
    lines(days,Tfer,  col='red')
    lines(days,TDETfe,col='blue')
    lines(days,TMIC*Fe_N,col='green')
    lines(days,TMES*Fe_N,col='tan')
    lines(days,TPHY*Fe_N,col='cyan',lty=2)
    legend('topright',legend=c('TFE','Tfer','TDETfe','MIC','MES','PHY'),
           col=c(1,'red','blue','green','tan','cyan'),
           lty=c(1,1,1,1,1,2))
  
  dev.off()
}

