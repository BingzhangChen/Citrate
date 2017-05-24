Model <- 'NPZDN2'
Stn   <- 'HOT'
DIR   <- paste0('~/Working/FlexEFT1D/DRAM_0.9/',Model,'/',Stn,'/')
setwd(DIR)
source('~/Working/FlexEFT1D/Rscripts/plot_1D.R')

#Plot time evolution of TN:
TN   <- getData(DIR,Stn,'TN')
plot(TN$days, TN$data, type = 'l')

#Plot PON:

plot_1D('PON',Model,Stn,title='',ZLIM=NULL, finalyr = F)
