#Plot station map:
source('~/Working/FlexEFT1D/Rscripts/Stnmap_Kv_S1K2.R')
Model <- 'NPZclosure'
Stn   <- 'K2'
DIR   <- paste0('~/Working/FlexEFT1D/DRAM/',Model,'/',Stn,'/')
setwd(DIR)
source('~/Working/FlexEFT1D/Rscripts/loglike_params_NPclos.R')  #Plot time-evolution of Log-Likelihoods and parameters 

#Check total concentration and variance:
source('~/Working/FlexEFT1D/Rscripts/plot_TN_TVAR.R')
plot_TN_TVAR('NPZclosure',Stns='K2')
#Plot an example of five years to show seasonal cycle:
if (model == 'NPclosure'){
   VARs    <- c('NO3','PHY','VNO3','VPHY','COVNP')
}else if (model == 'NPZclosure'){
   VARs    <- c('NO3','PHY','ZOO','VNO3','VPHY','VZOO','COVNP','COVNZ', 'COVPZ')
}
source('~/Working/FlexEFT1D/Rscripts/plot_stn_contour.R')
plot_stn('K2', VARs, Model='NPZclosure', finalyr = F, Dmax = -200)

#One plot for one station (only final year):
if (model == 'NPclosure'){
   VARs <- c('NO3','PHY_T','CHL_T','VNO3','VPHY', 'COVNP')
}else if (model == 'NPZclosure'){
   VARs <- c('NO3','PHY','CHL','VNO3','VPHY', 'VZOO', 'COVNP', 'COVPZ','COVNZ')
}
Stns <- c('K2')
plot_stn(Stns, VARs = c('NO3','PHY','ZOO'), 
         Model = 'NPZclosure', finalyr = T, Dmax = -180)

#Plot comparisons of vertical profiles between data and model based on best parameter
COLS     <- 2:3
Models   <- c('NPZclosure')
Modnames <- c('NPZ closure')
Stns     <- c('S1')
Nstn     <- length(Stns)
DIR      <- paste0('~/Working/FlexEFT1D/DRAM/',Models,'/',Stns,'/')
setwd(DIR)

source('~/Working/FlexEFT1D/Rscripts/plot_vertical_NChlNPP.R')
for (Stn in Stns){
    plot_v_n(Stn, Models, VARS=c('DIN','CHL','NPP'), BOTH=F)
}

#Compare single run with different beta values (0.1 and 2.0):
Models   <- c('NPZclosure_sRun')
DIR      <- paste0('~/Working/FlexEFT1D/DRAM/',Models,'/',Stns,'/')
setwd(DIR)

Stn <- 'S1'
system('ln -s S1.out.beta0.001 S1.out')

plot_stn(Stns, VARs = c('NO3','PHY','ZOO'), 
         Model = 'NPZclosure', finalyr = T, Dmax = -180)
system('mv S1_1D.pdf S1_beta0.001.pdf')

plot_v_n(Stn, Models, VARS=c('DIN','CHL','NPP'), BOTH=F)
system('mv HOTNPZclosure_sRunV_NChlPP.pdf HOTNPZclosure_NChlPP_beta0.1.pdf')
plot_stn('HOT', VARs, Model='NPZclosure', finalyr = F, Dmax = -500)
system('mv HOT_1D.pdf HOT1D_beta0.1.pdf')
system('unlink HOT.out')

#Change to beta = 0.1
system('mv S1.out S1.out.beta0.1')
system('ln -s S1.out.beta0.1 S1.out')

plot_stn(Stns, VARs = c('NO3','PHY','ZOO'), 
         Model = 'NPZclosure_sRun', finalyr = T, Dmax = -180)
system('mv S1_1D.pdf S1_beta0.001.pdf')

plot_v_n(Stn, Models, VARS=c('DIN','CHL','NPP'), BOTH=F)
system('mv HOTNPZclosure_sRunV_NChlPP.pdf HOTNPZclosure_NChlPP_beta0.1.pdf')
plot_stn('HOT', VARs, Model='NPZclosure', finalyr = F, Dmax = -500)
system('mv HOT_1D.pdf HOT1D_beta0.1.pdf')
system('unlink HOT.out')

#Change to beta = 2.0
system('ln -s HOT.out.beta2.0 HOT.out')
plot_v_n(Stn, Models, VARS=c('DIN','CHL','NPP'), BOTH=F)
system('mv HOTNPZclosure_sRunV_NChlPP.pdf HOTNPZclosure_NChlPP_beta2.0.pdf')
plot_stn('HOT', VARs, Model='NPZclosure', finalyr = F, Dmax = -500)
system('mv HOT_1D.pdf HOT1D_beta2.0.pdf')
system('unlink HOT.out')
#Plot for HOT:
Model    <- 'NPZDcont_sRun'
Modnames <- 'NPZDcont'
Stn   <- 'HOT'
DIR   <- paste0('~/Working/FlexEFT1D/DRAM/',Model,'/',Stn,'/')
setwd(DIR)
fname='param.nml'
if(file.exists(fname)) file.remove(fname)
file.create(fname)
fw=file(fname,open='wt')
writeLines('&parameters',con=fw) #write header to file
txt=paste0('mu0hat=',bestpar$mu0hat,',')
writeLines(txt,con=fw) #write to file
txt=paste0('KN=',    bestpar$KN,',')
writeLines(txt,con=fw) #write to file
txt=paste0('KPHY=',  bestpar$KPHY,',')
writeLines(txt,con=fw) #write to file
txt=paste0('wDET=',  bestpar$wDET,',')
writeLines(txt,con=fw) #write to file
txt=paste0('aI0=',   bestpar$aI0_C,',')
writeLines(txt,con=fw) #write to file
txt=paste0('alphaI=',bestpar$alphaI,',')
writeLines(txt,con=fw) #write to file
txt=paste0('mz=',    bestpar$mz,',')
writeLines(txt,con=fw) #write to file
txt=paste0('KFe=',   bestpar$KFe,',')
writeLines(txt,con=fw) #write to file
txt=paste0('VTR=',   bestpar$VTR,',')
writeLines(txt,con=fw) #write to file
txt='bot_bound=0,'
writeLines(txt,con=fw) #write to file
writeLines('/',con=fw)
close(fw)
system('./run')
system('./citrate > Out')

plot_v_n(Stn, Model, VARS=c('DIN','CHL','NPP','PON'), BOTH=F)
plot_stn('HOT',VARS, Model='NPZDcont_sRun', finalyr = T, Dmax = -150)

