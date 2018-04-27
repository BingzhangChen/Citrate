#Plot station map:
source('~/Working/FlexEFT1D/Rscripts/Stnmap_Kv_S1K2.R')
Model <- 'NPZclosure'
DIR   <- paste0('~/Working/FlexEFT1D/DRAM/',Model,'/S1/')
setwd(DIR)
source('~/Working/FlexEFT1D/Rscripts/plot_1D.R')
source('~/Working/FlexEFT1D/Rscripts/loglike_params_NPclos.R')  #Plot time-evolution of Log-Likelihoods and parameters 

#Check total concentration and variance:
source('~/Working/FlexEFT1D/Rscripts/plot_TN_TVAR.R')
plot_TN_TVAR('NPZclosure',Stns='HOT')
#Plot an example of five years to show seasonal cycle:
if (model == 'NPclosure'){
   VARs    <- c('NO3','PHY','VNO3','VPHY','COVNP')
}else if (model == 'NPZclosure'){
   VARs    <- c('NO3','PHY','ZOO','VNO3','VPHY','VZOO','COVNP','COVNZ', 'COVPZ')
}
NVar    <- length(VARs)
pdffile <- paste0('fullyear_example_1D.pdf')

pdf(pdffile, width = 9, height = 6, paper = 'a4')
op <- par(font.lab = 1,
            family ="serif", cex.axis=1.2, cex.lab=1.2,
            mar    = c(2,2,1.5,3.5),
            mgp    = c(2.3,1,0),
            mfcol  = c(ceiling(sqrt(NVar)),ceiling(sqrt(NVar))),
            oma    = c(4,4,1,0)) 

j <- 0
for (Stn in c('HOT','S1')){
   for (i in 1:NVar){
       VAR = VARs[i]
       plot_1D(VAR,Model,Stn,finalyr = F, Dmax = -500)
       j   = j + 1
       if (i == 1){
         mtext(paste0(letters[j],')',Stn),adj=0, outer=F)
       } else{
         mtext(paste0(letters[j],')'),adj=0, outer=F)
       }
   }
}
mtext('Depth (m)', side = 2, outer=TRUE, line=2)
mtext('Days     ', side = 1, outer=TRUE, line=1)
mtext('An example of modelled 5 year patterns at HOT and S1',side=1,outer=T, line=2,adj=0)
dev.off()

#One plot for one station (only final year):
if (model == 'NPclosure'){
   VARs <- c('NO3','PHY_T','CHL_T','VNO3','VPHY', 'COVNP')
}else if{
   VARs <- c('NO3','PHY_T','CHL_T','VNO3','VPHY', 'VZOO', 'COVNP', 'COV')
}
Stns <- c('HOT')

source('~/Working/FlexEFT1D/Rscripts/plot_stn_contour.R')
plot_stn(Stns, VARs, Model='NPclosure', finalyr = T, Dmax = -180)

#Plot comparisons of vertical profiles between data and model based on best parameter
COLS     <- 2:3
Models   <- c('NPZclosure')
Modnames <- c('NPZ closure')
Stns     <- c('HOT')
Nstn     <- length(Stns)
DIR      <- paste0('~/Working/FlexEFT1D/DRAM/',Models,'/',Stns,'/')
setwd(DIR)

source('~/Working/FlexEFT1D/Rscripts/plot_vertical_NChlNPP.R')
for (Stn in Stns){
    plot_v_n(Stn, Models, VARS=c('DIN','CHL','NPP'), BOTH=F)
}


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

