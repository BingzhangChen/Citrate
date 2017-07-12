#Compare PHY, ZOO, and DET between the two models within ML
setwd('~/Working/FlexEFT1D/DRAM_0.8')  
source('../Rscripts/get_MLD.R')
source('~/Working/FlexEFT1D/Rscripts/interp.R')
source('~/Working/FlexEFT1D/Rscripts/MLDmean.R')
COLS <- c('green','red')
Stns <- c('S1')
Models <- c('NPZDFix_sRun','EFTsimple')
Nstn <- length(Stns)
pdf('PHY_ZOO.pdf', width=2.5*Nstn,height=7,paper='a4')
 op <- par(font.lab = 1,
             family ="serif",
             mar    = c(2,2,1.5,0),
             mgp    = c(2.3,1,0),
             mfrow  = c(3,Nstn),
             oma    = c(4,4,0,0)) 
j=0
for (Var in c('PHY1','ZOO','DET')){
  for (Stn in Stns){
    if (Stn == 'S1'){
       if (Var == 'PHY1'){
         Ymax = 1
       } else if (Var == 'ZOO'){
         Ymax = 1
       } else if (Var == 'DET'){
         Ymax = 1
       }
    } else if (Stn == 'K2'){
       if (Var == 'PHY1'){
         Ymax = 8
       } else if (Var == 'ZOO'){
         Ymax = 5
       } else if (Var == 'DET'){
         Ymax = 8
       }
    }

    XLab = ''
    YLab = ''
  #Obtain MLD data:
    MLD=get_MLD(Stn)
    j  =j+1
    ii =0
    N  =length(MLD)
    #Plot the general frame:
    plot(1:N,1:N,type='n',ylim=c(0,Ymax),
         xlab=XLab,ylab=YLab,xaxt='n')
    axis(1, at=seq(0,360,by=60))  #plot the axis every 2 months
    mtext(letters[j],adj=0,cex=1)
    if (j <= Nstn) mtext(Stn,adj=.5)

    for (model in Models){
       ii      =  ii+1
       filedir =  paste('~/working/FlexEFT1D/DRAM_0.8/',
                          model,'/',Stn,'/',sep='')
 
       s_dat = MLDmean(filedir,Stn,Var, MLD)
       lines(1:length(s_dat),s_dat, lwd=1,col=COLS[ii])

       #Plot legend:
       if (Var == 'PHY1' && Stn == Stns[1]) {
          legend('topright',c('MONOD','PAHLOW'),col=COLS,lty=1)
       }
    }
  }
}
mtext('Date of the year',side=1,adj=.5,outer=TRUE,line=.5)
vars=c('PHY','ZOO','DET')
dis =c(2.7,1.5,.3)
for (i in 1:length(vars)){
    txt  <- bquote(.(vars[i]) ~ ' (mmol '*m^-3*')')
    mtext(txt,side=2,adj=dis[i]/3,outer=TRUE)
}
par(op)
dev.off()
