setwd('~/Working/FlexEFT1D/DRAM')  
source('~/Working/FlexEFT1D/Rscripts/plot_1D.R')
VALs     <- c('muN1','The1','QN_1')
Models   <- c('NPZDFix_sRun','EFTsimple_sRun')
Modnames <- c('MONOD','PAHLOW')

#Turn off pdf device:

while (length(dev.list()) > 0) dev.off()
PDFfile <- paste0('muthetaQN_NPZDFix_EFTsimple_1D.pdf')
pdf(PDFfile, width=9,height=6,paper='a4')
op <- par(font.lab = 1, las = 1,
             family ="serif",
             mar    = c(2,3,2.5,2),
             mgp    = c(2.3,1,0),
             oma    = c(4,4,0,0),
             mfrow  = c(3,4)   ) 
for (v in VALs){
  if (v == 'muN1'){
     ZLIM <- c(1E-4, .8)
  }else if (v == 'The1'){
     ZLIM <- c(0, 0.6)
  } else if (v =='QN_1'){
     ZLIM <- c(0, 0.2)
  }

  Stns <- c('HOT','S1')
  for (n in 1:length(Stns)){
    Stn <- Stns[n]
    for (m in 1:length(Models)){
      title <- Modnames[m]
      model <- Models[m]
      plot_1D(v,model,Stn,ZLIM=ZLIM,finalyr=T,BOTH=F, Dmax = -200)
      if (v == 'muN1') {
        mtext(title, side = 3, adj=0.5, line=1)
        if (Stn == 'HOT') Stn = 'ALOHA'
        mtext(Stn, side = 3,   adj=0, cex = .8)
        if (Stn == 'ALOHA') Stn = 'HOT'
      }
    }
  }
}
par(las=0)
mtext('Depth (m)', side = 2, outer = T, line = 1)
dev.off()


