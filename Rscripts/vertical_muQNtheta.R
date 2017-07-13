setwd('~/Working/FlexEFT1D/DRAM_0.9')  
source('~/Working/FlexEFT1D/Rscripts/plot_1D.R')
VALs     <- c('muN1','The1','QN_1')
Models   <- c('NPZDFix_sRun','EFTsimple_sRun')
#Modnames <- c('NPZDFix','Geider','Pahlow','NPZDFixFe','GeiderFe','PahlowFe')
Modnames <- c('MONOD','PAHLOW')

#Turn off pdf device:

while (length(dev.list()) > 0) dev.off()
for (v in VALs){
  if (v == 'muN1'){
     ZLIM <- c(1E-4, .8)
  }else if (v == 'The1'){
     ZLIM <- c(0, 0.4)
  } else if (v =='QN_1'){
     ZLIM <- c(0, 0.3)
  }
  PDFfile <- paste0(v,'NPZDFix_EFTsimple_1D.pdf')
  pdf(PDFfile, width=8,height=5,paper='a4')
  op <- par(font.lab = 1, las = 1,
               family ="serif",
               mar    = c(2,3,2.5,2),
               mgp    = c(2.3,1,0),
               oma    = c(4,4,0,0),
               mfcol  = c(2,2)   ) 

  Stns <- c('HOT','S1')
  for (n in 1:length(Stns)){
    Stn <- Stns[n]
    for (m in 1:length(Models)){
      title <- Modnames[m]
      model <- Models[m]
      plot_1D(v,model,Stn,ZLIM=ZLIM,finalyr=T,BOTH=F)
      mtext(title, side = 3, adj=0.5, line=1)
      if (Stn == 'HOT') Stn = 'ALOHA'
      mtext(Stn, side = 3, adj=0)
      if (Stn == 'ALOHA') Stn = 'HOT'
    }
  }
  dev.off()
}


