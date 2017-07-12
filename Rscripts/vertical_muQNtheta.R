setwd('~/Working/FlexEFT1D/DRAM_0.9/NPZDFix_sRun')  
source('~/Working/FlexEFT1D/Rscripts/plot_1D.R')
VALs     <- c('muN1','The1','QN_1')
Models   <- c('NPZDFix_sRun','EFTsimple_sRun')
#Modnames <- c('NPZDFix','Geider','Pahlow','NPZDFixFe','GeiderFe','PahlowFe')
Modnames <- c('MONOD','PAHLOW')

for (v in VALs){
  if (v == 'muN1'){
     ZLIM <- c(0, 1.4)
  }else if (v == 'The1'){
     ZLIM <- c(0, 0.4)
  } else if (v =='QN_1'){
     ZLIM <- c(0, 0.5)
  }
  PDFfile <- paste0(v,'_1D.pdf')
  pdf(PDFfile, width=8,height=5,paper='a4')
  op <- par(font.lab = 1, las = 1,
               family ="serif",
               mar    = c(2,2,1.5,2),
               mgp    = c(2.3,1,0),
               oma    = c(4,4,0,0),
               mfrow  = c(2,2)   ) 

  Stns <- c('HOT','S1')
  for (n in 1:length(Stns)){
    Stn <- Stns[n]
    for (m in 1:length(Models)){
      title <- Modnames[m]
      model <- Models[m]
      plot_1D(v,model,Stn,title=title,ZLIM=ZLIM,finalyr=T,BOTH=F)
    }
  }
  dev.off()
}


