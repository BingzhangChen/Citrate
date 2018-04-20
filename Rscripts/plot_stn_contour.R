source('~/Working/FlexEFT1D/Rscripts/plot_1D.R')
plot_stn <- function(Stns, VARs, Model='NPZDcont', finalyr = T, Dmax = -150){
   NVar = length(VARs)
   if (length(Stns) > 1){
      Both = T
   }else{
      Both = F
   }
   for (Stn in Stns){
     pdffile <- paste0(Stn,'_1D.pdf')
     pdf(pdffile, width = 8, height = 8, paper = 'a4')
     op <- par(font.lab = 1,
                 family ="serif", cex.axis=1.2, cex.lab=1.2,
                 mar    = c(2,2,1.5,3.5),
                 mgp    = c(2.3,1,0),
                 mfrow  = c(ceiling(NVar/2),2),
                 oma    = c(4,4,1,0)) 
     for (i in 1:NVar){
         VAR = VARs[i]
         plot_1D(VAR,Model,Stn,finalyr = finalyr,BOTH=Both, Dmax = Dmax)
         mtext(paste0(letters[i],')'),adj=0, outer=F)
     }
     mtext('Depth (m)', side = 2, outer=TRUE, line=2)
     if (Stn == 'HOT') Stn = 'ALOHA'
     mtext(paste('Modelled seasonal patterns at',Stn),side=1,outer=T, line=2,adj=0)
     dev.off()
   }
}

