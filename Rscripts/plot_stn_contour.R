source('~/Working/FlexEFT1D/Rscripts/plot_1D.R')
plot_stn <- function(Stns, VARs, Model='NPZDcont', BOTH = F, finalyr = T, Dmax = -150){
   NVar = length(VARs)
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

plot_full <- function(VARs, Stns = c('HOT', 'S1'), 
                      Model   = 'NPZclosure', Dmax = -500,
                      pdffile = 'fullyear_example_1D.pdf'){
   NVar   <- length(VARs)
   source('~/Working/FlexEFT1D/Rscripts/plot_1D.R')

   pdf(pdffile, width = 9, height = 6, paper = 'a4')
   op <- par(font.lab = 1,
               family ="serif", cex.axis=1.2, cex.lab=1.2,
               mar    = c(2,2,1.5,3.5),
               mgp    = c(2.3,1,0),
               mfcol  = c(ceiling(sqrt(NVar)),ceiling(sqrt(NVar))),
               oma    = c(4,4,1,0)) 
   
   j <- 0
   for (Stn in Stns){
      for (i in 1:NVar){
          VAR = VARs[i]
          plot_1D(VAR,Model,Stn,finalyr = F, Dmax = Dmax)
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
   mtext(paste0('An example of modelled 5 year patterns at ', 
                paste0(Stns, collapse=',')),side=1,outer=T, line=2,adj=0)
   dev.off()
}

