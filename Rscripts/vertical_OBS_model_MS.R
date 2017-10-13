setwd('~/Working/FlexEFT1D/DRAM')  
source('~/Working/FlexEFT1D/Rscripts/getData.R')
source('~/Working/FlexEFT1D/Rscripts/Varname.R')
COLS     <- c(3,2)
Models   <- c('NPZDFix_sRun','EFTsimple_sRun')
Modnames <- c('MONOD', 'PAHLOW')
Stns     <- c('S1','HOT')
Nstn     <- length(Stns)

for (Stn in Stns){
 fname <- paste0('Vertical_mod_obs_NChlNPP_',Stn,'.pdf')
 pdf(fname, width=2*3,height=8,paper='a4')
 
 op <- par(font.lab = 1,
              family ="serif",
              mar    = c(2,2,1.5,0.5),
              mgp    = c(2.3,1,0),
              oma    = c(4,4,0,0),
              mfcol  = c(4,3)   ) 
 
 j <- 0
 for (Var in c('DIN','CHL','NPP')){
    VARNAME = Varname(Var)
    #Read observational data
    dat <- paste('~/Working/FlexEFT1D/Forcing/',Stn,'_',Var,'.dat',sep='')
    dat <- read.table(dat,header=T)
    #Average into 4 seasons
    DOYs    <- seq(0,360,length.out=5)
    seasons <- c('Winter','Spring','Summer','Fall')
    for (i in 1:4){
      if (i == 4) {
         par(mar=c(4,2,1.5,0.5))
      }else{
         par(mar=c(2,2,1.5,0.5))
      }
      j    <- j+1
      dat_ <- dat[dat$DOY > DOYs[i] & dat$DOY <= DOYs[i+1],]
      xmax <- max(dat_[,3])
     
     # if (i==1) mtext(Stn,adj=.7)
      ii = 0
      for (model in Models){
         ii    <- ii+1
          DIR  <- paste0('~/Working/FlexEFT1D/DRAM/',model,'/',Stn,'/')
         #Get modeled data
         if (Var == 'CHL'){
          Chl  <- getData(DIR,Stn, 'CHL_T')
         }else if (Var == 'DIN'){
          Chl  <- getData(DIR,Stn, 'NO3')
         }else if(Var == 'NPP'){
          Chl  <- getData(DIR, Stn, 'NPP_T')
         }else if(Var == 'C_Chl'){
          Chl  <- getData(DIR, Stn,'The1')
         }else{
          Chl  <- getData(DIR,Stn,Var)
         }
          days    = Chl$days
         depth    = Chl$depth
          Chl     = Chl$data
          d_per_y = 360
          #Get the data of the final year
          w       = (nrow(Chl)-d_per_y+1):nrow(Chl)
          Chl     = Chl[w,]

          #Read model results:
          mdoy  = 1:nrow(Chl)
          k     = which(mdoy > DOYs[i] & mdoy <= DOYs[i+1])
          Chl_  = Chl[k,]
          #Calculate quantiles:
          cff   = sapply(1:ncol(Chl_),function(k)quantile(Chl_[,k], probs=c(0.025,0.5,0.975)))
          xmax  = max(xmax, t(cff)[,2])*1.2 #Determine maximal x value
          #Plot model and data:
          if (model == Models[1]){
             plot(dat_[,3],-dat_$Depth,xlim=c(0,xmax), ylim=c(-250,0),
                    xlab=VARNAME,ylab='',pch=16,cex=.5,cex.lab=1.2)
             mtext(paste(letters[j],')',seasons[i]),adj=0)
          }

          matlines(t(cff), depth, lty=c(2,1,2), lwd=c(.5,1.5,.5), col=COLS[ii])
      }
      if (Var == 'DIN' && i == 1){
         legend('topright',Modnames,col=COLS,lty=1)
      }
    }
  }
  mtext('Depth (m)',side = 2, outer=TRUE, line=1)
dev.off()
}
