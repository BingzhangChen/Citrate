#Plot the model output of surface Chl:
setwd('~/Working/FlexEFT1D/DRAM_0.9/')  
source('~/Working/FlexEFT1D/Rscripts/getData.R')
source('../Rscripts/get_obs_MLD.R')
source('../Rscripts/Hz.R')
COLS <- c('green','red')
LTYS <- c(rep(1,3),rep(3,3))
Stns <- c('S1','HOT')
Nstn <- length(Stns)
Models <- c('NPZDFix_sRun','EFTsimple_sRun')

pdf('Fig5_ML_mod_obs.pdf', width=2.5*Nstn,height=6,paper='a4')
   op <- par(font.lab = 1,
                family ="serif",
                mar    = c(4,4,1,3),
                oma    = c(2,2,0,0),
                mgp    = c(2.3,1,0),
                mfrow  = c(3,Nstn)) 
j=0
for (Var in c('TIN','CHL','NPP')){
  for (Stn in Stns){

    j       = j+1
    MLD     = get_MLD(Stn)
    Dat     = get_obs_MLD(Stn,Var)
    ymax    = max(Dat$dat)*2
    XLab    = ''
    YLab    = ''
    par(mar=c(2,2,2,.2))
    #Write out Ylab:
    if (Var == 'TIN'){
       if (Stn == Stns[1]) {
           YLab = expression(paste('TIN (µmol '*L^-1*')'))
           par(mar=c(2,4,2,.2))
       }
       if (Stn == 'HOT') ymax=0.5
    }else if(Var == 'CHL'){
       if (Stn == Stns[1]){
          YLab=expression(paste("Chl "*italic(a)*' (µg '*L^-1*')'))
          par(mar=c(2,4,2,.2))
       }
       if (Stn == 'K2') {
           ymax=20
       }else if(Stn == 'HOT'){
           ymax=.6
       }
    } else if(Var == 'NPP'){
       if (Stn==Stns[1]){
         YLab=expression(paste('NPP  (mg C '*m^-3*' '*d^-1*')'))
         par(mar=c(2,4,2,.2))
       }
       if (Stn=='K2') ymax=140
    }

    #Plot data points:
    plot(Dat$DOY,Dat$dat,
      xlab=XLab,
      ylab=YLab,
      xlim=c(0,360),ylim=c(0,ymax),cex=.5,cex.lab=1.2,pch=16,xaxt='n')
    axis(1, at=seq(0,360,by=60))
    mtext(letters[j],adj=0,cex=1)
    if (j <= Nstn) {
        if (Stn == 'HOT'){
          mtext('ALOHA',adj=.5)
        }else{
          mtext(Stn,adj=.5)
        }
    }

    ii <- 0
    for (model in Models){
       ii      <- ii+1
       DIR     <- paste('~/Working/FlexEFT1D/DRAM_0.9/',model,'/',Stn,'/',sep='')

       #Get modeled Chl data
       if (Var == 'CHL'){
        Chl  = getData(DIR,Stn,'CHL_T')
       }else if (Var == 'TIN'){
        Chl  = getData(DIR,Stn,'NO3')
       }else if(Var == 'NPP'){
        Chl  = getData(DIR,Stn,'NPP_T')
       }
        days = Chl$days
       depth = Chl$depth
        Chl  = Chl$data
        d_per_y = 360
        #Get the data of the final year
        w       = (nrow(Chl)-d_per_y+1):nrow(Chl)
        Chl     = Chl[w,]

        s_Chl   = numeric(length(MLD))

     #   if (Var == 'NPP'){
     #   #Calculate integrated values within the euphotic zone
     #      HZ <- Hz()
     #      for (k in 1:length(MLD)){
     #          s_Chl[k] <- sum(HZ*Chl[k,])
     #      }
     #   }else{
           
        #calculate the average data within the MLD
           for (k in 1:length(MLD)){
              w    = which(depth >= -abs(MLD[k]))
              xff  = mean(as.numeric(Chl[k,w]))
           s_Chl[k]= xff
           }
     #   }
        lines(1:length(s_Chl),s_Chl, col=COLS[ii])

        if (Var == 'TIN' && Stn == Stns[1]){
           legend('topright',c('MONOD','PAHLOW'),
                  col=COLS,lty=1  )
        }
    }
  }
}
mtext('Date of the year',side=1,line=1,outer=TRUE)
dev.off()
