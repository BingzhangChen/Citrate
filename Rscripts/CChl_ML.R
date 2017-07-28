#Compare C:Chl ratios and growth/grazing rates between the two models within ML
setwd('~/Working/FlexEFT1D/DRAM')  
source('~/Working/FlexEFT1D/Rscripts/getData.R')
source('~/Working/FlexEFT1D/Rscripts/get_obs_MLD.R')
source('~/Working/FlexEFT1D/Rscripts/Interpolate_WOA.R')
#Estimate C:Chl from remote sensing:
Chl_C    <- readnc('Chl_C')  #A global Chl:C data on the surface
COLS     <- c('green','red')
Models   <- c('NPZDFix_sRun','EFTsimple_sRun')
#Modnames <- c('NPZDFix','Geider','Pahlow','NPZDFixFe','GeiderFe','PahlowFe')
Modnames <- c('MONOD','PAHLOW')
Stns     <- c('S1','HOT')
Nstn     <- length(Stns)

pdf('Fig8_CChl_mu_mod.pdf', width=2.5*Nstn,height=8,paper='a4')

op <- par(font.lab = 1, las = 1,
             family ="serif",
             mar    = c(4,4,1,3),
             oma    = c(2,2,0,0),
             mgp    = c(2.3,1,0),
             mfrow  = c(3,Nstn)) 
j=0
for (Var in c('mu','C_Chl','QN')){
  if (Var == 'mu'){
    Ymax = 1.2
  } else if (Var == 'QN'){
    Ymax = .5
  } else if (Var == 'C_Chl'){
    Ymax = .5
  }
  for (Stn in Stns){
    XLab    = ''
    YLab    = ''
    par(mar=c(2,2,2,.2))
    #Write out Ylab:
    if (Var == 'mu' && Stn == Stns[1]){
       YLab = expression(paste('Growth rate ('*d^-1*')'))
       par(mar=c(2,4,2,.2))
    }else if(Var == 'g'     && Stn==Stns[1]){
       YLab = expression(paste('Grazing rate ('*d^-1*')'))
       par(mar=c(2,4,2,.2))
    }else if(Var == 'C_Chl' && Stn==Stns[1]){
       YLab = expression(paste("Chl:C "*' (gChl '*molC^-1*')'))
       par(mar=c(2,4,2,.2))
    }else if(Var == 'QN' && Stn==Stns[1]){
       YLab = expression(paste("N:C "*' (molN '*molC^-1*')'))
       par(mar=c(2,4,2,.2))
    }
   
  #Obtain MLD data:
    MLD=get_MLD(Stn)
    j  =j+1
    ii =0
    N  =length(MLD)
    #Plot the general frame:
    plot(1:N,1:N,type='n',ylim=c(0,Ymax),cex.lab=1.2,
         xlab=XLab,ylab=YLab,xaxt='n')
    axis(1, at=seq(0,360,by=60))  #plot the axis every 2 months
    mtext(letters[j],adj=0,cex=1)
    if (j <= Nstn) {
        if (Stn == 'HOT'){
          mtext('ALOHA',adj=.5)
        }else{
          mtext(Stn,    adj=.5)
        }
    }

  #Extract model data:
    for (model in Models){
       ii  <-  ii+1
       DIR <- paste('~/Working/FlexEFT1D/DRAM/',model,'/',Stn,'/',sep='')
 
       #Get modeled data
       if (Var == 'mu'){
           dat  = getData(DIR,Stn, 'muN1')
       }else if (Var == 'QN'){
           dat  = getData(DIR,Stn,'QN_1')
       }else if (Var == 'C_Chl'){
           dat  = getData(DIR, Stn, 'The1')
       }
        days    <- dat$days
       depth    <- dat$depth
        dat     <- dat$data
        d_per_y <- 360
        #Get the data of the final year
        w       <- (nrow(dat)-d_per_y+1):nrow(dat)
        dat     <- dat[w,]
        s_dat   <- numeric(length(MLD))
        for (k in 1:length(MLD)){
             w    <- which(depth >= -abs(MLD[k]))
             xff  <- mean(as.numeric(dat[k,w]))
          s_dat[k] <- xff
        }

        lines(1:length(s_dat),s_dat, lwd=1,col=COLS[ii])

        #Plot remote sensing estimates of C:Chl ratio on the figure
        if (Var == 'C_Chl'){
           theta <- get_theta(Stn)
           points(theta$time,theta$data,pch=2)
        } else if (Var=='QN'){
          #Plot measurements of N:C on the figure
           QN=get_obs_MLD(Stn,'QN')
           points(QN$DOY,QN$dat,pch=2)
        }
        if (Var == 'mu' && Stn == 'S1') {
           legend('topright',Modnames,col=COLS,lty=1)
        }
    }
  }
}
mtext('Date of the year',side=1,line=1,outer=TRUE)

dev.off()
