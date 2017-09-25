#Plot station map:
source('~/Working/FlexEFT1D/Rscripts/Stnmap_Kv_S1K2.R')

Model <- 'NPZDcont'
DIR   <- paste0('~/Working/FlexEFT1D/DRAM/',Model,'/BOTH_TD/')
setwd(DIR)
source('~/Working/FlexEFT1D/Rscripts/plot_1D.R')
source('~/Working/FlexEFT1D/Rscripts/loglike_params.R')  #Plot time-evolution of Log-Likelihoods and parameters of size model (Fig. 4 and 5)

setwd('/Users/apple/Working/FlexEFT1D/DRAM/NPZDcont_sRun/BOTH')
Model   <- 'NPZDcont_sRun'

#Calculate total nitrogen and iron:
#source('~/Working/FlexEFT1D/Rscripts/plot_TN.R')

#Plot an example of four years to show seasonal cycle:
VARs    <- c('NO3','CHL_T','Fer','ZOO','R_PMU','R_VAR')
NVar    <- length(VARs)
pdffile <- paste0('fullyear_example_1D.pdf')

pdf(pdffile, width = 9, height = 6, paper = 'a4')
op <- par(font.lab = 1,
            family ="serif", cex.axis=1.2, cex.lab=1.2,
            mar    = c(2,2,1.5,3.5),
            mgp    = c(2.3,1,0),
            mfcol  = c(NVar,2),
            oma    = c(4,4,1,0)) 

j <- 0
for (Stn in c('K2','S1')){
   for (i in 1:NVar){
       VAR = VARs[i]
       plot_1D(VAR,Model,Stn,finalyr = F, Dmax = -180)
       j   = j + 1
       if (i == 1){
         mtext(paste0(letters[j],')',Stn),adj=0, outer=F)
       } else{
         mtext(paste0(letters[j],')'),adj=0, outer=F)
       }
   }
}
mtext('Depth (m)', side = 2, outer=TRUE, line=2)
mtext('Fig. 3 An example of modelled 4 year patterns at K2 and S1',side=1,outer=T, line=2,adj=0)
dev.off()


#One plot for one station (only final year):
VARs <- c('NO3','Fer','CHL_T',
          'R_PMU','R_VAR','muAvg','dmudl','d2mu','The1','QN_1')
Stns = c('K2','S1')

source('~/Working/FlexEFT1D/Rscripts/plot_stn_contour.R')
plot_stn(Stns, VARS, Model='NPZDcont', finalyr = T, Dmax = -150)
plot_stn('HOT',VARS, Model='NPZDcont_sRun', finalyr = T, Dmax = -150)

COLS     <- 2:3
Models   <- c('NPZDcont')
Modnames <- c('NPZDcont')
Stns     <- c('K2','S1')
Nstn     <- length(Stns)
DIR   <- paste0('~/Working/FlexEFT1D/DRAM/',Models,'/BOTH_TD/')
setwd(DIR)

source('~/Working/FlexEFT1D/Rscripts/plot_vertical_NChlNPP.R')
for (Stn in Stns){
    plot_v_n(Stn, Models, VARS=c('DIN','CHL','NPP','PON'))
}

#Plot for HOT:
Model    <- 'NPZDcont_sRun'
Modnames <- 'NPZDcont'
Stn   <- 'HOT'
DIR   <- paste0('~/Working/FlexEFT1D/DRAM/',Model,'/',Stn,'/')
setwd(DIR)

plot_v_n(Stn, Model, VARS=c('DIN','CHL','NPP','PON'), BOTH=F)

source('~/Working/FlexEFT1D/Rscripts/plot_vertical_size.R')
#Plot vertical distributions of size
for (Stn in Stns){
    plot_v_size(Stn, Models)
}

#Iron diagnostics:
TEMPBOL = function(Ea = 0.65, kb= 8.62E-5, tC, Tr=15) exp(-(Ea/kb)*(1./(273.15 + tC)-1./(273.15 + Tr)))

#source: dust deposition:
     Fe_N     <- 0.0265 
     Ez       <- .6
     model    <- 'NPZDcont'
     Stn      <- 'K2'
     Var      <- 'd2mu'
     filedir  <- paste0('~/working/FlexEFT1D/DRAM/',model,'/BOTH_TD/')
      data    <- getData(filedir,Stn,Var)
      days    <- data$days
     depth    <- data$depth
     VAR =  getData(filedir,Stn,'R_VAR')$data
     PMU =  getData(filedir,Stn,'R_PMU')$data
  #Find the depth corresponding to maximal VAR
     w   = which.max(as.numeric(VAR[460,]))


     d2mu=  data$data
     d4mu=  getData(filedir,Stn,'d4mu')$data
     d2g1=  getData(filedir,Stn,'d2gd1')$data
     d2g2=  getData(filedir,Stn,'d2gd2')$data
     PVAR=  getData(filedir,Stn,'VAR')$data
     PPMU=  getData(filedir,Stn,'PMU')$data
     mu  =  getData(filedir,Stn,'muN1')$data
     PHY =  getData(filedir,Stn,'PHY1')$data
   
     plot(days, d2mu[,w], type = 'l', ylim = c(-1E-3,1E-3))
     points(days,d2g1[,w], type = 'l', col=2)
     points(days,d2g2[,w], type = 'l', col=3)
     points(days,mu[,w], type = 'l', col=3)
     points(days,d4mu[,w]*.1, type = 'l', col=4)

     #Plot VAR: at Day 460, downward diffusion too large!!
     ND=430
     #Plot vertical distribution at day 460
     par(mfrow=c(3,2))
     plot(as.numeric(VAR[ND,]), depth, type = 'b')
     points(as.numeric(VAR[ND+1,]), depth, type = 'b',col=2)
     plot(as.numeric(PHY[ND,]), depth, type = 'b')
     plot(as.numeric(PMU[ND,]), depth, type = 'b')
     points(as.numeric(PMU[ND+1,]), depth, type = 'b',col=2)
     plot(as.numeric(D_VAR[ND,]), depth, type = 'b',xlim=c(-0.01,0.05))
     points(as.numeric(PVAR[ND,]), depth, type = 'b',col=2)
     plot(as.numeric(D_PMU[ND,]), depth, type = 'b',xlim=c(-0.01,0.05))
     points(as.numeric(PPMU[ND,]), depth, type = 'b',col=2)

     plot(days, VAR[,w], type = 'l')
     plot(days, PVAR[,w], type = 'l')
     dVAR=  getData(filedir,Stn,'dVAR')$data
     plot(days, dVAR[,w], type = 'l')
     plot(days, PHY[,w], type = 'l', ylim=c(0,.02))
     points(days, PVAR[,w], type = 'l',col=2)
     points(days, PVAR[,w+1], type = 'l',col=3)
     points(days, PVAR[,w]/(PHY[,w]**2), type = 'l',col=3)

     D_VAR=  getData(filedir,Stn,'D_VAR')$data
     D_PMU=  getData(filedir,Stn,'D_PMU')$data
     plot(days, D_VAR[,w], type = 'l',ylim=c(0,8E-5))
     points(days, PVAR[,w], type = 'l',col=2)
     points(days, D_VAR[,w+1], type = 'l',col=2)

     plot(days, D_PMU[,w], type = 'l',ylim=c(0,8E-5))
     points(days, PPMU[,w], type = 'l',col=2)
     points(days, D_PMU[,w+1], type = 'l',col=2)

    #Dissolved Fe:
      Fer     <- data$data
      #Calculate the net changes of Fer:
      M       <- ncol(Fer)
      NFer    <- matrix(NA, nr = nrow(Fer)-1, nc = ncol(Fer))
      for (i in 1:M){
         NFer[,i] <- diff(Fer[,i])
      }
      NFer <- rbind(NFer,0)
    #Dust deposition:
      depo    <- getData(filedir,Stn,'DtDep')$data

    #Phyto. uptake
      PHY     <- getData(filedir,Stn,'PHY_T')$data
      muAvg   <- getData(filedir,Stn,'muAvg')$data

      phyV    <- PHY*muAvg*Fe_N

    #Zooplankton excretion:
      PP_NZ   <- getData(filedir,Stn,'Z2N')$data*Fe_N

    #Iron scavenging
      Fescav  <- getData(filedir,Stn,'Fescv')$data 

    #Detritus regeneration:
      TEMP    <- getData(filedir,Stn,'Temp')$data
      DETFe   <- getData(filedir,Stn,'DETFe')$data
      PP_FeD  <- DETFe * 0.1 * TEMPBOL(Ea=Ez,tC=TEMP) 

    #Diffusion:
      Fediff  <- getData(filedir,Stn,'D_Fe')$data

    #Net budget:
     delFe  <- depo - phyV + PP_NZ - Fescav + PP_FeD + Fediff

  #Plot different contributions within MLD:
source('~/Working/FlexEFT1D/Rscripts/get_MLD.R')
source('~/Working/FlexEFT1D/Rscripts/Hz.R')
MLD = get_MLD(Stn)
Hz  = Hz(hmax = 250)

INT_MLD = function(dat){
  #dat: a matrix of nrow = length(days), nc = length(depth)
  #Determine the depth to integrate:
  V <- numeric(nrow(dat))
  for (i in 1:nrow(dat)){
         ff <- i%%length(MLD)

         if (ff == 0) ff <- length(MLD)

         w  <- which(depth >= -abs(MLD[ ff ]))
      V[i]  <- sum(as.numeric(dat[i,w]*Hz[w]))
  }
  return(V)
}

depo_    = INT_MLD(depo)
phyV_    = INT_MLD(phyV)
PP_NZ_   = INT_MLD(PP_NZ)
Fescav_  = INT_MLD(Fescav)
PP_FeD_  = INT_MLD(PP_FeD)
Fediff_  = INT_MLD(Fediff)

plot(days,   depo_, ylim=c(-0.2,0.2), type = 'l', col=2)
points(days,-phyV_,col='blue',type='l')
points(days,PP_NZ_,col='green',type='l')
points(days,-Fescav_,lty=2,type='l')
points(days, PP_FeD_,lty=2,col=2,type='l')
points(days, Fediff_,lty=1,col='orange',type='l')
legend('topright',legend = c('Deposition', 'Phy. uptake','Zoo. excretion','Scavenging','Det. rem.','Diffusion'),
    col=c(2,'blue','green',1,2,'orange'),lty=c(1,1,1,2,2,1))
#Compare The net budget with net changes of Fe:
plot(days, Fer[,M], type = 'l')


plot(days, delFe[,M], type = 'l')

#add net changes of Fe:
points(days, NFer[,M], col = 2, type = 'l')
