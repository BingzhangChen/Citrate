Model <- 'NPZDcont'
DIR   <- paste0('~/Working/FlexEFT1D/DRAM/',Model,'/BOTH_TD/')
setwd(DIR)
source('~/Working/FlexEFT1D/Rscripts/plot_1D.R')
source('~/Working/FlexEFT1D/Rscripts/loglike_size.R')  #Plot Log-Likelihoods of size model
source('~/Working/FlexEFT1D/Rscripts/params.R')  #Plot Log-Likelihoods of size model

setwd('/Users/apple/Working/FlexEFT1D/DRAM_0.9/NPZDcont_sRun/BOTH')
Model   <- 'NPZDcont_sRun'

#Calculate total nitrogen and iron:


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
VARs <- c('NO3','Fer','CHL_T','PHY1','R_PMU','R_VAR','muAvg','dmudl','d2mu','TD_VAR')
NVar <- length(VARs)
#TD_VAR is the contribution of "trait diffusion" to changes of size variance

for (Stn in c('K2','S1')){
  if (Stn == 'K2'){
      FigNo = 9
   }else if (Stn == 'S1'){
      FigNo = 10
   }
   pdffile <- paste0('Fig_',FigNo,Stn,'_1D.pdf')
   pdf(pdffile, width = 8, height = 10, paper = 'a4')
   op <- par(font.lab = 1,
               family ="serif", cex.axis=1.2, cex.lab=1.2,
               mar    = c(2,2,1.5,3.5),
               mgp    = c(2.3,1,0),
               mfrow  = c(ceiling(NVar/2),2),
               oma    = c(4,4,1,0)) 
   for (i in 1:NVar){
       VAR = VARs[i]
       plot_1D(VAR,Model,Stn,finalyr = F, Dmax = -500)
       mtext(paste0(letters[i],')'),adj=0, outer=F)
   }
   mtext('Depth (m)', side = 2, outer=TRUE, line=2)
   mtext(paste('Fig.', FigNo,'. Modelled seasonal patterns at',Stn),side=1,outer=T, line=2,adj=0)
   dev.off()
}

COLS     <- 2:3
Models   <- c('NPZDcont')
Modnames <- c('NPZDcont')
Stns     <- c('K2','S1')
Nstn     <- length(Stns)
source('~/Working/FlexEFT1D/Rscripts/plot_vertical_NChlNPP.R')
for (Stn in Stns){
    plot_v_n(Stn, Models)
}
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
     Var      <- 'Fer'
     filedir  <- paste0('~/working/FlexEFT1D/DRAM/',model,'/BOTH_TD/')
      data    <- getData(filedir,Stn,Var)
      days    <- data$days
     depth    <- data$depth
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
