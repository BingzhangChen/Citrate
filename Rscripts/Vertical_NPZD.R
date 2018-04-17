setwd('~/Working/FlexEFT1D/DRAM')  
#system('scp -r bzchen@chula.yes.jamstec.go.jp:/data16/bzchen/FlexEFT1D/DRAM_0.9/NPPZDD ./')
#system('scp -r bzchen@chula.yes.jamstec.go.jp:/data16/bzchen/FlexEFT1D/DRAM_0.9/NPZDcont/S1 ~/Working/FlexEFT1D/DRAM_0.9/NPZDcont/')
#system('scp -r bzchen@chula.yes.jamstec.go.jp:/data16/bzchen/FlexEFT1D/DRAM_0.9/NPZDFix ./')
#system('scp -r bzchen@chula.yes.jamstec.go.jp:/data16/bzchen/FlexEFT1D/DRAM_0.9/EFTsimple_sRun ./')

source('~/Working/FlexEFT1D/Rscripts/getData.R')
source('~/Working/FlexEFT1D/Rscripts/get_obs_MLD.R')
source('~/Working/FlexEFT1D/Rscripts/Interpolate_WOA.R')

#Estimate C:Chl from remote sensing:
Chl_C = readnc('Chl_C')  #A global Chl:C data on the surface
#Obtain theta from satellite:
theta = get_theta(Stn)
Mo    = theta$time
theta = theta$data
COLS     <- c(1,2)
Models   <- c('NPZDN2')
Modnames <- c('')
Stns     <- c('HOT')
Nstn     <- length(Stns)
#Get enspar:
model = 'NPZDN2'
Stn   = 'HOT'
DIR   = paste0('~/working/FlexEFT1D/DRAM/',model,'/',Stn,'/')
setwd(DIR)
np    = 2  #The number of CPUs for paralell computing
EnsLen= 1  #The number of ensembles
enspar= paste0(DIR,'enspar')
enspar= read.table(enspar, header=T)
enssig= paste0(DIR,'enssig')
enssig= read.table(enssig, header=T)

#Get bestpar:
best    = which.max(enspar$LogL)
bestpar = enspar[best,]

source('~/Working/FlexEFT1D/Rscripts/plot_vertical_NChlNPP.R')
variables=c('DIN','CHL','NPP','PON','DIP')  #Diazotroph underestimated
for (Stn in Stns){
    plot_v_n(Stn, Models, VARS = variables,BOTH=F)
}

source('~/Working/FlexEFT1D/Rscripts/plot_vertical_size.R')
#Plot vertical distributions of size
for (Stn in Stns){
    plot_v_size(Stn, Models)
}


#PLot out vertical distributions of mu, theta, and QN
for (Stn in Stns){
   fname <- paste(Stn,'_',model,'_Vertical_mu_theta.pdf',sep='')
   pdf(fname, width=2*3,height=8,paper='a4')
   #Average into 4 seasons
   DOYs    <- seq(0,360,length.out=5)
   seasons <- c('Winter','Spring','Summer','Fall')
   
   op <- par(font.lab  = 1,
                family ="serif",
                mar    = c(2,2,1.5,0.5),
                mgp    = c(2.3,1,0),
                oma    = c(4,4,0,0),
                mfcol  = c(4,3) ) 
   
   j <- 0
   Xmax <- 1
   for (Var in c('muN','The','QN_')){
     if (substr(Var,1,2) == 'QN'){
       Varname  <- 'N:C'
       Xmax     <-  .4
     }else if (substr(Var,1,3) == 'muN'){
       Varname  <- bquote('µ ( '*d^-1*')')
       Xmax     <- .5
     }else if (substr(Var,1,3) == 'The'){
       Varname  <- expression(paste("Chl:C "*' (gChl '*molC^-1*')'))
       Xmax     <- 0.6
     }else if (Var == 'PHY1'){
       Varname  <- expression(paste("PHY "*'(mmol '*m^-3*')'))
     }else{
       Varname  <- bquote(.(Var) ~ ' (mmol '*m^-3*')')
     }
   
     for (i in 1:4){
         if (i == 4) {
            par(mar=c(4,2,1.5,0.5))
         }else{
            par(mar=c(2,2,1.5,0.5))
         }
         ii <- 0
         for (model in Models){
            ii    <- ii+1
             DIR  <- paste0('~/Working/FlexEFT1D/DRAM/',model,'/',Stn,'/')
            #Get modeled data
            if (model %in% c('NPZDFix','NPZDN2')){
               NPHY <- 1
            }else if (model %in% c('NPZD2sp','NPPZDD')){
               NPHY <- 2
            }
            for (kk in 1:NPHY){
               cff   = paste0(Var,kk)
               Chl   = getData(DIR,Stn,cff)
               days  = Chl$days
               depth = Chl$depth
             Chl     = Chl$data
             d_per_y = 360
             #Get the data of the final year
             w       = (nrow(Chl)-d_per_y+1):nrow(Chl)
             Chl     = Chl[w,]
   
             #Read model results:
             mdoy  = 1:nrow(Chl)
             k     = which(mdoy > DOYs[i] & mdoy <= DOYs[i+1])
             Chl_  = Chl[k,]
             #Calculate mean:
             cff  <- apply(Chl_,2,mean)
            
             if (kk==1 && ii==1) {
                 plot(cff,depth, xlim=c(0,Xmax),
                      type='n', xlab=Varname,ylab='')
                 j = j + 1
                 mtext(paste(letters[j],')',seasons[i]),adj=0)
             }
             if (Var=='muN' && i == 1 && length(Modnames) > 1){
                 legend('bottomright',Modnames,col=COLS,lty=1)
             }

             lines(cff, depth, lty=kk, lwd=1.5, col=COLS[ii])
             if (Var == 'The'){
                #Calculate Avg. from the summer:
                k     = which(Mo > DOYs[i] & Mo <= DOYs[i+1])
                cff   = as.numeric(theta[k])
                points(mean(cff), -5)
                segments(mean(cff)-2*sd(cff),-5, mean(cff)+2*sd(cff),-5)
             }else if(Var == 'QN_'){
                #Add data of PON:POC
                #Read observational data
                 dat <- paste0('~/Working/FlexEFT1D/',Stn,'/',Stn,'_N2C.dat')
                 dat <- read.table(dat,header=T)
                 dat <- dat[dat$DOY >= DOYs[i] & dat$DOY <= DOYs[i+1],]
                 points(dat[,3], -dat$Depth, pch=16, cex=.5)
             }

           }  #kk
         }  #model
       }  #i (season)
     }  #Var 
   mtext('Depth (m)',side = 2, outer=TRUE, line=1)
   dev.off()
}

#PLot out vertical distributions of PHY, ZOO, DET
for (Stn in Stns){
   fname <- paste(Stn,'_Vertical_PZD.pdf',sep='')
   pdf(fname, width=2*3,height=8,paper='a4')
   #Average into 4 seasons
   DOYs    <- seq(0,360,length.out=5)
   seasons <- c('Winter','Spring','Summer','Fall')
   
   op <- par(font.lab = 1,
                family ="serif",
                mar    = c(2,2,1.5,0.5),
                mgp    = c(2.3,1,0),
                oma    = c(4,4,0,0),
                mfcol  = c(4,3) ) 
   
   j <- 0
   Xmax <- 1
   for (Var in c('PHY','ZOO','DET')){
     if (Var == 'PHY'){
       Varname  <- expression(paste("PHY "*'(mmol '*m^-3*')'))
       Xmax     <- .4
     }else{
       Varname  <- bquote(.(Var) ~ ' (mmol '*m^-3*')')
     }
   
     for (i in 1:4){
         if (i == 4) {
            par(mar=c(4,2,1.5,0.5))
         }else{
            par(mar=c(2,2,1.5,0.5))
         }
         ii <- 0
         for (model in Models){
            ii    <- ii+1
             DIR  <- paste0('~/Working/FlexEFT1D/DRAM_0.9/',model,'/',Stn,'/')
         
            #Get modeled data
            NPHY <- 1
            if (Var=='PHY' && model == 'NPZD2sp'){
               NPHY <- 2
            }
            for (kk in 1:NPHY){
              if (Var=='PHY'){
               cff   = paste0(Var,kk)
              }else{
               cff   = Var
              }
               Chl   = getData(DIR,Stn,cff)
               days  = Chl$days
               depth = Chl$depth
             Chl     = Chl$data
             d_per_y = 360
             #Get the data of the final year
             w       = (nrow(Chl)-d_per_y+1):nrow(Chl)
             Chl     = Chl[w,]
   
             #Read model results:
             mdoy  = 1:nrow(Chl)
             k     = which(mdoy > DOYs[i] & mdoy <= DOYs[i+1])
             Chl_  = Chl[k,]
             #Calculate mean:
             cff  <- apply(Chl_,2,mean)
            
             if (kk==1 && ii==1) {
                 plot(cff,depth, xlim=c(0,Xmax),
                      type='n', xlab=Varname,ylab='')
                 mtext(paste(letters[j],')',seasons[i]),adj=0)
                 j = j + 1
             }
             if (Var=='PHY' && i == 1 && length(Modnames) > 1){
                 legend('topright',Modnames,col=COLS,lty=1)
             }

             lines(cff, depth, lty=kk, lwd=1.5, col=COLS[ii])
           }  #kk
         }  #model
       }  #i (season)
     }  #Var 
   mtext('Depth (m)',side = 2, outer=TRUE, line=1)
   dev.off()
}
#PLot out vertical distributions of nutrient and light limitation:
for (Stn in Stns){
   fname <- paste(Stn,'_V_Lno3_SI.pdf',sep='')
   pdf(fname, width=2*2,height=8,paper='a4')
   #Average into 4 seasons
   DOYs    <- seq(0,360,length.out=5)
   seasons <- c('Winter','Spring','Summer','Fall')
   
   op <- par(font.lab = 1,
                family ="serif",
                mar    = c(2,2,1.5,0.5),
                mgp    = c(2.3,1,0),
                oma    = c(4,4,0,0),
                mfcol  = c(4,2) ) 
   
   j <- 0
   Xmax <- 1
   for (Var in c('Lno','SI_')){
     if (Var == 'Lno'){
         Varname = 'Nutrient limitation index'
     } else if (Var == 'SI_'){
         Varname = 'Light limitation index'
     }
     for (i in 1:4){
         if (i == 4) {
            par(mar=c(4,2,1.5,0.5))
         }else{
            par(mar=c(2,2,1.5,0.5))
         }
         ii <- 0
         for (model in Models){
            ii    <- ii+1
             DIR  <- paste0('~/Working/FlexEFT1D/DRAM_0.9/',model,'/',Stn,'/')
            #Get modeled data
            if (model == Models[1]){
               NPHY <- 1
            }else if (model == Models[2]){
               NPHY <- 2
            }
            for (kk in 1:NPHY){
               cff   = paste0(Var,kk)
               Chl   = getData(DIR,Stn,cff)
               days  = Chl$days
               depth = Chl$depth
             Chl     = Chl$data
             d_per_y = 360
             #Get the data of the final year
             w       = (nrow(Chl)-d_per_y+1):nrow(Chl)
             Chl     = Chl[w,]
   
             #Read model results:
             mdoy  = 1:nrow(Chl)
             k     = which(mdoy > DOYs[i] & mdoy <= DOYs[i+1])
             Chl_  = Chl[k,]
             #Calculate mean:
             cff   = apply(Chl_,2,mean)
            
             if (kk==1 && ii==1) {
                 plot(cff,depth, xlim=c(0,Xmax),
                      type='n', xlab=Varname,ylab='')
                 j = j + 1
                 mtext(paste(letters[j],')',seasons[i]),adj=0)
             }
             if (Var=='Lno' && i == 1 && length(Modnames) > 1){
                 legend('topright',Modnames,col=COLS,lty=1)
             }

             lines(cff, depth, lty=kk, lwd=1.5, col=COLS[ii])
           }  #kk
         }  #model
       }  #i (season)
     }  #Var 
   mtext('Depth (m)',side = 2, outer=TRUE, line=1)
   dev.off()
}

