setwd('~/Working/FlexEFT1D/DRAM') 
source('~/Working/FlexEFT1D/Rscripts/getData.R')
source('~/Working/FlexEFT1D/Rscripts/get_obs_MLD.R')
source('~/Working/FlexEFT1D/Rscripts/Interpolate_WOA.R')

#Obtain theta from satellite:
Chl_C = readnc('Chl_C', ROMS=F)
theta = get_theta('HOT')
Mo    = theta$time
theta = theta$data
COLS     <- c(3,2)
Models   <- c('NPZDFix_sRun','EFTsimple_sRun')
Modnames <- c('MONOD', 'PAHLOW')
Stns     <- c('HOT','S1')
Nstn     <- length(Stns)

#PLot out vertical distributions of QN
fname <- paste('Fig10_Vertical_QN.pdf',sep='')
pdf(fname, width=2*2,height=6,paper='a4')
#Average into 4 seasons
DOYs    <- seq(0,360,length.out=5)
seasons <- c('Winter','Spring','Summer','Fall')

op <- par(font.lab  = 1,
             family ="serif",
             mar    = c(2,2,1.5,0.5),
             cex.lab= 1.2,
            cex.axis= 1.2,
             mgp    = c(2.3,1,0),
             oma    = c(4,4,0,0),
             mfcol  = c(4,2) ) 
j <- 0
for (Stn in Stns){
   Xmax <- 1
   for (Var in c('QN_')){
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
             NPHY <- 1
            for (kk in 1:NPHY){
               cff   = paste0(Var,kk)
               Chl   = getData(DIR,Stn,cff)
               days  = Chl$days
               depth = Chl$depth
             Chl     = Chl$data
             Dmax    = -200
             Kdep    = which(depth >= Dmax)
            depth    = depth[Kdep]
             Chl     = Chl[,Kdep]

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
                 if (i == 1){
                     mtext(Stn,adj=1)
                     if (Stn == 'S1') legend('topright', Modnames,lty=1, col=COLS)
                 }
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
}
mtext('Depth (m)',side = 2, outer=TRUE, line=1)
dev.off()
