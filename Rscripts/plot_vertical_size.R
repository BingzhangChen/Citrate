plot_v_size <- function(Stn,Models,BOTH=T){
   if (Stn == 'K2'){
      FigNo = 8
   }else if (Stn == 'S1'){
      FigNo = 9
   }

   fname <- paste0(Models,collapse='_')
   fname <- paste0('Fig_',FigNo,Stn,fname,'Vertical_mod_obs_size.pdf')
   Dmax  <- 150
   pdf(fname, width=2*3,height=8,paper='a4')
   
   op <- par(font.lab = 1,
                family ="serif",
                mar    = c(2,2,1.5,0.5),
                mgp    = c(2.3,1,0),
                oma    = c(4,4,4,0),
                mfcol  = c(4,4)   ) 
   
   j <- 0
   #Read observational data
   DAT = paste0(Stn,'_size_Perc.dat')
   DAT = read.table(DAT,header=T)
   
   for (K in 1:4){
     Var = paste0('CHLs',K)
     dat = DAT[,c(1,2,K+2)] 
     if (K == 1){
       Varname  <- 'Perc. > 10 µm'
     }else if (K == 2){
       Varname  <- 'Perc. 3~10 µm'
     }else if (K == 3){
       Varname  <- 'Perc. 1~3 µm'
     }else{
       Varname  <- 'Perc. < 1 µm'
     }
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
       dat_ <- dat_[dat_$Depth <= Dmax,]
       xmax <- max(dat_[,3])
       plot(dat_[,3],-dat_$Depth,xlim=c(0,1),ylim=c(-150,0),
            xlab=Varname,ylab='',pch=16,cex=.5,cex.lab=1.2)
       mtext(paste(letters[j],')',seasons[i]),adj=0)
   
        # if (i==1) mtext(Stn,adj=.7)
         ii = 0
         for (model in Models){
            ii    <- ii+1
            if (BOTH){
             DIR  <- paste0('~/working/FlexEFT1D/DRAM/',model,'/BOTH_TD/')
            } else{
             DIR  <- paste0('~/working/FlexEFT1D/DRAM/',model,'/',Stn,'/')
            }

            #Get modeled data
             Chl     = getData(DIR,Stn,Var)
             days    = Chl$days
            depth    = Chl$depth
             Chl     = Chl$data
             d_per_y = 360
             #Get the data of the final year
             w   <- (nrow(Chl)-d_per_y+1):nrow(Chl)
             Chl <- Chl[w,]
             Chl <- sapply(1:ncol(Chl), function(x)as.numeric(as.character(Chl[,x])))
             #Read model results:
             mdoy  <- 1:nrow(Chl)
             k     <- which(mdoy > DOYs[i] & mdoy <= DOYs[i+1])
             Chl_  <- Chl[k,]
             #Calculate quantiles:
             cff   <- sapply(1:ncol(Chl_),function(k)quantile(Chl_[,k], probs=c(0.025,0.5,0.975)))
             matlines(t(cff), depth, lty=c(2,1,2), lwd=c(.3,1.5,.3), col=COLS[ii])
         }
         
         if (K == 1 && i == 1 && length(Modnames) > 1){
            legend('topright',Modnames,col=COLS,lty=1)
         }
       }
     }   
     mtext('Depth (m)',side = 2, outer=TRUE, line=1)
    # mtext(Sys.time(), side = 3, outer=TRUE, line=2)
     mtext(paste('Fig.', FigNo,'. Model fittings to vertical profiles of four size fractions at',Stn),
           side=1,outer=T, line=2,adj=0)
   dev.off()
}
