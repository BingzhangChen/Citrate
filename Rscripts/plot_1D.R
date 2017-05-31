library(plot3D)
source('~/Working/FlexEFT1D/Rscripts/getData.R')
plot_1D <- function(Var,model,Stn,title='',ZLIM=NULL, finalyr = F, BOTH = T, Dmax = -200){
    if (BOTH){
     filedir  <- paste0('~/working/FlexEFT1D/DRAM_0.9/',model,'/BOTH/')
    } else{
     filedir  <- paste0('~/working/FlexEFT1D/DRAM_0.9/',model,'/',Stn,'/')
    }
      data    <- getData(filedir,Stn,Var)
      days    <- data$days
     depth    <- data$depth
      data    <- data$data
      Kdep    <- which(depth >= Dmax)
     depth    <- depth[Kdep]
      data    <- data[,Kdep]
     if (finalyr) {
        #Take the final year   
      fday <- days[length(days)]  #Final day
      cff  <- which( (days > fday-360) & (days <= fday))
      days <- days[cff]%%360
      days[days==0] <- 360
      data <- data[cff,]
    }
    ZLIM = as.double(quantile(unlist(data),probs=c(0.01,0.95)))
    data[data < ZLIM[1]] <- ZLIM[1]
    data[data > ZLIM[2]] <- ZLIM[2]
    image2D(as.matrix(data), x=days, y=-depth,
           xlab='',ylab='',main=title,adj=0,xaxt='n',cex.axis=1.2,cex.lab=1.2)
    #Add Variable name and units:
    if (Var == 'NO3'){
      Varname  <- bquote('TIN (mmol '*m^-3*')')
    }else if (Var == 'Fer'){
      Varname  <- expression(paste("Fer "*'(µmol '*m^-3*')'))
    }else if (Var == 'CHL_T'){
      Varname  <- expression(paste("Chl "*'(mg '*m^-3*')'))
    }else if (Var == 'NPP'){
      Varname  <- bquote(.(Var)~' (mg C '*m^-3*' '*d^-1*')')
    }else if (substr(Var,1,2) == 'mu'){
      K        <- substr(Var,4,4)
#      Varname  <- bquote('µ'~ .(K) ~' ( '*d^-1*')')
      Varname  <- bquote('µ ( '*d^-1*')')
    }else if (Var == 'The1'){
      Varname  <- expression(paste("Chl:C "*' (gChl '*molC^-1*')'))
    }else if (Var == 'QN_1'){
      Varname  <- expression(paste("N:C "*' (molN '*molC^-1*')'))
    }else if (Var == 'PHY1'){
      Varname  <- expression(paste("PHY "*'(mmol '*m^-3*')'))
    }else if (Var == 'R_PMU'){
      Varname  <- expression(paste("Mean size "*'( '*µm*')'))
    }else if (Var == 'R_VAR'){
      Varname  <- expression(paste("Ln size variance "*'('*µm^3*')'^2*')'))
    }else if (Var == 'dmudl'){
      Varname  <- expression(paste("dµ/dL ("*d^-1*' (ln '*µm^3*')'^-1*' )'))
    }else if (Var == 'd2mu'){
      Varname  <- expression(paste("d2µ/dL2 ("*d^-2*' (ln '*µm^3*')'^-2*')'))
    }else if (Var == 'TD_VAR'){
      Varname  <- expression(paste("Trait diffusion ("*d^-2*' (ln '*µm^3*')'^-2*')'))
    }else if (Var == 'd2gdl'){
      Varname  <- expression(paste("d2g/dL2 ("*d^-2*' (ln '*µm^3*')'^-2))
    }else{
      Varname  <- bquote(.(Var) ~ ' (mmol '*m^-3*')')
    }

    mtext(Varname, side = 3, adj=1)

    #Plot nutricline
    #if(Var == 'NO3') {
    #  NO3cln <- numeric(nrow(data))
    #  for (i in 1:nrow(data)){
    #      x  <- data[i,]
    #      z  <- which.min(x>1)
    #  NO3cln[i] <- depth[z]
    #  }
    #  lines(days,-NO3cln,lwd=2,col='green')
    #}

   if (!finalyr) {
      axis(1, at=seq(0,max(days),by=360))
      abline(  v=seq(0,max(days),by=360),lty=3)
   }else {
      axis(1, at=seq(15,360,by=90), labels=c('Jan','Apr','Jul','Oct'))
   }
}


