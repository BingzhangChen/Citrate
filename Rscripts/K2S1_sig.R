#Read size fractionated data:

#Read observational data
X = c(1, 3, 10)
X = log(X)

setwd('~/Working/FlexEFT1D/DRAM_0.9')

get_musig <- function(Stn){
  #Get total Chl
  TCHL = paste0('~/Working/FlexEFT1D/Forcing/',Stn,'_size.dat')
  TCHL = read.table(TCHL,header=T)
  TCHL1= TCHL[, -c(1,2)]
  TCHL$TChl = apply(TCHL1, 1, sum)
  
  for (i in 3:6){
      TCHL[,i] = TCHL[,i]/TCHL$TChl
  }
  TCHL$R_PMU = NA
  TCHL$R_VAR = NA
  DAT        = TCHL
  DAT        = DAT[DAT$TChl > 0, ]
  
  for (i in 1:nrow(DAT)){
     if (any(DAT[i, 3:6] <= 0))  next 
     #Calculate pnorm (percentage < threshold):
     Y    = numeric(length(X))
     Y[1] = DAT[i, 6]
     Y[2] = DAT[i, 6] + DAT[i, 5]
     Y[3] = Y[2]      + DAT[i, 4]
     Yq   = qnorm(Y)
     LM   = lm(X ~ Yq)
     LM   = as.numeric(coef(LM))
     DAT$R_PMU[i]  = exp(LM[1])  #Real ESD, NOT logESD
     DAT$R_VAR[i]  = LM[2]
     #Plot an example:
     #if (i == 1){
     #   X1 = exp(X)
     #   pdf(paste0(Stn,'Fit_example.pdf'), width=4,height=4,page='a4')
     #   plot(X1, Y, ylim = c(0,1), 
     #               xlab ='ESD (µm)',
     #               ylab ='Cumulative probability')
     #   newx = seq(0.5, 25, by = 0.1)
     #   newy = pnorm(log(newx), LM[1], LM[2])
     #   lines(newx,newy)
     #   dev.off()
     #}
  }
  #Calculate Shannon index:
  DAT$SWH  <- sapply(1:nrow(DAT), function(i)shannon(DAT[i,3:6]))
  return(DAT)
}

#Plotting:
musig <- function(Stn){

  DAT <- get_musig(Stn)
   
  #Plot TChl and percentages:
  plot(DAT$TChl, DAT$SIZE10, pch=16, col=3,
       xlab = 'Total Chl (µg/L)',
       ylab = '> 10 µm %')

  source('~/Working/FlexEFT1D/Rscripts/Varname.R')
  #Average into 4 seasons
  DOYs    <- seq(0,360,length.out=5)
  seasons <- c('Winter','Spring','Summer','Fall')
  
  fname <- paste0(Stn,'_obs_size.pdf')
  pdf(fname, width=2*2,height=8,paper='a4')
  op <-  par(font.lab = 1,
               family ="serif",
               mar    = c(2,2,1.5,0.5),
               mgp    = c(2.3,1,0),
               oma    = c(4,4,4,0),
               mfcol  = c(4,2)   ) 
  j <- 0
  for (VAR in c('R_PMU', 'R_VAR')){
     VARNAME <- Varname(VAR)
     xmax    <- quantile(DAT[,VAR], probs=0.99, na.rm=T)
     for (i in 1:4){
          if (i == 4) {
             par(mar=c(4,2,1.5,0.5))
          }else{
             par(mar=c(2,2,1.5,0.5))
          }
          j    <- j+1
          dat_ <- DAT[DAT$DOY > DOYs[i] & DAT$DOY <= DOYs[i+1],]
          plot(dat_[,VAR],-dat_$Depth,xlim=c(0,xmax), ylim=c(-150,0),
               xlab=VARNAME,ylab='',pch=16,cex=.5,cex.lab=1.2)
          mtext(paste(letters[j],')',seasons[i]),adj=0)
     }
  }
  dev.off()
} 

K2 <- get_musig('K2')
S1 <- get_musig('S1')

#Plot out Shannon index vs. size diversity for K2 and S1:
pwd1   = 'npacific'
setwd(paste0('~/Roms_tools/Run/',pwd1))

pdf('Shannon_size_diversity_K2S1.pdf', width = 5, height = 5)
op <- par(font.lab  = 1,
             family ="serif",
             mar    = c(4,4,.5,0.5),
             mgp    = c(2.3,1,0),
             oma    = c(4,4,0,0),
             mfcol  = c(1,1), pch=1, 
             cex.lab=1.2,cex.axis=1.2 ) 

plot(K2$SWH, K2$R_VAR, pch=16, xlim=c(0.8,1.5), ylim=c(0,3),
           xlab = 'Shannon index',
           ylab = expression("Size diversity ((Ln "*µm^3*')'^2*")"))
points(S1$SWH, S1$R_VAR)
legend('topleft',legend=c('K2','S1'),
        pch=c(16,1))

dev.off()

