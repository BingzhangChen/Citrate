#This R script plots out relationships between mu and size variance at both stations:
#Read data:
source('~/Working/FlexEFT1D/Rscripts/plot_1D.R')
StnData <- function(Stn, model = 'NPZDcont'){
 filedir  <- paste0('~/working/FlexEFT1D/DRAM/',model,'/BOTH_TD/')
    muAvg <- get_Final_Data(filedir,Stn,'muAvg', finalyr = T, Dmax = -100)
    Days  <- muAvg$days
    Season<- integer(length(Days))
    Season[Days <= 90] = 1
    Season[Days > 90 & Days <= 180] = 2
    Season[Days >180 & Days <= 270] = 3
    Season[Days > 270] = 4
    muAvg <- muAvg$data
    muAvg <- muAvg[, ncol(muAvg)]
    VAR   <- get_Final_Data(filedir,Stn,'R_VAR', finalyr = T, Dmax = -100)
    VAR   <- VAR$data
    VAR   <- VAR[, ncol(VAR)]
    d2mu  <- get_Final_Data(filedir,Stn,'d2mu', finalyr = T, Dmax = -100)
    d2mu  <- d2mu$data
    d2mu  <- d2mu[, ncol(d2mu)]
    d4mu  <- get_Final_Data(filedir,Stn,'d4mu', finalyr = T, Dmax = -100)
    d4mu  <- d4mu$data
    d4mu  <- d4mu[, ncol(d4mu)]
    d2gd1 <- get_Final_Data(filedir,Stn,'d2gd1', finalyr = T, Dmax = -100)
    d2gd1 <- d2gd1$data
    d2gd1 <- d2gd1[, ncol(d2gd1)]
    d2gd2 <- get_Final_Data(filedir,Stn,'d2gd2', finalyr = T, Dmax = -100)
    d2gd2 <- d2gd2$data
    d2gd2 <- d2gd2[, ncol(d2gd2)]
    muN1  <- get_Final_Data(filedir,Stn,'muN1', finalyr = T, Dmax = -100)
    muN1  <- muN1$data
    muN1  <- muN1[, ncol(muN1)]
    PHY1  <- get_Final_Data(filedir,Stn,'PHY1', finalyr = T, Dmax = -100)
    PHY1  <- PHY1$data
    PHY1  <- PHY1[, ncol(PHY1)]
    PMU   <- get_Final_Data(filedir,Stn,'R_PMU', finalyr = T, Dmax = -100)
    PMU  <- PMU$data
    PMU  <- PMU[, ncol(PMU)]

    #Diffusion:
    D_VAR  <- get_Final_Data(filedir,Stn,'D_VAR', finalyr = T, Dmax = -100)
    D_VAR  <- D_VAR$data
    D_VAR  <- D_VAR[, ncol(D_VAR)]
    D_PMU  <- get_Final_Data(filedir,Stn,'D_PMU', finalyr = T, Dmax = -100)
    D_PMU  <- D_PMU$data
    D_PMU  <- D_PMU[, ncol(D_PMU)]
    D_P1  <- get_Final_Data(filedir,Stn,'D_P1', finalyr = T, Dmax = -100)
    D_P1  <- D_P1$data
    D_P1  <- D_P1[, ncol(D_P1)]

    dLdt  <- (D_PMU - PMU*D_P1)/PHY1

    dVdt  <- (D_VAR - (PMU^2+VAR)*D_P1)/PHY1 -2*PMU*dLdt
    return(list(Days=Days,muAvg=muAvg, 
                VAR =VAR, d2mu=d2mu,  d4mu=d4mu,
                d2g1=d2gd1,d2g2=d2gd2, mu =muN1,
                
                Season=Season))
}
S1     = StnData('S1')
muAvg1 = S1$muAvg
  VAR1 = S1$VAR
Season = S1$Season
K2     = StnData('K2')
muAvg2 = K2$muAvg
  VAR2 = K2$VAR

pdf('Fig13_mu_Var.pdf', width = 9, height = 6, paper = 'a4')
op <- par(font.lab = 1,
            family ="serif", cex.axis=1.2, cex.lab=1.2,
            mar    = c(2,2,1.5,3.5),
            mgp    = c(2.3,1,0),
            mfcol  = c(3,1),
            oma    = c(4,4,1,0)) 


    plot(muAvg1, VAR1, xlim=range(c(muAvg1, muAvg2)), ylim=range(c(VAR1,VAR2)),
         xlab = bquote('µ ( '*d^-1*')'),
         ylab = expression(paste("Ln size variance "*'('*µm^3*')'^2*')')), 
         pch  = 16, col = Season)
    points(muAvg2, VAR2, pch=1, col=Season)
    legend('topright', legend=c('Winter', 'Spring', 'Summer', 'Fall'),
           pch=16, col=1:4)

    #Decompose different contributions to variance changes:
   #At S1:
   VTR <- 0.1
   #Contributions for d2mu:
   D2mu <- VAR1*(VAR1 - 5*VTR)*S1$d2mu
   D2G1 <- -VAR1**2*S1$d2g1
   D2G2 <- -VAR1**2*S1$d2g2
   D4mu <- VAR1**2*VTR*S1$d4mu
   VMU  <- 2*VTR*S1$mu

   #Contributions from diffusion:

   N_var<- D2mu+D2G1+D2G2+D4mu+VMU  #Net effect
   plot(S1$Days, D2mu, type = 'l', xaxt = 'n',
        ylim = range(c(D2mu,D2G1,D2G2,D4mu,VMU)),
        xlab = '', 
        ylab = 'Effects on size diversity')
   lines(S1$Days, D2G1, col=2)
   lines(S1$Days, D2G2, col=2, lty = 2)
   lines(S1$Days, D4mu, col=3, lty = 2)
   lines(S1$Days, VMU,  col=4, lty = 3)
   lines(S1$Days, N_var,col=5, lty = 4, lwd=2)
   axis(1, at=seq(15,360,by=90),
        labels=c('Jan','Apr','Jul','Oct'))
dev.off()
