#This R script plots out relationships between mu and size variance at both stations:
#Read data:
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
    return(list(muAvg=muAvg, VAR=VAR, d2mu=d2mu,d4mu=d4mu,
                d2g1 =d2gd1,d2g2=d2gd2,mu =muN1,Season=Season))
}
S1     = StnData('S1')
muAvg1 = S1$muAvg
  VAR1 = S1$VAR
Season = S1$Season
K2     = StnData('K2')
muAvg2 = K2$muAvg
  VAR2 = K2$VAR

    plot(muAvg1, VAR1, xlim=range(c(muAvg1, muAvg2)), ylim=range(c(VAR1,VAR2)),
         xlab = bquote('µ ( '*d^-1*')'),
         ylab = expression(paste("Ln size variance "*'('*µm^3*')'^2*')')), 
         pch  = 16, col = Season)
    points(muAvg2, VAR2, pch=1, col=Season)
    legend('topright', legend=c('Winter', 'Spring', 'Summer', 'Fall'),
           pch=16, col=1:4)

    #Decompose different contributions to variance changes:
   #At S1:

