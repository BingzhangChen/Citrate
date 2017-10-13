#This R script plots out relationships between mu and size variance at both stations:
#Read data:
source('~/Working/FlexEFT1D/Rscripts/plot_1D.R')
StnData <- function(Stn, model = 'NPZDcont'){
 filedir  <- paste0('~/working/FlexEFT1D/DRAM/',model,'/BOTH_TD/')
 get_SUR  <- function(x){
    DAT   <- get_Final_Data(filedir,Stn, x, finalyr = T, Dmax = -100)
    DAT   <- DAT$data
    DAT   <- DAT[, ncol(DAT)]
    return(DAT)
 }
    muAvg <- get_Final_Data(filedir,Stn,'muAvg', finalyr = T, Dmax = -100)
    Days  <- muAvg$days
    Season<- integer(length(Days))
    Season[Days <= 90] = 1
    Season[Days > 90 & Days <= 180] = 2
    Season[Days >180 & Days <= 270] = 3
    Season[Days > 270] = 4
    muAvg <- muAvg$data
    muAvg <- muAvg[, ncol(muAvg)]
    VAR   <- get_SUR('R_VAR')
    d2mu  <- get_SUR('d2mu')
    d4mu  <- get_SUR('d4mu')
    d2gd1 <- get_SUR('d2gd1')
    d2gd2 <- get_SUR('d2gd2')
    muN1  <- get_SUR('muN1')
    P     <- get_SUR('PHY1')

    #Diffusion:
    DY <- get_SUR('D_VAR')
    DX <- get_SUR('D_PMU')
    DP <- get_SUR('D_P1')
    Y  <- get_SUR('VAR')
    X  <- get_SUR('PMU')
    dVdt  <- DY/P + (2*X^2/P^3 - Y/P^2)*DP - 2*X/P^2*DX
    #dLdt  <- (D_PMU - PMU*D_P1)/PHY1
    #dVdt  <- (D_VAR - (PMU^2+VAR)*D_P1)/PHY1 -2*PMU*dLdt
    return(list(Days=Days,muAvg=muAvg, 
                VAR =VAR, d2mu=d2mu,  d4mu=d4mu,
                d2g1=d2gd1,d2g2=d2gd2, mu =muN1,
                dVdt=dVdt,
                Season=Season))
}
S1     = StnData('S1')
muAvg1 = S1$muAvg
  VAR1 = S1$VAR
Season = S1$Season
K2     = StnData('K2')
muAvg2 = K2$muAvg
  VAR2 = K2$VAR

pdf('Fig12_mu_Var.pdf', width = 5, height = 7)
op <- par(font.lab = 1,
            family ="serif", cex.axis=1.2, cex.lab=1.2,
            mar    = c(3.5,4,1.5,.5),
            mgp    = c(2.3,1,0),
            mfcol  = c(3,1),
            oma    = c(4,4,1,0)) 

plot(muAvg1, VAR1, xlim=range(c(muAvg1, muAvg2)), ylim=range(c(VAR1,VAR2)),
     xlab = bquote(µ[com] * ' ('*d^-1*')'),
     ylab = expression(paste("Size variance "*'(ln '*µm^3*')'^2)), 
     pch  = 16, col = Season)
points(muAvg2, VAR2, pch=1, col=Season)
legend('topright', legend=c('Winter', 'Spring', 'Summer', 'Fall','S1','K2'),
       pch=c(rep(16,5),1), col=c(1:4,1,1))
mtext('a)',adj=0, outer=F)

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

N_var <- D2mu+D2G1+D2G2+D4mu+VMU+S1$dVdt  #Net effect
plot(S1$Days, D2mu, type = 'l', xaxt = 'n',
     ylim = c(-0.2,0.5),
     xlab = '', 
     ylab = 'Effects on size diversity')
lines(S1$Days, D2G1, col='gray', lty = 2, lwd = .8)
lines(S1$Days, D2G2, col='gray0',lty = 2, lwd = .8)
lines(S1$Days, D4mu, col='tan',  lty = 2, lwd = .8)
lines(S1$Days, VMU,  col='green')
lines(S1$Days, S1$dVdt, col='blue', lty = 4)  #Contribution by diffusion
lines(S1$Days, N_var,col=6, lty = 4)
DVDt <- c(diff(VAR1),NA)  #Net changes in VAR based on modeled Variance
lines(S1$Days, DVDt, col=2, lty = 1)
abline(h=0)
legend('topleft', 
       legend=c('Competition','MIC grazing','MES grazing','d4µ/dL4'),
       col=c(1,'gray','gray0','tan'),
       lty=c(1,    2,    2,     2  ),
       lwd=c(1,.8,.8,.8))
legend('topright', 
       legend=c('Trait diffusion','Diffusion','Net effect','Net changes'),
       col=c('green','blue',6,2),
       lty=c(      1,  4,   4,1),
       lwd=c(1,1,1,1))

axis(1, at=seq(15,360,by=90),
     labels=c('Jan','Apr','Jul','Oct'))

mtext('b) S1',adj=0, outer=F)

#Plot K2:
D2mu <- VAR2*(VAR2 - 5*VTR)*K2$d2mu
D2G1 <- -VAR2**2*K2$d2g1
D2G2 <- -VAR2**2*K2$d2g2
D4mu <- VAR2**2*VTR*K2$d4mu
VMU  <- 2*VTR*K2$mu

N_var <- D2mu+D2G1+D2G2+D4mu+VMU+K2$dVdt  #Net effect
plot(K2$Days, D2mu, type = 'l', xaxt = 'n',
     ylim = range(c(D2mu,D2G1,D2G2,D4mu,VMU)),
     xlab = '', 
     ylab = 'Effects on size diversity')
lines(K2$Days, D2G1, col='gray', lty = 2, lwd = .8)
lines(K2$Days, D2G2, col='gray0',lty = 2, lwd = .8)
lines(K2$Days, D4mu, col='tan',  lty = 2, lwd = .8)
lines(K2$Days, VMU,  col='green')
lines(K2$Days, K2$dVdt, col='blue', lty = 4, lwd = 1)  #Contribution by diffusion
lines(K2$Days, N_var,col=6, lty = 4, lwd=1)
DVDt <- c(diff(VAR2),NA)  #Net changes in VAR based on modeled Variance
lines(K2$Days, DVDt, col=2, lty = 1, lwd=1)
abline(h=0)
axis(1, at=seq(15,360,by=90),
     labels=c('Jan','Apr','Jul','Oct'))

mtext('c) K2',adj=0, outer=F)

dev.off()
