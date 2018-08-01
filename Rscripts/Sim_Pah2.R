# Maximal growth rate (Umax) and half-saturation constant size dependent
source('~/Working/FlexEFT1D/Rscripts/mu_EFT.R')
nls.control(maxiter = 80, tol = 1e-03,
    minFactor = 1/(2048*5),
    printEval = FALSE, warnOnly = TRUE)

param_DRAM=TRUE

if (param_DRAM) {
   #Read best parameters from DRAM output:
   params = read.table('~/Working/FlexEFT1D/DRAM_0.8/EFTsimple/S1/enspar',header=T)
   params = params[which.max(params$LogL),]
   sigma  = read.table('~/Working/FlexEFT1D/DRAM_0.8/EFTsimple/S1/enssig',header=T)
}else{
   NPar   = 5
   params = data.frame(matrix(NA,nr=1,nc=NPar))
   colnames(params)=c("mu0hat","aI0" ,"A0N", "Q0N", "wDET")
   params$mu0hat=5
   params$aI0   =0.59
   params$A0N   =0.23
   params$Q0N   =0.04
   params$wDET  =1
}

#Compare relationships between growth vs. nutrient and growth vs. light
#for two models
#Plot out:

NO3   <- seq(1E-4, 20, 0.05)
PAR   <- seq(.01,400,1)
dat1  <- expand.grid(NO3=NO3,PAR=PAR)
N     <- nrow(dat1)

#Fit the growth curves of flexible model with Monod equation:
#Light is saturating (400 W m-2)
#Use Monod model to fit Pahlow model:

#Flexible:
MU1  <- function(NO3,par_,temp=15)mu_EFT(temp=temp,NO3 = NO3, par_ = par_,aI=params$aI0,Q0N=params$Q0N,mu0=params$mu0hat,A0N=params$A0N)
MU2  <- sapply(1:N,function(x)MU1(par_=dat1[x,'PAR'],NO3=dat1[x,'NO3'])$mu) 

dat1$mu = MU2
Fit2 <- nls( mu ~ rmax * NO3/(NO3 + KN)*(1-exp(-aI0*PAR/rmax)),
                data  = dat1, 
                start = list(rmax = 3, KN = 0.1, aI0 = 0.05),
                algorithm = "port" )
#rmax=coef(Fit2)[1]
#  KN=coef(Fit2)[2]
#aI0 =coef(Fit2)[3]
#the_EFT <- sapply(1:N,function(x)MU1(par_=dat1[x,'PAR'],NO3=dat1[x,'NO3'])$Theta) 
#themax  <- max(the_EFT)

#Consistent with the ms
rmax=2.44
KN  =1.33
aI0 =0.11
themax=.61

#MONOD model:
MONOD <- function(NO3,par_,temp=15){
         mu_fix(NO3=NO3,par_=par_,mu0=rmax,KN0=KN,Temp=temp,
                Q0N=params$Q0N,aI0C=aI0, thetamax=themax)
}

MU_sim2 <- predict(Fit2,newdata=data.frame(NO3=NO3,PAR=max(dat1$PAR)))

#Fit3 <- nls( mu ~ rmax * NO3/(NO3 + KN)*(1-exp(-aI0*PAR/rmax/(NO3/(NO3+KN)))),
#                data  = dat1, 
#                start = list(rmax = 3, KN = 0.1, aI0 = 0.05),
#                algorithm = "port" )
#rmax=coef(Fit3)[1]
#  KN=coef(Fit3)[2]
#aI0 =coef(Fit3)[3]

pdf('~/Working/FlexEFT1D/DRAM/N_PAR_EFT_FIX_One.pdf',
                          width=3.5, height=8,paper='a4')
op <- par( font.lab = 1,
             family ="serif",
             mar    = c(4,4,2,0.3),
             mgp    = c(2.3,1,0),
             mfrow  = c(2,1),
             cex.lab= 1.2, cex=1,
             lwd    = 1.5,
             pch    = 16,
             cex.axis=1) 

MU2 = dat1[dat1$PAR>=399,]$mu

plot(NO3, MU2, cex=.3, col=2,
               xlab = expression(paste('Nitrogen (µmol '*L^-1*')')),
               ylab = expression(paste('Growth rate ('*d^-1*')')))
points(NO3,MU_sim2,type='l',col=3)
mtext('a) Growth ~ nutrient',adj=0)

#plot P-I relationships under nutrient replete conditions:

MU_EFT <- dat1[dat1$NO3==max(dat1$NO3),]$mu
     
par(mar=c(4,4,1,0))
plot(PAR, MU_EFT, cex=.3, col=2,
     xlab = expression(paste('PAR (W '*m^-2*')')),
     ylab = expression(paste('Growth rate ('*d^-1*')'))
)
mtext('b) Growth ~ light',adj=0)
MU_sim <- predict(Fit2,newdata=data.frame(NO3=max(dat1$NO3),PAR=PAR))
points(PAR,MU_sim,type='l',col=3)
legend('bottomright',pch=c(16,NA),lty=c(NA,1),col=c(2,3),
        legend=c('PAHLOW','MONOD'))
dev.off()

#Change parameters in the single run namelist:
setwd('~/Working/FlexEFT1D/DRAM/NPZDFix_sRun/HOT_A2/')
fname='param.nml'
mz   =0.15
gmax =1
if(file.exists(fname)) file.remove(fname)
file.create(fname)
fw=file(fname,open='wt')
writeLines('&parameters',con=fw) #write header to file
txt=paste0('mu0hat=',rmax,',')
writeLines(txt,con=fw) #write to file
txt=paste0('KN=',KN,',')
writeLines(txt,con=fw) #write to file
txt=paste0('wDET=',params$wDET,',')
writeLines(txt,con=fw) #write to file
txt=paste0('aI0_C=',aI0,',')
writeLines(txt,con=fw) #write to file
txt=paste0('Q0N=',params$Q0N,',')
writeLines(txt,con=fw) #write to file
txt=paste0('mz=',mz,',')
writeLines(txt,con=fw) #write to file
txt=paste0('gmax=',gmax,',')
writeLines(txt,con=fw) #write to file
writeLines('/',con=fw)
close(fw)
system('./run')
system('./NPZDFix > Out')

setwd('~/Working/FlexEFT1D/DRAM/EFTsimple_sRun/HOT/')
fname='param.nml'
if(file.exists(fname)) file.remove(fname)
file.create(fname)
fw=file(fname,open='wt')
writeLines('&parameters',con=fw) #write header to file
txt=paste0('mu0hat=',params$mu0hat,',')
writeLines(txt,con=fw) #write to file
txt=paste0('A0N=', params$A0N,',')
writeLines(txt,con=fw) #write to file
txt=paste0('wDET=',params$wDET,',')
writeLines(txt,con=fw) #write to file
txt=paste0('Q0N=',params$Q0N,',')
writeLines(txt,con=fw) #write to file
txt=paste0('aI0=',params$aI0,',')
writeLines(txt,con=fw) #write to file
txt=paste0('mz=',mz,',')
writeLines(txt,con=fw) #write to file
txt=paste0('gmax=',gmax,',')
writeLines(txt,con=fw) #write to file
writeLines('/',con=fw)
close(fw)
system('./run')
system('./EFTsimple > Out')

#Check the interactions between light and nutrient:
pdf('~/Working/FlexEFT1D/DRAM/Light_NO3_interaction_rev.pdf',
                    width=8, height=9,paper='a4')
mumax = 3.5

op <- par( font.lab = 1,
             family ="serif",
             mar    = c(3.5,4,1,0.2),
             mgp    = c(2.3,1,0),
             mfrow  = c(3,3),
             cex.lab= 1.2, cex=1,
             lwd    = 1.5, las=1,
             cex.axis=1) 

PAR     <- c(5, 30, 400)
NO3     <- seq(1E-4, 20, 0.01)
MU_FIX  <- matrix(NA, nr = length(NO3), nc = length(PAR) )
MU_EFT  <- MU_FIX

THE_FIX <- MU_FIX
THE_EFT <- MU_FIX
 QN_FIX <- MU_FIX
 QN_EFT <- MU_FIX
for (i in 1:length(PAR)){

    MU_FIX[,i]  <- MONOD(par_=PAR[i],  NO3=NO3)$mu 
    MU_EFT[,i]  <- sapply(1:length(NO3), function(x)MU1(par_=PAR[i],NO3=NO3[x])$mu)
    if(i==1) {
       plot(NO3, MU_FIX[,i], type = 'l', ylim=c(0,mumax), lty=i,col='green',
            xlab = expression(paste('Nitrogen (µmol '*L^-1*')')),
            ylab = expression(paste('Growth rate ('*d^-1*')')))
    }else{
       points(NO3, MU_FIX[,i], type = 'l',lty=i,col=3)
    }
    points(NO3,MU_EFT[,i],type='l', lty=i,col=2)
}
legend('topleft', legend=paste0('PAR=',PAR), lty=1:length(PAR), cex=.8 )

mtext('a',adj=0)

#Chl:C ratio:
for (i in 1:length(PAR)){

    THE_FIX[,i]  <- MONOD(par_=PAR[i],  NO3=NO3)$Theta 
    THE_EFT[,i]  <- sapply(1:length(NO3),function(x)MU1(par_=PAR[i],NO3=NO3[x])$Theta) 
    if(i==1) {
       plot(NO3, THE_FIX[,i], type = 'l', ylim=c(0,themax), lty=i,col=3,
            xlab = expression(paste('Nitrogen (µmol '*L^-1*')')),
            ylab = expression(paste('Chl:C (molC '*gChl^-1*')')))
    }else{
       points(NO3, THE_FIX[,i], type = 'l',lty=i,col=3)
    }
    points(NO3,THE_EFT[,i],type='l', lty=i,col=2)
}

mtext('b',adj=0)

#QN ratio:
for (i in 1:length(PAR)){

    QN_FIX[,i]  <- MONOD(par_=PAR[i],  NO3=NO3)$QN 
    QN_EFT[,i]  <- sapply(1:length(NO3),function(x) MU1(par_=PAR[i],NO3=NO3[x])$QN) 
    if(i==1) {
       plot(NO3, QN_FIX[,i], type = 'l', 
            ylim = c(params$Q0N,4.1*params$Q0N), lty=i,col=3,
            xlab = expression(paste('Nitrogen (µmol '*L^-1*')')),
            ylab = expression(paste('N:C (mol '*mol^-1*')')))
    }else{
       points(NO3, QN_FIX[,i], type = 'l',lty=i,col=3)
    }
    points(NO3,QN_EFT[,i],type='l', lty=i,col=2)
}

mtext('c',adj=0)
legend('bottomright',legend=c('Monod','Pahlow'),lty=1,col=c(3,2), cex=.8 )
NO3     <- c(.1, 1, 10)
PAR     <- seq(0.1, 300, .1)
N       <- length(PAR)
MU_FIX  <- matrix(NA, nr = length(PAR), nc = length(NO3) )
MU_EFT  <- MU_FIX
THE_FIX <- MU_FIX
THE_EFT <- MU_FIX
 QN_FIX <- MU_FIX
 QN_EFT <- MU_FIX

par(mar    = c(4,4,1.5,0.2))
for (i in 1:length(NO3)){

    MU_FIX[,i]  <- MONOD(par_=PAR,  NO3=NO3[i])$mu 
    MU_EFT[,i]  <- sapply(1:N, function(x)MU1(par_=PAR[x],NO3=NO3[i])$mu) 
    if(i==1) {
       plot(PAR, MU_FIX[,i], type = 'l', ylim=c(0,mumax), lty=i,col=3,
            xlab = expression(paste('PAR (W '*m^-2*')')),
            ylab = expression(paste('Growth rate ('*d^-1*')')))
    }else{
       points(PAR, MU_FIX[,i], type = 'l',lty=i,col=3)
    }
    points(PAR,MU_EFT[,i],type='l', lty=i,col=2)
}
mtext('d',adj=0)
legend('topleft',legend=paste0('DIN=',NO3),lty=1:length(NO3), cex=.8 )

for (i in 1:length(NO3)){

    THE_FIX[,i]  <- MONOD(par_=PAR,  NO3=NO3[i])$Theta 
    THE_EFT[,i]  <- sapply(1:N, function(x)MU1(par_=PAR[x],NO3=NO3[i])$Theta) 
    if(i==1) {
       plot(PAR, THE_FIX[,i], type = 'l', ylim=c(0,themax), lty=i,,col=3,
            xlab = expression(paste('PAR (W '*m^-2*')')),
            ylab = expression(paste('Chl:C (molC '*gChl^-1*')')))
    }else{
       points(PAR, THE_FIX[,i], type = 'l',lty=i,col=3)
    }
    points(PAR,  THE_EFT[,i],type='l', lty=i,col=2)
}
mtext('e',adj=0)

for (i in 1:length(NO3)){

    QN_FIX[,i]  <- MONOD(par_=PAR,  NO3=NO3[i])$QN 
    QN_EFT[,i]  <- sapply(1:N, function(x)MU1(par_=PAR[x],NO3=NO3[i])$QN) 

    if(i==1) {
       plot(PAR, QN_FIX[,i], type = 'l', 
            ylim = c(params$Q0N,4.1*params$Q0N),
             lty = i, col = 3,
            xlab = expression(paste('PAR (W '*m^-2*')')),
            ylab = expression(paste('N:C (mol '*mol^-1*')')))
    }else{
       points(PAR, QN_FIX[,i], type = 'l',lty=i,col=3)
    }
    points(PAR,  QN_EFT[,i],type='l', lty=i,col=2)
}
mtext('f',adj=0)

#Plot mu, theta, and QN at different temperatures
#4 combinations of light and nutrients
SST = 1:30
NT  = length(SST)
env = expand.grid(PAR=c(5,400), DIN=c(0.1,10))
N   = nrow(env)
MU_FIX = matrix(NA, nr = NT, nc = N)
MU_EFT = MU_FIX
THE_FIX <- MU_FIX
THE_EFT <- MU_FIX
 QN_FIX <- MU_FIX
 QN_EFT <- MU_FIX

for (i in 1:N){
    MU_FIX[,i]  <- MONOD(par_=env[i,]$PAR, 
                          NO3=env[i,]$DIN,
                         temp=SST)$mu 

    MU_EFT[,i]  <- sapply(1:NT, function(x)MU1(temp=SST[x], 
                                               par_=env[i,]$PAR,
                                               NO3 =env[i,]$DIN)$mu) 

    if(i==1) {
       plot(SST, MU_FIX[,i], type = 'l', ylim=c(0,mumax), lty=i,col=3,
            xlab = 'Temperature (ºC)',
            ylab = expression(paste('Growth rate ('*d^-1*')')))
    }else{
       points(SST, MU_FIX[,i], type = 'l',lty=i,col=3)
    }
    points(SST,MU_EFT[,i],type='l', lty=i,col=2)
}
mtext('g',adj=0)
legend('topleft',legend=paste0('DIN=',env$DIN,', PAR=',env$PAR),lty=1:nrow(env), cex=.8 )

for (i in 1:N){
    THE_FIX[,i]  <- MONOD(par_=env[i,]$PAR, 
                           NO3=env[i,]$DIN,
                          temp=SST)$Theta 

    THE_EFT[,i]  <- sapply(1:NT, function(x)MU1(temp=SST[x], 
                                               par_=env[i,]$PAR,
                                               NO3 =env[i,]$DIN)$Theta) 

    if(i==1) {
       plot(SST, THE_FIX[,i], type = 'l', ylim=c(0,0.65), lty=i,col=3,
            xlab = 'Temperature (ºC)',
            ylab = expression(paste('Chl:C (molC '*gChl^-1*')')))
    }else{
       points(SST, THE_FIX[,i], type = 'l',lty=i,col=3)
    }
    points(SST,THE_EFT[,i],type='l', lty=i,col=2)
}
mtext('h',adj=0)

for (i in 1:N){
    QN_FIX[,i]  <- MONOD(par_=env[i,]$PAR, 
                           NO3=env[i,]$DIN,
                          temp=SST)$QN 

    QN_EFT[,i]  <- sapply(1:NT, function(x)MU1(temp=SST[x], 
                                               par_=env[i,]$PAR,
                                               NO3 =env[i,]$DIN)$QN) 

    if(i==1) {
       plot(SST, QN_FIX[,i], type = 'l',  lty=i,col=3,
            ylim = c(params$Q0N,4.1*params$Q0N),
            xlab = 'Temperature (ºC)',
            ylab = expression(paste('N:C (mol '*mol^-1*')')))
    }else{
       points(SST, QN_FIX[,i], type = 'l',lty=i,col=3)
    }
    points(SST, QN_EFT[,i], type='l', lty=i,col=2)
}
mtext('i',adj=0)

dev.off()

#Plot 1D results:
source('~/Working/FlexEFT1D/Rscripts/surface_OBS_model.R')
source('~/Working/FlexEFT1D/Rscripts/MLD_PHYZOO.R')
source('~/Working/FlexEFT1D/Rscripts/CChl_ML.R')
source('~/Working/FlexEFT1D/Rscripts/vertical_OBS_model_MS.R')
source('~/Working/FlexEFT1D/Rscripts/vertical_muQNtheta.R')
source('~/Working/FlexEFT1D/Rscripts/Vertical_QN.R')

pdf('Growth.pdf')
op <- par(font.lab = 1,
             family ="serif",
             mar    = c(2,2,1.5,1.5),
             mgp    = c(2.3,1,0),
             oma    = c(4,4,0,0),
             mfcol  = c(2,2)   ) 
 plot_1D('muN1',model,Stn)
 plot_1D('muN2',model,Stn)
 plot_1D('PHY1',model,Stn)
 plot_1D('PHY2',model,Stn)
dev.off()


model <- 'EFT2sp'
Stn   <- 'S1'
setwd('~/Working/FlexEFT1D/DRAM_0.9/EFT2sp/BOTH/')  
pdf('Growth.pdf')
op <- par(font.lab = 1,
             family ="serif",
             mar    = c(2,2,1.5,1.5),
             mgp    = c(2.3,1,0),
             oma    = c(4,4,0,0),
             mfcol  = c(2,2)   ) 
 plot_1D('muN1',model,Stn)
 plot_1D('muN2',model,Stn)
 plot_1D('PHY1',model,Stn)
 plot_1D('PHY2',model,Stn)
dev.off()


#Test the models with data:

source('~/Working/FlexEFT1D/Rscripts/PhyQN.R')
