#Flexible:
MU_P  <- function(NO3,par_,temp)mu_EFT(temp=temp,NO3=NO3,par_ = par_,
                  aI=params$aI0,Q0N=params$Q0N,mu0=params$mu0hat,A0N=params$A0N)


MONOD <- function(NO3,par_,temp){
         mu_fix(Temp=temp,NO3=NO3,par_=par_,mu0=rmax,KN0=KN,
                Q0N=params$Q0N,aI0C=aI0, thetamax=themax)
}

setwd('/Users/apple/Working/Draft/Opt_model/')
#Get data file:
file = 'Phyto_N_C.csv'
dat  = read.csv(file)

#Calculate integrated light (convert to W m-2):
dat$PAR = dat$Light*1E-6*3600*24*dat$LD/.4
dat$lnPAR = log(dat$PAR)
IDs     = unique(dat$ID)
Nsp     = length(IDs)

#Run a LME model: 
library(lme4)
QN_lme  = lmer(QN ~ lnPAR + (lnPAR|ID), data=dat)
df      = coef(summary(QN_lme))
df2     = coef(QN_lme)$ID

#Predictions of the Markus model:
newpar   = seq(min(dat$PAR), max(dat$PAR), 0.01)
newQN    = sapply(1:length(newpar),function(x)MU_P(par_=newpar[x],NO3=20,temp=20)$QN) 
newtheta = sapply(1:length(newpar),function(x)MU_P(par_=newpar[x],NO3=20,temp=20)$Theta) 

pdf('Fig3_QN_theta_PAR.pdf', width=7, height=9, paper='a4')
op <- par(font.lab = 1,
            family ="serif",
            mgp    = c(2.3,1,0),
            oma    = c(4,4,0,0))
zones=matrix(1:4, ncol=2, byrow=TRUE)
layout(zones, widths=c(4/5,1/5),heights=c(1,1))
par(mar=c(2,4,2,0))
QNLim  = c(.1, .56)
plot(dat$lnPAR, dat$QN, type = 'n', ylim=QNLim,
     xlab = expression(paste('Log PAR (W '*m^-2*')')),
     ylab = 'N:C molar ratio')
mtext('a)', cex=1.2,adj=0, line = .5)
COLORS = c(1:6,1:3)
txt    = character(Nsp)
for (i in 1:Nsp){
    cff = dat[dat$ID == IDs[i],]
    points(cff$lnPAR, cff$QN, pch = 15 + i, col = COLORS[i])
    txt[i] = paste(cff$Ref[1],cff$Species[1])
    newx   = seq(min(cff$lnPAR), max(cff$lnPAR), 0.01)
    newy   = df2[i,1]+df2[i,2]*newx
    #Plot individual lines of the LME model:
    lines(newx,newy, lwd=.4, col=COLORS[i],lty=3)
}
legend('topright',txt,col=COLORS, pch=16:(15+Nsp) )

#Plot the mean line of LME model:
abline(df[1,1],df[2,1], lwd=2)

#Plot the Pahlow model prediction:
lines(log(newpar),newQN, lwd=2, lty=4, col=2)

#Plot redfield ratio:
abline(h=16/106,  lty=2, col='blue')

#Plot histogram of QN:
yhist = hist(c(dat$QN,max(QNLim)), plot=FALSE, breaks=25)
yhist$counts[length(yhist$counts)] = 0
par(mar=c(2,0,2,.1))
barplot(yhist$counts, 
        axes=FALSE, space=0, horiz=TRUE)

par(mar=c(4,4,1.5,0))
plot(dat$lnPAR, dat$theta, type = 'n',
     xlab = expression(paste('Log PAR (W '*m^-2*')')),
     ylab = 'Chl:C (g:mol)')
mtext('b)', cex=1.2,adj=0, line = .5)

#LME of theta (2nd order)
theta_lme  = lmer(theta ~ lnPAR + I(lnPAR^2) + (0+lnPAR|ID)+(0+I(lnPAR^2)|ID)+(1|ID),
                  data=dat)
df         = coef(summary(theta_lme))
df2        = coef(theta_lme)$ID

for (i in 1:Nsp){
    cff = dat[dat$ID == IDs[i],]
    points(cff$lnPAR, cff$theta, pch = 15 + i, col = COLORS[i])
    newx   = seq(min(cff$lnPAR), max(cff$lnPAR), 0.01)
    newy   = df2[i,1]+df2[i,2]*newx + df2[i,3]*newx**2
    #Plot individual lines of the LME model:
    lines(newx,newy, lwd=.4, col=COLORS[i],lty=3)
}

#Plot the mean line of LME model:
lme.theta = df[1,1] + df[2,1]*log(newpar) + df[3,1]*log(newpar)**2
lines(log(newpar),lme.theta, lwd=2)

#Plot the Pahlow model prediction:
lines(log(newpar),newtheta, lwd=2, lty=4, col=2)

#Plot classic ratio of 50:1:
abline(h=12/50,  lty=2, col='blue')
legend('topright',c('PAHLOW','LME'),col=c(2,1),lty=c(4,1),lwd=2)
#Plot histogram of QN:
yhist = hist(dat$theta, plot=FALSE, breaks=25)
par(mar=c(4,0,1.5,.1))
barplot(yhist$counts, axes=FALSE, space=0, horiz=TRUE)

dev.off()


