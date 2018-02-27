setwd('~/Working/FlexEFT1D/VerticalNPP')  
#Process PON and N:C data at HOT:
PON <- '~/Working/FlexEFT1D/HOT/HOT_POC.csv'
PON <- read.csv(PON)
PON$PON = PON$pn
PON$N2C = PON$pn/PON$pc
PON <- PON[PON$date > 0 & PON$pn > 0, ]

gfile = function(VAR){
   #Remove invalid date:
   PON <- PON[PON$date > 0 & PON[,VAR] > 0, ]
   #Calculate DOY:
   DOY <- function(date){
       date <- as.character(date)
       N    <- nchar(date)
       yy   <- as.integer(substr(date, N-1,N))
       if (yy > 20) {
         Year <- yy + 1900
       } else{
         Year <- yy + 2000
       }
       day   <- as.integer(substr(date, N-3,N-2))
       mo    <- as.integer(substr(date, 1,  N-4))
       Date  <- paste0(day,'/',mo,'/',Year)
       Date  <- format(Date, format = "%d/%m/%y")
       return(as.numeric(strftime(Date, format = "%j")))
   }
   PON$DOY <- sapply(1:nrow(PON),function(i)DOY(PON$date[i]))
   cff     <- data.frame(DOY=PON$DOY, Depth=PON$press, Data=PON[,VAR])
   fname   <- paste0('~/Working/FlexEFT1D/HOT/HOT_',VAR,'.dat')
   write.table(cff, file = fname, row.names = F)
}

gfile('N2C')

#Process CTD data at HOT:
CTD  <- '~/Working/FlexEFT1D/HOT/HOT_CTD.csv'
CTD  <- read.csv(CTD)
#1988/10/1 the first day
CTD$Date <- as.Date(CTD$julian, origin = "1988-10-01")
CTD$Date <- as.character(format(CTD$Date,'%m/%d/%Y'))
CTD$Depth<- CTD$press
CTD$Sigma_Theta <- CTD$sigma
CTD$Station     <- 'HOT'
fname           <- '~/Working/FlexEFT1D/Obs_data/HOT_obs.csv'
write.csv(CTD, fname,row.names=F)


get.sum.vdata <- function(Stn,Var){
   #Read observational data
   dat <- paste0('~/Working/FlexEFT1D/',Stn,'/',Stn,'_',Var,'.dat')
   dat <- read.table(dat,header=T)
   
   #Select only summer data (July ~ Sep):
   DOYs    <- seq(0,360,length.out=5)
   
   dat <- dat[dat$DOY >= DOYs[3] & dat$DOY <= DOYs[4],]
   dat$x    <- dat[,3]
   l1       <- loess(x ~ -Depth, data=dat, span=0.5)
   newDepth <- 1:200
   l2       <- predict(l1, newdata = data.frame(Depth=newDepth)) 
   return(list(Depth=newDepth,dat=as.numeric(l2)))
}
#Get TIN data:
TIN_sum = get.sum.vdata('HOT','TIN')

#Get DIP data:
DIP_sum  = get.sum.vdata('HOT','DIP')

#Plot vertical N:P ratio:
pdf('Summer_NtoP.pdf',width=4,height=4,paper='a4')
p <- par(font.lab  = 1,
             family ="serif",
             mar    = c(4,4,1.5,0.5),
             mgp    = c(2.3,1,0),
             oma    = c(4,4,0,0),
             mfcol  = c(1,1)   )
plot(TIN_sum$dat/DIP_sum$dat,-TIN_sum$Depth,type = 'l',
      xlab = 'N:P', ylab = 'Depth (m)') 
dev.off()


plot.obs <- function(Stn, Var){
   #Read observational data
   dat <- paste0('~/Working/FlexEFT1D/',Stn,'/',Stn,'_',Var,'.dat')
   dat <- read.table(dat,header=T)
   
   #Select only summer data (July ~ Sep):
   DOYs    <- seq(0,360,length.out=5)
   
   dat <- dat[dat$DOY >= DOYs[3] & dat$DOY <= DOYs[4],]
   X.range  <- quantile(dat[,3], prob=c(0.01,0.99), na.rm = T)

   if (Var == 'TIN' || Var == 'DOC' || Var == 'DIP'){
      Varname  <- bquote(.(Var)~' (mmol '*m^-3*')')
      if (Var == 'DIP'){
        X.range  <- c(0, .4)
      }else{
        X.range  <- c(0, 5)
      }
   }else if (Var == 'CHL'){
     Varname  <-  bquote(.(Var)~' (mg '*m^-3*')')
      X.range <- c(0, .5)
   }else if (Var == 'NPP'){
     Varname  <- bquote(.(Var)~' (mg C '*m^-3*' '*d^-1*')')
      X.range <- c(0, 10)
   }else if (Var == 'N2C'){
     Varname  <- 'N:C (mol:mol)'
   }else if (Var %in% c('PRO','SYN','EUK')){
     Varname  <- bquote(.(Var)~' ('*10^5*' cells '*mL^-3*')')
   }



   if (Var == 'N2C') X.range = c(0.1,0.2)
   plot(dat[,3],-dat$Depth, xlim=X.range, ylim=c(-200,0),
            xlab=Varname,ylab='Depth (m)',pch=16,cex=.5,cex.lab=1.2)
   if (nrow(dat) > 20) {
      dat$x    <- dat[,3]
      l1       <- loess(x ~ -Depth, data=dat, span=0.5)
      newDepth <- 1:200
      l2       <- predict(l1, newdata = data.frame(Depth=newDepth), 
                                   se = TRUE)
      lines(l2$fit, -newDepth, lwd=2)
   }
   if (Stn == 'HOT') Stn = 'ALOHA'
   if (Var == 'N2C') abline(v=16/106)
   mtext(paste(Stn,Var),adj=0)
}

fname <- paste0('Vertical_summer_ChlNPP3.pdf')
pdf(fname, width=9,height=9,paper='a4')

op <- par(font.lab  = 1,
             family ="serif",
             mar    = c(4,4,1.5,0.5),
             mgp    = c(2.3,1,0),
             oma    = c(4,4,0,0),
             mfcol  = c(3,3)   ) 

plot.obs('HOT','TIN')
plot.obs('HOT','CHL')
plot.obs('HOT','NPP')

plot.obs('BATS','TIN')
plot.obs('BATS','CHL')
plot.obs('BATS','NPP')

plot.obs('S1','TIN')
plot.obs('S1','CHL')
plot.obs('S1','NPP')
dev.off()

fname <- paste('Vertical_summer_N2C.pdf',sep='')
pdf(fname, width=7,height=4,paper='a4')
op <- par(font.lab  = 1,
             family ="serif",
             mar    = c(4,4,1.5,0.5),
             mgp    = c(2.3,1,0),
             oma    = c(4,4,0,0),
             mfcol  = c(1,2)   ) 
plot.obs('HOT','N2C')
plot.obs('S1', 'N2C')
dev.off()


fname <- paste('Vertical_summer_HOT_DOC.pdf',sep='')
pdf(fname, width=6,height=9,paper='a4')

op <- par(font.lab  = 1,
             family ="serif",
             mar    = c(4,4,1.5,0.5),
             mgp    = c(2.3,1,0),
             oma    = c(4,4,0,0),
             mfcol  = c(2,2)   ) 
plot.obs('HOT','DOC')
plot.obs('HOT','PRO')
plot.obs('HOT','SYN')
plot.obs('HOT','EUK')
dev.off()

newDepth <- 1:200
#Get Theta from Fennel and Boss (2003):
Fennel.theta <- 'theta_Fennel03.csv'
F.theta      <- read.csv(Fennel.theta)
l1 <- loess(theta ~ Depth, data=F.theta, span=0.5)
l2 <- predict(l1, newdata = data.frame(Depth=newDepth)) 
l2 <- as.numeric(l2)

l2[is.na(l2) & newDepth < min(F.theta$Depth)] <- l2[!is.na(l2)][1]
l2[is.na(l2) & newDepth > max(F.theta$Depth)] <- l2[!is.na(l2)][length(l2[!is.na(l2)])]

#Get Chl at HOT
HOT.CHL <- get.sum.vdata(Stn='HOT',Var='CHL')$dat
HOT.CHL[is.na(HOT.CHL)] <- HOT.CHL[!is.na(HOT.CHL)][1]

#Calculate Phytoplankton Nitrogen at HOT:
HOT.Cphy <- HOT.CHL/l2 *106/16*12   #Unit: µg C/L

fname <- paste('Vertical_summer_HOT_mu.pdf',sep='')
pdf(fname, width=5,height=3,paper='a4')

op <- par(font.lab  = 1, las=1, lwd=1.5,
             family ="serif",
             mar    = c(4,4,1.5,1.5),
             mgp    = c(2.3,1,0),
             oma    = c(4,4,0,0),
             mfcol  = c(1,2)   ) 
plot(HOT.Cphy,-newDepth, 
     xlab = 'Phyto carbon (µg/L)',
     ylab = 'Depth (m)',
     type = 'l')

plot(HOT.muP,-newDepth, 
     xlab =  expression(paste('Growth rate ('*d^-1*')')),
     ylab = 'Depth (m)',
     type = 'l')

dev.off()

#Get NPP at HOT
HOT.NPP <- get.sum.vdata(Stn='HOT',Var='NPP')$dat

#elements that are not NA
k <- which(!is.na(HOT.NPP))

HOT.NPP[is.na(HOT.NPP) & newDepth < k[1] ] <- HOT.NPP[!is.na(HOT.NPP)][1]

#Calculate phytoplankton specific growth rate:
HOT.muP <- HOT.NPP/HOT.Cphy

#Calculate Theta from S1:

source('~/Working/FlexEFT1D/Rscripts/vertical_PHYmu.R')

Stn   <- 'S1'
#Obtain theta and QN:

get.vprof <- function(Stn,VAR,model='EFTsimple'){
   DIR   <- paste0('~/Working/FlexEFT1D/DRAM_0.8/',model,'/',Stn,'/')
   Dat   <- getData(DIR, Stn,VAR)
   days  <- Dat$days
   depth <- Dat$depth
   Dat   <- Dat$data
   d_per_y <- 360
   #Get the data of the final year
   w       <- (nrow(Dat)-d_per_y+1):nrow(Dat)
   Dat     <- Dat[w,]
    #Read model results:
   mdoy    <- 1:nrow(Chl)
   #Average into 4 seasons
   DOYs    <- seq(0,360,length.out=5)
   
   #Select only summer data
   k     <- which(mdoy > DOYs[3] & mdoy <= DOYs[3+1])
   Dat_  <- Dat[k,]
   Dat_  <- apply(Dat_,2,mean)
   dat   <- data.frame(Depth=-depth,x=Dat_)
   l1       <- loess(x ~ -Depth, data=dat, span=0.5)
   l2       <- predict(l1, newdata = data.frame(Depth=newDepth)) 
   l2       <- as.numeric(l2)
   l2[is.na(l2)] <- l2[!is.na(l2)][1]
   return(list(Depth=newDepth,dat=as.numeric(l2)))
}

#Get theta at S1 from model
S1.theta <- get.vprof(Stn='S1',VAR='The1')$dat

#Get Chl at S1
S1.CHL <- get.sum.vdata(Stn='S1',Var='CHL')$dat

#Get Carbon
S1.Carbon <- S1.CHL / S1.theta*12

#Get NPP at S1
S1.NPP <- get.sum.vdata(Stn='S1',Var='NPP')$dat
S1.NPP[S1.NPP < 0] <- NA
S1.mu <- S1.NPP/S1.Carbon

plot(S1.mu, -newDepth, type = 'l')
plot(S1.Carbon, -newDepth, type = 'l')

#Check model results:
#Download DRAM results:
system('scp -r bzchen@chula.yes.jamstec.go.jp:/data16/bzchen/FlexEFT1D/DRAM_0.9/NPZDFix/ ~/Working/FlexEFT1D/DRAM_0.9/')
#Read enspar:
setwd('~/Working/FlexEFT1D/DRAM_0.9/NPZDFix')
enspar.HOT <- read.table('HOT/enspar',header=T)
enspar.S1  <- read.table('S1/enspar',header=T)

source('~/Working/FlexEFT1D/DRAM_0.9/Vertical_NPZD.R')

#Plot vertical patterns of µ, g, theta, QN, 


#Plot phytoplankton traits:
#Test the relationship between theta and 
