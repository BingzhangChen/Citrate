setwd('~/Working/FlexEFT1D/BATS')
depth_total    <- 500 #Total depth (m) of the station

#Read NPP data:
file <-  'bats_NPP.csv'
dat  <-  read.csv(file)
ND   <- which(colnames(dat) == 'Temp')
for (k in ND:ncol(dat)){
    x        <- dat[,k] 
    x        <- as.numeric(as.character(x))
    x[x < 0] <- NA
    dat[,k]  <- x
}

#Convert cast time to DOY:
dat$date  <- as.character(dat$Date)
dat$date  <- as.Date(dat$date,'%Y%m%d')
dat$DOY   <- as.numeric(strftime(dat$date, format = "%j")) 
dat0      <- dat  #Save original data

#Retrieve only upper 500 meters
dat       <- dat[dat$Depth <= depth_total, ]

dat1      <- dat[,c('DOY','Depth','PP')]
colnames(dat1)[3] <- 'NPP' 
dat1 <- na.omit(dat1)

#write out csv file:
file <- 'BATS_NPP.dat'
write.table(dat1,file,row.names=F)

#Read Chl data:
dat  <-  read.csv('bats_Chl.csv')
#Convert cast time to DOY:
dat$date  <- as.character(dat$yyyymmdd)
dat$date  <- as.Date(dat$date,'%Y%m%d')
dat$DOY   <- as.numeric(strftime(dat$date, format = "%j")) 
Chl0      <- dat  #Save original data

#Retrieve only upper 500 meters
dat       <- dat[dat$Depth <= depth_total, ]

dat1      <- dat[,c('DOY','Depth','p16')]
colnames(dat1)[3] <- 'Chl' 
dat1[dat1[,3] < 0, 3] <- NA
dat1 <- na.omit(dat1)

#write out csv file:
file <- 'BATS_CHL.dat'
write.table(dat1,file,row.names=F)

#For DIN, DIP, PON, POP
dat <- read.csv('bats_bot.csv')

#Convert cast time to DOY:
dat$date  <- as.character(dat$yyyymmdd)
dat$date  <- as.Date(dat$date,'%Y%m%d')
dat$DOY   <- as.numeric(strftime(dat$date, format = "%j")) 
Bottle    <- dat  #Save original data

#Retrieve only upper 500 meters
dat0      <- dat[dat$Depth <= depth_total, ]

for (var in c('TIN','PON','DIP','POP')){
   dat <- dat0 
   if (var == 'TIN'){
      dat[,var] <- dat$NO3
   } else if (var == 'PON'){
      dat[,var] <- dat$PON/14  #Unit: µg/kg!!!
   } else if (var == 'POP'){
      dat[,var] <- dat$POP
   } else if (var == 'DIP'){
      dat[,var] <- dat$PO4 #Normal DIP
      dat$SRP[dat$SRP < 0] <- NA
      #Use low concentration DIP if available
      for (i in 1:nrow(dat)){
          if (!is.na(dat$SRP[i])) dat[i,var] <- dat$SRP[i]/1E3
      }
   } else{
     stop('Variable name incorrect!')
   }
   dat  <- dat[!is.na(dat[,var]),]
   dat  <- dat[!is.na(dat$DOY),  ]
   dat  <- dat[dat[,var] > 0,]
   dat1 <- dat[,c('DOY','Depth',var)]
   for (j in 1:ncol(dat1)){
       dat1[,j] <- as.numeric(dat1[,j])
   }

   #write out csv file:
   file <- paste0('BATS','_',var,'.dat')
   write.table(dat1,file,row.names=F)
}

#For PON and POP, produce interpolated bottom values for boundary condition:
for (VAR %in% c('PON','POP')){
   file <- paste0('BATS','_',VAR,'.dat')
    dat <- read.table(file, header = T)
   DOYu <- unique(dat$DOY)
   Bdat <- numeric(length(DOYu))
   Bdat <- NA
   data <- data.frame(DOY = DOYu, dat = Bdat)
   #First, obtain the data at bottom depth for each profile:
   for (i in 1:length(DOYu)){
       cff <- dat[dat$DOY == DOYu[i],]
       if (max(cff$Depth) > 400) {
          w       <- which.max(cff$Depth)
          data[i,2] <- cff[w,3]
       }
   }
   
   LOESS <- loess(dat ~ DOY, data)
   # Interpolate for 12 months:
   newM  <- seq(0.5,11.5,length.out=12)
   newM  <- newM*30
   newY  <- predict(LOESS, DOY = newM)
}

#For HOT, PON = 0.04, POP = .002

