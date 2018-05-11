#Plot LogLikehood of two models at two stations in one plot:
pwd1  <- getwd()
model <- 'NPZclosure'
stn   <- 'S1'
setwd(paste0('~/Working/FlexEFT1D/DRAM/',model, '/',stn, '/'))
Nstn       <- 1
burnin     <- 100
NDTYPE     <- 3  #The number of obs. types
np         <- 5  #The number of CPUs for paralell computing
EnsLen     <- 100  #The number of ensembles
enssig     <- read.table('enssig',header=T)
enspar     <- read.table('enspar',header=T)

#Convert params to original unit:
enspar[,2:ncol(enspar)] <- exp(enspar[,2:ncol(enspar)])


#Get bestpar:
best    = which.max(enspar$LogL)
bestpar = enspar[best,]

#Best SSqE:
bestsig = enssig[which.max(enssig$LogL),]

NPAR    = ncol(enspar)-1
runno   = nrow(enssig)/np   #Number of runs for each processor
intv    = runno/EnsLen      #Number of runs each ensemble

#To break out into several processes:
loglike = matrix(NA,nr  = runno, nc = np)
ssqe    = array(NA, dim = c(runno, NDTYPE*Nstn, np))
params  = array(NA, dim = c(runno, NPAR       , np))
for (j in 1:np){
  w    = c()
  for (i in 1:EnsLen){
     w = c(w,((j-1)*intv+(i-1)*np*intv+1):((j-1)*intv+(i-1)*np*intv+intv))
  }
  loglike[,j]=enssig[w,1]
  ssqe[,,j]  =as.matrix(enssig[w,(ncol(enssig)-NDTYPE*Nstn+1):ncol(enssig)])
  params[,,j]=as.matrix(enspar[w,                           2:ncol(enspar)])
}
fx = burnin:runno
pdf('LogLike.pdf',
      width=3*2,height=4*2,paper='a4')
op <- par(font.lab = 1,
           family ="serif",
           mar    = c(2,4,2,.3),
           mgp    = c(2,1,0),
           mfrow  = c(2,2),
           oma    = c(3,0,0,0))
plot(fx/1E+3,loglike[fx,1],
      ylim=range(loglike[fx,]),
      xlab='',
      ylab='Log-likelihood',
      lwd =0.4,type='n')
for (i in 1:np){
   points((fx)/1E+3,loglike[fx,i],type='l')
}

ii = 1
mtext(letters[ii],adj=0,cex=.8)

varnames <- c('DIN','Chl','NPP')

for (i in 1:(NDTYPE*Nstn)){

    plot(fx/1E+3,ssqe[fx,i,1],
         ylim=range(ssqe[fx,i,]),
         xlab='',
         ylab=paste('SSqE_',varnames[i],sep=''),
         lwd =0.4,type = 'n')
    for (iii in 1:np){
       points((fx)/1E+3,ssqe[fx,i,iii],type='l')
    }

    ii=ii+1
    mtext(letters[ii],adj=0,cex=1)
}
xlab  = expression(paste('Run number ('*10^3*')',sep=''))
mtext(xlab,side=1,outer=TRUE,cex=.8)
dev.off()

paramlab     <- character(NPAR)
paramlab[1]  <- expression(paste(italic(µ)['0']))
paramlab[2]  <- expression(paste(italic(I)[opt]))
paramlab[3]  <- expression(paste(italic(alpha)))
paramlab[4]  <- expression(paste(italic(K)['N']))
paramlab[5]  <- expression(paste(italic(D)['P']))
paramlab[6]  <- expression(paste(italic(W)))
paramlab[7]  <- expression(paste(italic(beta)))

if (model == 'NPZclosure'){
   paramlab[8]  <- expression(paste(italic(g)['max']))
   paramlab[9]  <- expression(paste(italic(m)['z']))
}

#Plot time evolution of different parameters:
pdf('Params.pdf',
      width=4*2,height=3*2,paper='a4')
op <- par(font.lab = 1,
           family ="serif",
           mar    = c(2,4,1,.3),
           mgp    = c(2,1,0),
           oma    = c(3,0,0,0),
           mfrow  = c(ceiling(NPAR/2),2), cex.lab=1.2, cex.axis=1.2)

for (i in 1:NPAR){
    plot(fx/1000, params[fx,i,1],
         ylim=range(params[fx,i,]),
         xlab='',
         ylab=paramlab[i],
         lwd =0.4,type='n')
    mtext(letters[i],adj=0,cex=.8)
    for (k in 1:np){
       points(fx/1000, params[fx,i,k], type='l')
    }
}
xlab  = expression(paste('Run number ('*10^3*')',sep=''))
mtext(xlab,side=1,outer=TRUE,line=1,cex=.8)
mtext('Time evolution of fitted model parameters.',side=1,outer=T, line=2,adj=0)
dev.off()

##Plot pairs of parameter values together with LogLike and their correlations with Loglike:

enspar1   <- enspar[enspar$LogL > 100, ]
#Convert Loglike to vectors
Loglike2    <- quantile(enspar1$LogL, probs=c(0.1,0.5,0.9))
enspar1$Col <- cut(enspar1$LogL, breaks = c(-Inf,Loglike2,+Inf), 
                 labels = jet.col(n=(length(Loglike2)+1)),
                 right  = FALSE)

Rcol <- unique(enspar1$Col[order(enspar1$LogL)])
pdf('params_loglike.pdf', width = 9, height = 6, paper = 'a4')
op <- par(font.lab = 1,
            family ="serif", cex.axis=1, cex.lab=1.2,
            mar    = c(4,4,1.5,0.5),
            mgp    = c(2.3,1,0),
            mfcol  = c(2,2),
            oma    = c(4,4,1,0)) 

for (i in 2:(NPAR-1)){
    for (j in (i+1):(NPAR+1)){
       plot(enspar1[,i],enspar1[,j], 
            xlab=paramlab[i-1],ylab=paramlab[j-1], pch=16, cex=.5,col=enspar1$Col)

       #Plot legend:
       if (i == 2 && j == 3){
          legend('topleft',pch=16, cex=.5,
              legend=paste0('Loglike < ', c(Loglike2,max(enspar1[,1]))),
              border=Rcol,fill=Rcol,col=Rcol)
       }
    }
}
dev.off()
