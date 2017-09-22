#Plot LogLikehood of two models at two stations in one plot:
pwd1 <- getwd()
setwd('~/Working/FlexEFT1D/DRAM/NPZDcont/BOTH_TD/')  
stns    = c('K2','S1')
Nstn    = 2
burnin  = 100
NDTYPE  = 9  #The number of obs. types
np      = 5  #The number of CPUs for paralell computing
EnsLen  = 100  #The number of ensembles
enssig  = read.table('enssig',header=T)
enspar  = read.table('enspar',header=T)

#Get bestpar:
best    = which.max(enspar$LogL)
bestpar = enspar[best,]

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
pdf('Fig4.LogLike_size.pdf',
      width=3*2,height=4*2,paper='a4')
op <- par(font.lab = 1,
           family ="serif",
           mar    = c(2,4,2,.3),
           mgp    = c(2,1,0),
           oma    = c(3,0,0,0))
nf <- layout(matrix(c(rep(1,3),2:19), 7, 3, byrow = T), respect = F)
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

varnames <- c('DIN','Chl','NPP','PON','Fer','P10','P03','P01','P_1')
varnames <- rep(varnames,2)
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
    if (ii%%9 == 2){
       if (ii > 6) {
          j <- 2
       } else{
          j <- 1
       }
       mtext(paste(letters[ii],')',stns[j]),adj=.6)
    }else{
       mtext(letters[ii],adj=0,cex=1)
    }
}
xlab  = expression(paste('Run number ('*10^3*')',sep=''))
mtext(xlab,side=1,outer=TRUE,cex=.8)
mtext('Fig. 4',side=1,outer=T, line=2,adj=0)
dev.off()

paramlab = character(NPAR)
source('~/Working/FlexEFT1D/Rscripts/param_lab.R')
#Plot time evolution of different parameters:
pdf('Fig5.Params.pdf',
      width=3*2,height=3*2,paper='a4')
op <- par(font.lab = 1,
           family ="serif",
           mar    = c(2,4,1,.3),
           mgp    = c(2,1,0),
           oma    = c(3,0,0,0),
           mfrow  = c(3,3), cex.lab=1.2, cex.axis=1.2)

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
mtext('Fig. 5 Time evolution of fitted model parameters.',side=1,outer=T, line=2,adj=0)
dev.off()



##Plot pairs of parameter values together with LogLike:
#burnin = 1E3
#fx = burnin:runno
##Convert Loglike to vectors
#Loglike1 = as.vector(loglike[fx,]) 
#
##Convert params to matrix:
#PARS = matrix(NA, nr = length(fx)*np, nc = NPAR)
#
#for (i in 1:NPAR){
#  PARS[,i] = as.vector(params[fx, i, ])
#}
#
#PARS = as.data.frame(PARS)
#colnames(PARS) = names(enspar)[2:(NPAR+1)]
#
#Loglike2 = quantile(Loglike1,probs=c(0.1,0.25,0.5,0.75,0.9))
#PARS$Colour <- cut(Loglike1, breaks = c(-Inf,Loglike2,+Inf), 
#                 labels = jet.col(n=6),
#                 right  = FALSE)
#
#pdf('params_loglike.pdf', width = 9, height = 6, paper = 'a4')
#op <- par(font.lab = 1,
#            family ="serif", cex.axis=1.2, cex.lab=1.2,
#            mar    = c(2,2,1.5,3.5),
#            mgp    = c(2.3,1,0),
#            mfcol  = c(4,4),
#            oma    = c(4,4,1,0)) 
#
#for (i in 1:(NPAR-1)){
#    for (j in (i+1):NPAR){
#        plot(PARS[,i],PARS[,j], xlab=paramlab[i],ylab=paramlab[j], pch=16, cex=.2,col=PARS$Colour)
#       # plot(PARS[,i],PARS[,j], pch=16, cex=.2,col=jet.col(8)[8])
#    #Plot legend:
#    if (i == 1 && j == 2){
#    legend('bottomleft',pch=16,legend=c(Loglike2,max(Loglike1)),
#           border=jet.col(6),fill=jet.col(6),col=jet.col(n=6))
#    }
#    }
#}
#
#
#
#dev.off()
#setwd(pwd1)
