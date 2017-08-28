#Plot LogLikehood of two models at two stations in one plot:
pwd1 <- getwd()
setwd('~/Working/FlexEFT1D/DRAM/NPZDcont/BOTH_TD/')  
#pdffile <- paste('LogLike.pdf',sep='')
#pdf(pdffile, width=5,height=5,paper='a4')
pdf('Fig3.LogLike_size.pdf',
      width=4*2,height=5*2,paper='a4')
op <- par(font.lab = 1,
           family ="serif",
           mar    = c(2,4,2,.3),
           mgp    = c(2,1,0),
           oma    = c(3,0,0,0))
nf <- layout(matrix(c(rep(1,4),2:17), 5, 4, byrow = T), respect = TRUE)
ii <- 1
stns <- c('K2','S1')
#for (stn in c('S1','K2')){
#  for (model in c('EFTdiscrete','EFTcont')){ 
Nstn    = 2
burnin  = 2
NDTYPE  = 9  #The number of obs. types
np      = 4  #The number of CPUs for paralell computing
EnsLen  = 2  #The number of ensembles
#filedir = paste('~/working/FlexEFT1D/DRAM_0.2/',stn,'/',model,sep='')
#setwd(filedir)
enssig  = read.table('enssig',header=T)
runno   = nrow(enssig)/np   #Number of runs for each processor
intv    = runno/EnsLen      #Number of runs each ensemble

#To break out into several processes:
loglike = matrix(NA,nr  = runno, nc = np)
ssqe    = array(NA, dim = c(runno, NDTYPE*Nstn, np))
for (j in 1:np){
  w    = c()
  for (i in 1:EnsLen){
     w = c(w,((j-1)*intv+(i-1)*np*intv+1):((j-1)*intv+(i-1)*np*intv+intv))
  }
  loglike[,j]=enssig[w,1]
  ssqe[,,j]  =as.matrix(enssig[w,(ncol(enssig)-NDTYPE*Nstn+1):ncol(enssig)])
}

loglike1 = apply(loglike[burnin:runno,],1,mean)
ssqe    =    ssqe[burnin:runno,,]

plot((burnin:runno)/1E+3,loglike[burnin:runno,1],
      ylim=range(loglike[burnin:runno,]),
      xlab='',
      ylab='Log-likelihood',
      lwd =0.4,type='n')
for (i in 1:np){
   points((burnin:runno)/1E+3,loglike[burnin:runno,i],type='l')
}

      #ii=ii+1
      mtext(letters[ii],adj=0,cex=.8)

      varnames <- c('TIN','Chl','NPP','PON','P10','P03','P01','P_1')
      varnames <- rep(varnames,2)
      for (i in 1:ncol(ssqe)){

          var1 <- paste0()
          plot(runno/1E+3,ssqe[,i],
               xlab='',
               ylab=paste('SSqE_',varnames[i],sep=''),
               lwd =0.4,type = 'l')
          ii=ii+1
          if (ii%%8 == 2){
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
mtext('Fig. 3',side=1,outer=T, line=2,adj=0)
dev.off()


setwd(pwd1)
