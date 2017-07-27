#Plot LogLikehood of two models at two stations in one plot:
pwd1 <- getwd()
setwd('~/Working/FlexEFT1D/DRAM_0.9/NPZDcont/BOTH/')  
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
burnin  = 10
NDTYPE  = 8
#filedir = paste('~/working/FlexEFT1D/DRAM_0.2/',stn,'/',model,sep='')
#setwd(filedir)
enssig  <- read.table('enssig',header=T)
enssig  <- enssig[burnin:nrow(enssig),]
runno   <- enssig[,1]
loglike <- enssig[,2]
sigma   <- enssig[,(2+1):(2+NDTYPE*Nstn)]
ssqe    <- enssig[,(ncol(enssig)-NDTYPE*Nstn+1):ncol(enssig)]

plot(runno/1E+3,loglike,
      xlab='',
      ylab='Log-likelihood',
      lwd =0.4,type='l')

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
