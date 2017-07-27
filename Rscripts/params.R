pwd1 <- getwd()
setwd('~/Working/FlexEFT1D/DRAM_0.9/NPZDcont/BOTH/')  
#Open enspar 
enspar = read.table('enspar1',header=T)
Z      = enspar[, 3:ncol(enspar)]
nms    = names(enspar)[3:ncol(enspar)]
burnin = 2000
enspar = enspar[burnin:nrow(enspar),]
runno  = enspar[,1]
paramlab = character(length(nms))
source('~/Working/FlexEFT1D/Rscripts/param_lab.R')
#Plot time evolution of different parameters:
pdf('Fig4.Params.pdf',
      width=3*2,height=5*2,paper='a4')
op <- par(font.lab = 1,
           family ="serif",
           mar    = c(2,4,1,.3),
           mgp    = c(2,1,0),
           oma    = c(3,0,0,0),
           mfrow  = c(5,3), cex.lab=1.2, cex.axis=1.2)

for (i in 1:length(nms)){
    plot(runno/1000, enspar[, i+2],
            xlab='',
            ylab=paramlab[i],
            lwd =0.4,type='l')
    mtext(letters[i],adj=0,cex=.8)
}
xlab  = expression(paste('Run number ('*10^3*')',sep=''))
mtext(xlab,side=1,outer=TRUE,line=1,cex=.8)
mtext('Fig. 4 Time evolution of fitted model parameters.',side=1,outer=T, line=2,adj=0)
dev.off()
setwd(pwd1)
