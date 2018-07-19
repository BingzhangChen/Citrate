beta = c(0.01, 0.1,0.5, 1, 2)
A    = 2
igm  = 1
imz  = 2
iKN  = 3
NPAR = iKN
MAX  = numeric(NPAR)
MIN  = MAX
MIN[igm] = 0.01
MAX[igm] = 5
MIN[imz] = 0.01
MAX[imz] = 0.8
MIN[iKN] = 0.01
MAX[iKN] = 3
NB   = length(beta)
NPS  = 20
gm   = MIN[igm] + (MAX[igm] - MIN[igm])/(NPS-1)*(0:(NPS-1))
mz   = MIN[imz] + (MAX[imz] - MIN[imz])/(NPS-1)*(0:(NPS-1))
KN   =(MIN[iKN] + (MAX[iKN] - MIN[iKN])/(NPS-1)*(0:(NPS-1))) * A

#Read stabfile:
pdf('stab.pdf', width=8,height=8,paper='a4')
op <- par(font.lab = 1,
             family ="serif",
             mar    = c(4,4,1.,0.5),
             mgp    = c(2.3,1,0),
             oma    = c(4,4,4,0),
             mfcol  = c(3, 3)) 
 
for (i in 1:NB){
    for (j in 1:NPS){
         pf   = paste0('KN',   sprintf('%1.2f',KN[j]),
                       'beta', sprintf('%1.2f',beta[i]))
         out  = paste0(pf, '.out')
         stb  = paste0(pf, '.stb')
         dat  = read.table(stb, header = T)
         #Convert to matrix
         stb  = matrix(dat[,3], nr = NPS, nc = NPS)

         #Be careful of all are stable (=1)
         if (all(stb == 1)) {
            image(mz, gm, stb, col = gray(c(1)), main=pf)
         }else if(all(stb == 0)){
            image(mz, gm, stb, col = gray(c(0)), main=pf)
         }else{
            image(mz, gm, stb, col = gray(c(0,1)), main=pf)
         }
    }
}
dev.off()

#Plot time series of individual simulations
plotT(beta = 0.1, gm = 0.80, mz = 0.09)
plotT <- function(beta = 0.01, gm = 0.01, mz = 0.01,KN = 0.96){

    b = sprintf('%1.2f', beta)
    g = sprintf('%1.2f', gm)
    m = sprintf('%1.2f', mz)
    K = sprintf('%1.2f', KN)

    pl = paste0('KN',K,'beta',b)
    file = paste0(pl, '.out') 

    if (!file.exists(file)) stop("File does not exist!")
    dat  = read.table(file, header = T)
    for (i in 1:ncol(dat)) dat[,i] = as.numeric(as.character(dat[,i]))
    dat  = na.omit(dat)
    f    = dat[abs(dat$Gm - as.numeric(g)) < 0.005 
             & abs(dat$mz - as.numeric(m)) < 0.005,]
    file1 = paste0(pl,'Gm',g,'mz',m,'.pdf')
    pdf(file1, width=5,height=8,paper='a4')
    op <- par(font.lab = 1,
                 family ="serif",
                 mar    = c(4,4,1.5,0.5),
                 mgp    = c(2.3,1,0),
                 oma    = c(4,4,4,0),
                 mfcol  = c(3, 1)) 
    
    cff = c(f$NO3, f$PHY, f$ZOO)
    if (all(is.na(cff))) {
      #plot nothing
        plot(1:10, type = 'n')         
    }else{
        plot(f$Days, f$NO3, 
             ylim = c(min(cff, na.rm = T), max(cff, na.rm = T)),
             xlab = 'Days', ylab = 'N, P, Z', type = 'l')
        points(f$Days, f$PHY, type = 'l', col = 2)
        points(f$Days, f$ZOO, type = 'l', col = 3)
        legend('bottomright', c('NO3', 'PHY', 'ZOO'), lty = 1, col=1:3)
    }
    
    cff = c(f$VNO3, f$VPHY, f$VZOO)
    if (all(is.na(cff))) {
      #plot nothing
        plot(1:10, type = 'n')         
    }else{
        plot(f$Days, f$VNO3, 
             ylim = c(min(cff, na.rm = T), max(cff, na.rm = T)),
             xlab = 'Days', ylab = 'VAR of N, P, Z', type = 'l')
        points(f$Days, f$VPHY, type = 'l', col = 2)
        points(f$Days, f$VZOO, type = 'l', col = 3)
        legend('bottomright', c('VNO3', 'VPHY', 'VZOO'), lty = 1, col=1:3)
        abline(0,0, lty = 3)
    }   
    cff = c(f$COVNP, f$COVNZ, f$COVPZ)
    if (all(is.na(cff))) {
      #plot nothing
      plot(1:10, type = 'n')         
    }else{
      plot(f$Days, f$COVNP, 
           ylim = c(min(cff, na.rm = T), max(cff, na.rm = T)),
           xlab = 'Days', ylab = 'Covariances', type = 'l')
      points(f$Days, f$COVNZ, type = 'l', col = 2)
      points(f$Days, f$COVPZ, type = 'l', col = 3)
      abline(0,0, lty = 3)
      legend('bottomright', c('COVNZ', 'COVPZ', 'COVNP'), lty = 1, col=1:3)
    }
    dev.off()
}
