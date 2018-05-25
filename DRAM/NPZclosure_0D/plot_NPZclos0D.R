beta = c(0.01, .1, 1, 2, 5)
NB   = length(beta)
NPS  = 50
gm   = 0.1 + 4.9/(NPS-1)*c(0:(NPS-1))
stable = matrix(T, nr = NB, nc = NPS)
for (i in 1:NB){
  for (j in 1:NPS){
    b = sprintf('%1.2f', beta[i])
    g = sprintf('%1.2f', gm[j])
    file = paste0('beta',b,'Gm',g)
    dat  = read.table(file, header = T)
    b    = na.omit(dat)
    if (nrow(b) < nrow(dat)) stable[i,j] = F
    pdffile = paste0(file,'.pdf')
    pdf(pdffile, width=5,height=8,paper='a4')
    op <- par(font.lab = 1,
                 family ="serif",
                 mar    = c(4,4,1.5,0.5),
                 mgp    = c(2.3,1,0),
                 oma    = c(4,4,4,0),
                 mfcol  = c(3, 1)) 
    
    cff = c(dat$NO3, dat$PHY, dat$ZOO)
    if (all(is.na(cff))) {
      #plot nothing
        plot(1:10, type = 'n')         
    }else{
        plot(dat$Days, dat$NO3, 
             ylim = c(min(cff, na.rm = T), max(cff, na.rm = T)),
             xlab = 'Days', ylab = 'N, P, Z', type = 'l')
        points(dat$Days, dat$PHY, type = 'l', col = 2)
        points(dat$Days, dat$ZOO, type = 'l', col = 3)
        legend('bottomright', c('NO3', 'PHY', 'ZOO'), lty = 1, col=1:3)
    }
    
    cff = c(dat$VNO3, dat$VPHY, dat$VZOO)
    if (all(is.na(cff))) {
      #plot nothing
        plot(1:10, type = 'n')         
    }else{
        plot(dat$Days, dat$VNO3, 
             ylim = c(min(cff, na.rm = T), max(cff, na.rm = T)),
             xlab = 'Days', ylab = 'VAR of N, P, Z', type = 'l')
        points(dat$Days, dat$VPHY, type = 'l', col = 2)
        points(dat$Days, dat$VZOO, type = 'l', col = 3)
        legend('bottomright', c('VNO3', 'VPHY', 'VZOO'), lty = 1, col=1:3)
        abline(0,0, lty = 3)
    }   
    cff = c(dat$COVNP, dat$COVNZ, dat$COVPZ)
    if (all(is.na(cff))) {
      #plot nothing
      plot(1:10, type = 'n')         
    }else{
      plot(dat$Days, dat$COVNP, 
           ylim = c(min(cff, na.rm = T), max(cff, na.rm = T)),
           xlab = 'Days', ylab = 'Covariances', type = 'l')
      points(dat$Days, dat$COVNZ, type = 'l', col = 2)
      points(dat$Days, dat$COVPZ, type = 'l', col = 3)
      abline(0,0, lty = 3)
      legend('bottomright', c('COVNZ', 'COVPZ', 'COVNP'), lty = 1, col=1:3)
    }
    dev.off()
  }
}
