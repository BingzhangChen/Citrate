beta = c(0.01, 0.1,0.5, 1, 2)
igm  = 1
imz  = 2
NPAR = imz
MAX  = numeric(NPAR)
MIN  = MAX
MIN[igm] = 0.01
MAX[igm] = 5
MIN[imz] = 0.01
MAX[imz] = 0.8
NB   = length(beta)
NPS  = 50
gm   = MIN[igm] + (MAX[igm] - MIN[igm])/(NPS-1)*(0:(NPS-1))
mz   = MIN[imz] + (MAX[imz] - MIN[imz])/(NPS-1)*(0:(NPS-1))
stab = array(T, c(NB, NPS, NPS))
for (i in 1:NB){
    for (j in 1:NPS){
        for (u in 1:NPS){
            file = paste0('beta', sprintf('%1.2f', beta[i]),
                          'Gm',   sprintf('%1.2f',gm[j]),
                          'mz',   sprintf('%1.2f',mz[u]))
            dat  = read.table(file, header = T)
            for (k in 1:ncol(dat)) dat[,k] = as.numeric(as.character(dat[,k]))
            cff  = na.omit(dat)
            if (nrow(dat) > nrow(cff) || 
                any(c(dat$NO3, dat$PHY, dat$ZOO, dat$VNO3, dat$VPHY, dat$VZOO) <= 0, na.rm = T) )
                stab[i,j,u] = F
        }
    }
}

