#Plot examples of different ecotypes of different traits 
#Growth ~ temperature curve for 4 ecotypes with different Topt
source('~/Working/FlexEFT1D/Rscripts/NPZDcont.R')

Topt = c(0, 10, 20, 30)
NT   = length(Topt)

tC   = seq(0,30,.1)  #Environmental temperature
N    = length(tC)

#Assume size is 1 µm
LNV  = log(pi/6*1^3)

#Assume Iopt is 200 umol photons m-2 s-1
Iopt = log(200)
gr   = matrix(NA, nr = NT, nc = N)
for (i in 1:NT){
    gr[i,] = mu3(x=Iopt, Z = Topt[i]+273.15, t = tC+273.15, L=LNV)$mu
}
umin = min(gr)
umax = max(gr)
plot(tC, gr[1,], 
     ylim = c(umin, umax),
     type = 'l',
     xlab = 'Temperature (ºC)',
     ylab = 'Growth rate (d–1)')
for (i in 2:NT){
  points(tC,gr[i,],type = 'l', col=i)
}

#Plot growth rates of different Topt under the same temperature
Topt = seq(0,30,.1)  
NT   = length(Topt)

tC   = c(0, 10, 20, 30)#Environmental temperature
N    = length(tC)
gr   = matrix(NA, nr = NT, nc = N)
for (i in 1:N){
    gr[,i] = mu3(x=Iopt, Z = Topt+273.15, t = tC[i]+273.15, L=LNV)$mu
}
umin = min(gr)
umax = max(gr)
plot(Topt, gr[,1], 
     ylim = c(umin, umax),
     type = 'l',
     xlab = 'Optimal temperature (ºC)',
     ylab = 'Growth rate (d–1)')
for (i in 2:N){
  points(Topt,gr[,i],type = 'l', col=i)
}



