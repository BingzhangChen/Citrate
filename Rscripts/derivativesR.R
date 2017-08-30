mu0 = expression( mu0hat * exp(alphamu*L + betamu*L**2))
DD <- function(expr, name, order = 1) {
        if(order < 1) stop("'order' must be >= 1")
        if(order == 1) D(expr, name)
        else DD(D(expr, name), name, order - 1)
}
dmu0dL = D(mu0, 'L')
d2mu0dL2 = D(dmu0dL,'L')
d3mu0dL3 = D(d2mu0dL2,'L')
d4mu0dL4 = D(d3mu0dL3,'L')
eval(dmu0dL, list(mu0hat=1, alphamu=0.1,betamu=-0.01,L=1))
eval(d2mu0dL2, list(mu0hat=1, alphamu=0.1,betamu=-0.01,L=1))
eval(d4mu0dL4, list(mu0hat=1, alphamu=0.1,betamu=-0.01,L=1))

fN = expression(NO3/(NO3 + K0N* exp(alphaK*L)))
dfNdl = D(fN,'L')
d2fNdl2=D(dfNdl,'L')
d3fNdl3=D(d2fNdl2,'L')
d4fNdl4=D(d3fNdl3,'L')
eval(fN, list(K0N=1, NO3=1,alphaK=0.3, L=1))
eval(d2fNdl2, list(K0N=1, NO3=1,alphaK=0.3, L=1))
eval(d3fNdl3, list(K0N=1, NO3=1,alphaK=0.3, L=1))
eval(d4fNdl4, list(K0N=1, NO3=1,alphaK=0.3, L=1))

aI_mu0 = expression(aI0 * exp(alphaI*L)/(mu0hat*exp(alphamu*L+betamu*L**2)))
daI_mu0dl = D(aI_mu0, 'L')
d2aI_mu0dl2 = D(daI_mu0dl, 'L')
d3aI_mu0dl3 = D(d2aI_mu0dl2, 'L')
d4aI_mu0dl4 = D(d3aI_mu0dl3, 'L')

eval(daI_mu0dl, list(aI0=0.05, PAR=100,alphaI=-0.1, L=1,mu0hat=1, alphamu=0.1,betamu=-0.01))

eval(d2aI_mu0dl2, list(aI0=0.05, PAR=100,alphaI=-0.1, L=1,mu0hat=1, alphamu=0.1,betamu=-0.01))

eval(d3aI_mu0dl3, list(aI0=0.05, PAR=100,alphaI=-0.1, L=1,mu0hat=1, alphamu=0.1,betamu=-0.01))

eval(d4aI_mu0dl4, list(aI0=0.05, PAR=100,alphaI=-0.1, L=1,mu0hat=1, alphamu=0.1,betamu=-0.01))

aI_mu02 = expression(aI0 * exp(alphaI*L)/(mu0hat*exp(alphamu*L+betamu*L**2))**2)
daI_mu02dl = D(aI_mu02,'L')
d2aI_mu02dl2 = D(daI_mu02dl, 'L')

eval(daI_mu02dl, list(aI0=0.05, PAR=100,alphaI=-0.1, L=1,mu0hat=1, alphamu=0.1,betamu=-0.01))
eval(d2aI_mu02dl2, list(aI0=0.05, PAR=100,alphaI=-0.1, L=1,mu0hat=1, alphamu=0.1,betamu=-0.01))

aI_mu03 = expression(aI0 * exp(alphaI*L)/(mu0hat*exp(alphamu*L+betamu*L**2))**3)
daI_mu03dl = D(aI_mu03,'L')

eval(daI_mu03dl, list(aI0=0.05, PAR=100,alphaI=-0.1, L=1,mu0hat=1, alphamu=0.1,betamu=-0.01))


SI = expression(1-exp(-aI0 * exp(alphaI*L)*PAR/(mu0hat*exp(alphamu*L+betamu*L**2))))

SI = expression(1-exp(-aI0 * exp(alphaI*L)*PAR/mu0hat))
dSIdL = D(SI,'L')
d2SIdL2 = D(dSIdL,'L')
d3SIdL3 = D(d2SIdL2,'L')
d4SIdL4 = D(d3SIdL3,'L')
eval(SI, list(aI0=0.05, PAR=100,alphaI=-0.1, L=1,mu0hat=1, alphamu=0.1,betamu=-0.01))
eval(dSIdL, list(aI0=0.05, PAR=100,alphaI=-0.1, L=1,mu0hat=1, alphamu=0.1,betamu=-0.01))
eval(d2SIdL2, list(aI0=0.05, PAR=100,alphaI=-0.1, L=1,mu0hat=1, alphamu=0.1,betamu=-0.01))
eval(d3SIdL3, list(aI0=0.05, PAR=100,alphaI=-0.1, L=1,mu0hat=1, alphamu=0.1,betamu=-0.01))
eval(d4SIdL4, list(aI0=0.05, PAR=100,alphaI=-0.1, L=2,mu0hat=1, alphamu=0.1,betamu=-0.01))

g = expression(P/(P**2 + Kp**2))

dgdL = D(g, 'P')


f = expression((Kp^2-P^2)/(P*(Kp^2+P^2)))
dfdL = D(f,'P')


f = expression((P-x)/(x+P))
dfdL = D(f,'x')

mu= expression( mu0hat * exp(alphamu*L + betamu*L**2) * (1-exp(-aI0 * exp(alphaI*L)*PAR/mu0hat)) *(NO3/(NO3 + K0N* exp(alphaK*L))) )
dmudL   = D(mu, 'L')
d2mudL2 = D(dmudL,'L')
d3mudL3 = D(d2mudL2,'L')
d4mudL4 = D(d3mudL3,'L')
eval(dmudL, list(mu0hat=1, alphamu=0.1,betamu=-0.01,L=1))
eval(d2mudL2, list(mu0hat=1, alphamu=0.1,betamu=-0.01,L=2, aI0=0.05,alphaI=-.1,PAR=100, NO3=.2, K0N=.2, alphaK=.27 ))
eval(d3mudL3, list(mu0hat=1, alphamu=0.1,betamu=-0.01,L=2, aI0=0.05,alphaI=-.1,PAR=100, NO3=.2, K0N=.2, alphaK=.27 ))
eval(d4mudL4, list(mu0hat=1, alphamu=0.1,betamu=-0.01,L=2, aI0=0.05,alphaI=-.1,PAR=100, NO3=.2, K0N=.2, alphaK=.27))

#Write out analytic equation:
mu = expression( 1/(I/(a0*exp(k*x)*exp(2*x)) + 1/(mu0*exp(b*x)) -2/(a0*exp(k*x)*exp(x)) + 1/(I*a0*exp(k*x))) * exp(alphamu*L + betamu*L**2) * NO3/(NO3 + K0N* exp(alphaK*L)) * exp(0.0633*(t - 15)) *(1-((t-Topt)/w)**2) )

#Here x = log(Iopt)
#Parameterize light values:
#Read Edward data:
setwd('~/Working/FlexEFT1D/Rscripts')
dat = read.csv('Edwards2015.csv')
dat = dat[dat$I_opt < 1000,]
dat$x = log(dat$I_opt)
dat$lnmu=log(dat$mu_max)
dat$lna =log(dat$alpha)
lm1=lm(lnmu~x,dat)  #mu0 = 0.21  d-1 at Iopt of 1 umol photons m-2 s-1, b = 0.21
lm2=lm(lna ~x,dat)  #a0  = 0.149 umol photons-1 m2 s d-1 at Iopt of 1 umol photons m-2 s-1,k = -0.47
cor.test(dat$I_opt,dat$mu_max)
plot(dat$x, log(dat$mu_max))
plot(dat$x, dat$lna)

eval(mu, list(I=200, a0=0.149, k = -0.47, mu0=.21, b = 0.21, x = log(1), alphamu=0.2, betamu = 0, L=0, NO3=10, K0N=.2, alphaK=.27, t=15, Topt=15, w=5))
