beta = c(0.01, .1, 1, 2, 5)
NB   = length(beta)
NPS  = 50
gm   = 0.1 + 4.9/(NPS-1)*c(0:(NPS-1))
stable = matrix(T, nr = NB, nc = NPS)

