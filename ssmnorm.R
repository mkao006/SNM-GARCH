########################################################################
## This script contains function related to the distribution function
## of the skewed scale normal mixture (SSMN)
########################################################################

## The variance of the mixture normal distribution is sqrt(mix$pr *
## mix$pt^2)

## NOTE (Michael): Error function taken from the VGAM package
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

ierf <- function(x){
    qnorm((1 + x)/2)/sqrt(2)
}

dsmnorm <- function(x, varmix = disc(1, 1)){
    n = length(x)
    n.mix = length(varmix$pt)
    d = 1/sqrt(2 * pi * matrix(rep(varmix$pt, each = n), nc = n.mix)) *
        exp(-0.5 * x^2/matrix(rep(varmix$pt, each = n), nc = n.mix))
    c(d %*% varmix$pr)
}

dssmnorm <- function(x, varmix = disc(1, 1), xi = 1){
    z = x * (xi * Heaviside(-x) + 1/xi * Heaviside(x))
    (2/(xi + 1/xi)) * dsmnorm(z, varmix = varmix)
}

## xi = 1.5
## m1 = sqrt(2/pi)
## m2 = 1
## mu = m1 * (xi - 1/xi)
## sigma = sqrt((m2 - m1^2) * (xi^2 + 1/xi^2) + 2 * m2^2 - m1)
## curve(dmsnorm(x), -3, 3, n = 500)
## curve(dnorm(x), add = TRUE, col = "red")
## curve(dsnorm(x, xi = 1.5), add = TRUE, col = "red")
## curve(dmsnorm((x - mu), xi = 1.5)*sigma,
##       add = TRUE, col = "blue")
## curve(dmsnorm2(x, xi = 1.5), add = TRUE, col = "green")
## curve(dmsnorm(x, xi = 3, varmix = disc(c(0.2, 0.8), c(0.5, 0.5))),
##       add = TRUE)


psmnorm <- function(q, varmix = disc(1, 1)){
  z = outer(q, 1/sqrt(2 * varmix$pt^2))
  d = 0.5 * (1 + erf(z))
  c(d %*% varmix$pr)
}

## NOTE (Michael): Old and deprecated
## pssmnorm <- function(q, varmix = disc(1, 1), xi = 1){
##   Xi = xi^sign(q)
##   g = 2/(xi + 1/xi)
##   Probability = Heaviside(q) - sign(q) * g * Xi *
##     psmnorm(q = -abs(q)/Xi, varmix = varmix)
##   Probability
## }
##
## pssmnorm2 <- function(q, varmix = disc(1, 1), xi = 1){
##   Xi = xi^sign(-q)
##   g = 1/(xi + 1/xi)
##   qntl = g * xi^-1 * (erf(outer(q * Xi, 1/sqrt(2 * varmix$pt^2))) + 1)
##   qntl %*% varmix$pr
## }

pssmnorm <- function(q, varmix = disc(1, 1), xi = 1){
  g = 2/(xi + 1/xi)
  H = xi * Heaviside(-q) + 1/xi * Heaviside(q)
  P = g * 1/sqrt(2) * (erf(outer(q * H, 1/sqrt(2 * varmix$pt^2)))/sqrt(2 * H^2) +
  1/sqrt(2 * xi^2))
  c(P %*% varmix$pr)
}

## curve(pssmnorm(x), -5, 5, ylim = c(0, 1))
## curve(pssmnorm3(x), add = TRUE, col = "red", lty = 2)
## curve(pssmnorm(x, xi = 1.5), add = TRUE, n = 1000)
## curve(pssmnorm3(x, xi = 1.5), add = TRUE, col = "red", n = 1000, lty = 2)
## curve(pssmnorm(x, xi = 5), add = TRUE, n = 1000)
## curve(pssmnorm3(x, xi = 5), add = TRUE, col = "red", n = 1000, lty = 2)
## curve(pssmnorm3(x, varmix = disc(pt = c(0.5, 0.8), pr = c(0.5, 0.5)), xi = 1.5),
##       add = TRUE, col = "blue", n = 1000, lty = 2)
## curve(pssmnorm3(x, varmix = disc(pt = c(1.5, 3.8), pr = c(0.5, 0.5)), xi = 1.5),
##       add = TRUE, col = "blue", n = 1000, lty = 2)
## abline(h = 1/(1.5^2 + 1), col = "green", lty = 3)


qsmnorm <- function(p, varmix = disc(1, 1)){
  outer(ierf(2 * p - 1), sqrt(2 * varmix$pt^2)) %*% varmix$pr
}

## NOTE (Michael): This is a close solution but something is wrong
##                with it.

## qssmnorm <- function(p, varmix = disc(1, 1), xi = 1){
## ##    a = 1/(1 + 1/xi^2)
##     a = 0
##     H = xi * Heaviside(-p, a = a) + 1/xi * Heaviside(p, a = a)
##     F = (ierf(p * (xi + 1/xi) * H) - 1/xi) *
##         outer(1/H, sqrt(2 * varmix$pt^2))
##     F %*% varmix$pr
## }

## curve(qnorm(x), -1, 1, ylim = c(-5, 5), n = 100000)
## curve(qsnorm(x, xi = 1), col = "blue", add = TRUE, n = 100000)
## curve(qssmnorm(x), add = TRUE, col = "red", n = 100000)

## curve(qssmnorm(x, , varmix = disc(pt = c(1.5, 3.8), pr = c(0.5, 0.5))),
##       add = TRUE, col = "red", n = 100000)
## curve(qsmnorm(x,varmix = disc(pt = c(1.5, 3.8), pr = c(0.5, 0.5))),
##       add = TRUE, col = "blue", n = 100000, lty = 2)





## This is a bad implementation, but we will work with it for now.
qssmnorm <- function(p, varmix = disc(1, 1), xi = 1){
  prob = double(length(p))
  for(i in 1:length(p)){
    prob[i] = uniroot(function(x) pssmnorm(x, varmix = varmix, xi = xi) - p[i],
          interval = c(-20, 20))$root
  }
  prob
}


## NOTE (Michael): The following quantile function is instable and
## requires resolve
##
## qmsnorm <- function(p, varmix = disc(1, 1), xi = 1){
##   g = 2/(xi + 1/xi)
##   sig = sign(p - 1/2)
##   Xi = xi^sig
##   p = (Heaviside(p - 1/2) - sig * p)/(g * Xi)
##   Quantile = (-sig * qmnorm(p, varmix = varmix))
##   Quantile
## }
## curve(qnorm(x), 0, 1)
## curve(qmsnorm(x), add = TRUE, col = "red", lty = 3)
## curve(qmnorm(x, varmix = disc(pt = c(0.5, 2), pr = c(0.5, 0.5))),
##       add = TRUE, col = "blue")
## curve(qmsnorm(x, xi = 2), add = TRUE, col = "blue")
## curve(qmsnorm(x, varmix = disc(pt = c(0.5, 2), pr = c(0.5, 0.5)), xi = 2),
##       add = TRUE, col = "blue")



rsmnorm <- function(n, mix = disc(1, 1)){
    ## Generate a sample of length N from U(0, 1) to determine which
    ## distribution should the sample drawn from
    w <- runif(n)

    ## Select the right component distribution and generate the random
    ## sample
    breaks = cumsum(mix$pr)

    ## TODO: Try vectorize this
    m <- double(n)
    for(i in 1:n){
        m[i] = sum(w[i] > breaks) + 1
    }
    sim <- rnorm(n, mean = 0, sd = mix$pt[m])
    sim
}

## NOTE (Michael): This is incorrect, it is a mixture of skewed rather
##                 than the skewed mixture
rssmnorm <- function(n, mix = disc(1, 1), xi = 1){
    ## Generate a sample of length N from U(0, 1) to determine which
    ## distribution should the sample drawn from
    w <- runif(n)

    ## Select the right component distribution and generate the random
    ## sample
    breaks = cumsum(mix$pr)

    ## TODO: Try vectorize this
    m <- double(n)
    for(i in 1:n){
        m[i] = sum(w[i] > breaks) + 1
    }
    sim <- rsnorm(n, mean = 0, sd = mix$pt[m], xi = xi)
    sim
}


## A function for calculating the moments of a skewed scaled normal
## mixture
mssmnorm <- function(varmix = disc(1, 1), xi = 1, moments = 1){
  if(moments == 1){
    m = (sqrt(2 * varmix$pt/pi) * (xi^2 - xi^-2)/(xi + xi^-1)) %*% varmix$pr
  } else if(moments == 2){
    m = (varmix$pt * (xi^3 + xi^-3)/(xi + xi^-1)) %*% varmix$pr
  }
  as.numeric(m)
}
  
