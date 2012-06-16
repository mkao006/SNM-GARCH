########################################################################
## This is a rewrite of the log.mgarch function in which the mean
## estimated seperately and the density is reparametised
########################################################################

## Load library
library(compiler)
library(abind)
library(fGarch)


## Parameters:
## param[1] = omega
## param[2] = beta
## param[3] = alpha
## param[4] = gamma
## param[5] = sigma1

## TODO (Michael): Write print, predict function for the mgarch class.
##
## TODO (Michael): Handle the class in a better way. Don't like coercing
##                 to numeric.
## TODO (Michael): Improve the efficiency of this code (e.g. remove outer)
##
## NOTES (Michael): The use of heaviside function is so that there are
##                  not undefined points as would have happen if we used
##                  ifelse or something similar. However, will need to
##                  use something else as the heaviside gives equal
##                  weight to the left and the right hand side which is
##                  appropriate consider our non-linear function.


## Single point estimation initialisation (This is a special case when
## smoothing window takes the value of 1)
##
## logd.mgarch <- function(xt, beta, pt, which){
##   ## Initialise variables and list
##   xt <- as.numeric(xt)
##   T <- length(xt)
##   lpt <- length(pt)
##   lb <- length(beta)
##   dl <- vector("list", length = 4)
##   names(dl) <- c("ld", "db", "dt")
##
##   ## Calculate the conditional variance
##   betaSum <- as.numeric(filter(xt[-T]^2, beta[2], "recursive"))
##     sigma.t <-
##         sqrt(c(beta[5]^2,
##                beta[1] * (1 - cumprod(rep.int(beta[2], T - 1)))/
##                (1 - beta[2]) +
##                cumprod(rep.int(beta[2], T - 1)) * beta[5]^2 +
##                beta[3] * betaSum))
##
##   ## Calculate the density
##   if(which[1] == 1){
##     dl$ld = log(2) - log(beta[4] + 1/beta[4]) -
##       0.5 * log(2 * pi) - log(outer(sigma.t, pt)) -
##       ((xt * beta[4])^2 * Heaviside(-xt) +
##        (xt / beta[4])^2 * Heaviside(xt))/(2 * outer(sigma.t^2, pt^2))
##   }
##
##   ## Calculate the derivatives
##   if(which[2] == 1){
##
##     ## Piece wise Analytical derivative
##     dldsigma = -1/sigma.t +
##       ((xt * beta[4])^2 * Heaviside(-xt) +
##        (xt / beta[4])^2 * Heaviside(xt))/outer(sigma.t^3, pt^2)
##
##     sig.vec <- 2 * sigma.t
##     convFilter <- 1:(T - 2) * cumprod(c(1, rep(beta[2], T - 3)))
##     myconvolve <- beta[3] * convolve(convFilter, rev(xt^2),
##                                      type = "open")[1:(T - 2)]
##     cp.beta <- c(0, 0, myconvolve)
##
##     ## beta[1] - Mean of the conditional variance equation
##     dsigmadalpha0 <-
##       c(0, ((1 - cumprod(rep.int(beta[2], T - 1))))/(1 - beta[2]))
##
##     ## beta[2] - Coefficient for lagged variance
##     dsigmadalpha1 <-
##       (beta[1] * cumsum(c(0, 0, 1:(T - 2) * beta[2]^(0:(T - 3))))) +
##         c(0, (1:(T - 1)) * beta[2]^(0:(T - 2)) * beta[5]^2) +
##          cp.beta
##
##     ## beta[3] - Coefficient for lagged sqaured observation
##     dsigmadbeta1 <- c(0, betaSum)
##
##     ## beta[5] - The initial variance
##     dsigmadsigma <-
##       c(1, 2 * beta[5] * beta[2]^(1:(T - 1)))
##
##     ## beta[4] - Skewness parameter
##     dldxi = -((1 - 1/beta[4]^2)/(beta[4] + 1/beta[4])) +
##       (xt^2 / beta[4]^3 * Heaviside(xt) -
##        xt^2 * beta[4] * Heaviside(-xt))/outer(sigma.t^2, pt^2)
##
##     ## Combine everything into an array.
##     ##
##     ## TODO (Michael):Improve this part, it feels a little bit
##     ##                clumsy. Also try to avoid abind.
##     dldparams <- array(c(dldsigma * dsigmadalpha0/sig.vec,
##                          dldsigma * dsigmadalpha1/sig.vec,
##                          dldsigma * dsigmadbeta1/sig.vec,
##                          dldsigma * dsigmadsigma/c(1, sig.vec[-1])),
##                            dim = c(T, lpt, lb - 1))
##     dl$db <- abind(dldparams, dldxi, along = 3)
##     print(head(dl$db))
##   }
##   if(which[3] == 1){
##     dl$dt = -1/matrix(rep(pt, each = T), nc = lpt) +
##       ((xt * beta[4])^2 * Heaviside(-xt) +
##        (xt / beta[4])^2 * Heaviside(xt))/outer(sigma.t^2, pt^3)
##   }
##   dl
## }


## Smoothing window estimation initilisation (5 parameters)
## logd.mgarch <- function(xt, beta, pt, which){
##   k = attr(xt, "smoothWindow")
##   T <- length(xt)
##   lpt <- length(pt)
##   lb <- length(beta)
##   dl <- vector("list", length = 4)
##   names(dl) <- c("ld", "db", "dt")
## 
##   ## Calculate the conditional variance
##   betaSum <- c(filter(xt[c(k:(T - 1))]^2, beta[2], "recursive"))
##     sigma.t <-
##         sqrt(c(rep(beta[5]^2, k),
##                beta[1] * (1 - cumprod(rep.int(beta[2], T - k)))/
##                (1 - beta[2]) +
##                cumprod(rep.int(beta[2], T - k)) * beta[5]^2 +
##                beta[3] * betaSum))
## 
##   ## Calculate the density
##   if(which[1] == 1){
##     dl$ld = log(2) - log(beta[4] + 1/beta[4]) -
##       0.5 * log(2 * pi) - log(outer(sigma.t, pt)) -
##       ((xt * beta[4])^2 * Heaviside(-xt) +
##        (xt / beta[4])^2 * Heaviside(xt))/(2 * outer(sigma.t^2, pt^2))
##   }
## 
##   ## Calculate the derivatives
##   if(which[2] == 1){
## 
##     ## Piece wise Analytical derivative
##     dldsigma = -1/sigma.t +
##       ((xt * beta[4])^2 * Heaviside(-xt) +
##        (xt / beta[4])^2 * Heaviside(xt))/outer(sigma.t^3, pt^2)
## 
##     convFilter <- 1:(T - k - 1) * beta[2]^(0:(T - k - 2))
##     betaSum2 <- c(0, beta[3] * convolve(convFilter, rev(xt[c(k:(T - 1))]^2),
##                                         type = "open")[1:(T - k - 1)])
## 
##     ## beta[1] - Mean of the conditional variance equation
##     dsigmadalpha0 <-
##       c(rep(0, k), cumsum(beta[2]^(0:(T - k - 1))))
## 
##     ## beta[2] - Coefficient for lagged variance
##     dsigmadalpha1 <-
##       c(rep(0, k), ((beta[1] * cumsum(c(0:(T - k - 1) *
##                                         beta[2]^(-1:(T - k - 2))))) +
##         1:(T - k) * beta[2]^(0:(T - k - 1)) * beta[5]^2 + betaSum2))
## 
##     ## beta[3] - Coefficient for lagged sqaured observation
##     dsigmadbeta1 <- c(rep(0, k), betaSum)
## 
##     ## beta[5] - Initial variance
##     dsigmadsigmac <- c(rep(1, k), 2 * beta[5] * beta[2]^(1:(T - k)))
## 
##     ## beta[4] - Skewness parameter
##     dldgamma = -((1 - 1/beta[4]^2)/(beta[4] + 1/beta[4])) +
##       (xt^2 / beta[4]^3 * Heaviside(xt) -
##        xt^2 * beta[4] * Heaviside(-xt))/outer(sigma.t^2, pt^2)
## 
##     ## Combine everything
##     dl$db <- array(c(dldsigma * c(dsigmadalpha0/
##                          c(rep(1, k), 2 * sigma.t[-c(1:k)])),
##                      dldsigma * c(dsigmadalpha1/
##                          c(rep(1, k), 2 * sigma.t[-c(1:k)])),
##                      dldsigma * c(dsigmadbeta1/
##                          c(rep(1, k), 2 * sigma.t[-c(1:k)])),
##                      dldsigma * c(dsigmadsigmac/
##                          c(rep(1, k), 2 * sigma.t[-c(1:k)])),
##                      dldgamma),
##                    dim = c(T, lpt, lb))
##   }
##   if(which[3] == 1){
##     dl$dt = -1/matrix(rep(pt, each = T), nc = lpt) +
##       ((xt * beta[4])^2 * Heaviside(-xt) +
##        (xt / beta[4])^2 * Heaviside(xt))/outer(sigma.t^2, pt^3)
##   }
##   dl
## }

## Conditional variance estimation
## logd.mgarch <- function(xt, beta, pt, which){
##   sigma1 = attr(xt, "sigma1")
##   T <- length(xt)
##   lpt <- length(pt)
##   lb <- length(beta)
##   dl <- vector("list", length = 4)
##   names(dl) <- c("ld", "db", "dt")
## 
##   ## Calculate the conditional variance
##   ## sigma.t <- cond.sd(xt, beta)
##   betaSum <- as.numeric(filter(xt[-T]^2, beta[2], "recursive"))
##     sigma.t <-
##         sqrt(c(sigma1^2,
##                beta[1] * (1 - cumprod(rep.int(beta[2], T - 1)))/
##                (1 - beta[2]) +
##                cumprod(rep.int(beta[2], T - 1)) * sigma1^2 +
##                beta[3] * betaSum))
##   ## attr(xt, "init.sigma")[2] <- sigma.t[2]
## 
##   ## Calculate the density
##   if(which[1] == 1){
##     dl$ld = log(2) - log(beta[5] + 1/beta[5]) -
##       0.5 * log(2 * pi) - log(outer(sigma.t, pt)) -
##       ((xt * beta[5])^2 * Heaviside(-xt) +
##        (xt / beta[5])^2 * Heaviside(xt))/(2 * outer(sigma.t^2, pt^2))
##   }
## 
##   ## Calculate the derivatives
##   if(which[2] == 1){
## 
##     ## Piece wise Analytical derivative
##     dldsigma = -1/sigma.t +
##       ((xt * beta[5])^2 * Heaviside(-xt) +
##        (xt / beta[5])^2 * Heaviside(xt))/outer(sigma.t^3, pt^2)
## 
##     convFilter <- 1:(T - 2) * cumprod(c(1, rep(beta[2], T - 3)))
##     myconvolve <- beta[3] * convolve(convFilter, rev(xt^2),
##                                      type = "open")[1:(T - 2)]
##     cp.beta <- c(0, myconvolve)
## 
##     ## d[1] = theoretical back estimate
##     ## d[2] = constrained back esitmate
##     ## d[3] = unconditional estimate
##     ## beta[1] - Mean of the conditional variance equation
##     d1 = (1 - beta[2] - beta[3])^-1
##     dsigmadalpha0 <-
##       c(d1, ((1 - cumprod(rep.int(beta[2], T - 1))))/((1 - beta[2])))
## 
##     ## beta[2] - Coefficient for lagged variance
##     d2 = (1 - beta[2] - beta[3])^-2
##     dsigmadalpha1 <-
##       c(d2, ((beta[1] * cumsum(c(0, 1:(T - 2) * beta[2]^(0:(T - 3))))) +
##         1:(T - 1) * beta[2]^(0:(T - 2)) * sigma1^2 + cp.beta))
## 
##     ## beta[3] - Coefficient for lagged sqaured observation
##     d3 = (1 - beta[2] - beta[3])^-2
##     dsigmadbeta1 <- c(d3, betaSum)
## 
##     ## beta[5] - Skewness parameter
##     dldxi = -((1 - 1/beta[5]^2)/(beta[5] + 1/beta[5])) +
##       (xt^2 / beta[5]^3 * Heaviside(xt) -
##        xt^2 * beta[5] * Heaviside(-xt))/outer(sigma.t^2, pt^2)
## 
##     ## Combine everything into an array.
##     ##
##     ## TODO (Michael):Improve this part, it feels a little bit
##     ##                clumsy. Also try to avoid abind.
##     dldparams <- array(c(dldsigma * dsigmadalpha0/c(1, 2 * sigma.t[-1]),
##                          dldsigma * dsigmadalpha1/c(1, 2 * sigma.t[-1]),
##                          dldsigma * dsigmadbeta1/c(1, 2 * sigma.t[-1])),
##                          dim = c(T, lpt, lb - 1))
##     dl$db <- abind(dldparams, dldxi, along = 3)
##   }
##   if(which[3] == 1){
##     dl$dt = -1/matrix(rep(pt, each = T), nc = lpt) +
##       ((xt * beta[5])^2 * Heaviside(-xt) +
##        (xt / beta[5])^2 * Heaviside(xt))/outer(sigma.t^2, pt^3)
##   }
##   dl
## }

## Back estimation initialisation
## logd.mgarch <- function(xt, beta, pt, which){
##   sigma1 = attr(xt, "sigma1")
##   T <- length(xt)
##   lpt <- length(pt)
##   lb <- length(beta)
##   dl <- vector("list", length = 4)
##   names(dl) <- c("ld", "db", "dt")
## 
##   ## Calculate the conditional variance
##   ## sigma.t <- cond.sd(xt, beta)
##   betaSum <- as.numeric(filter(xt[-T]^2, beta[2], "recursive"))
##     sigma.t <-
##         sqrt(c(sigma1^2,
##                beta[1] * (1 - cumprod(rep.int(beta[2], T - 1)))/
##                (1 - beta[2]) +
##                cumprod(rep.int(beta[2], T - 1)) * sigma1^2 +
##                beta[3] * betaSum))
##   ## attr(xt, "init.sigma")[2] <- sigma.t[2]
## 
##   ## Calculate the density
##   if(which[1] == 1){
##     dl$ld = log(2) - log(beta[5] + 1/beta[5]) -
##       0.5 * log(2 * pi) - log(outer(sigma.t, pt)) -
##       ((xt * beta[5])^2 * Heaviside(-xt) +
##        (xt / beta[5])^2 * Heaviside(xt))/(2 * outer(sigma.t^2, pt^2))
##   }
## 
##   ## Calculate the derivatives
##   if(which[2] == 1){
## 
##     ## Piece wise Analytical derivative
##     dldsigma = -1/sigma.t +
##       ((xt * beta[5])^2 * Heaviside(-xt) +
##        (xt / beta[5])^2 * Heaviside(xt))/outer(sigma.t^3, pt^2)
## 
##     convFilter <- 1:(T - 2) * cumprod(c(1, rep(beta[2], T - 3)))
##     myconvolve <- beta[3] * convolve(convFilter, rev(xt^2),
##                                      type = "open")[1:(T - 2)]
##     cp.beta <- c(0, myconvolve)
## 
##     ## beta[1] - Mean of the conditional variance equation
##     ## d1 = -0.5 * ((sigma.t[2]^2 - beta[1] - beta[3] * xt[1]^2) * beta[2])^-0.5
##     d1 = 0
##     dsigmadalpha0 <-
##       c(d1, ((1 - cumprod(rep.int(beta[2], T - 1))))/((1 - beta[2])))
## 
##     ## beta[2] - Coefficient for lagged variance
##     ## d2 = -0.5 * (sigma.t[2]^2 - beta[1] - beta[3] * xt[1]^2)^0.5 * beta[2]^-1.5
##     d2 = 0
##     dsigmadalpha1 <-
##       c(d2, ((beta[1] * cumsum(c(0, 1:(T - 2) * beta[2]^(0:(T - 3))))) +
##         1:(T - 1) * beta[2]^(0:(T - 2)) * sigma1^2 + cp.beta))
## 
##     ## beta[3] - Coefficient for lagged sqaured observation
##     ## d3 = -0.5 * xt[1]^2 * ((sigma.t[2]^2 - beta[1] - beta[3] * xt[1]^2) *
##     ##                        beta[2])^-0.5
##     d3 = 0
##     dsigmadbeta1 <- c(d3, betaSum)
## 
##     ## beta[5] - Skewness parameter
##     dldxi = -((1 - 1/beta[5]^2)/(beta[5] + 1/beta[5])) +
##       (xt^2 / beta[5]^3 * Heaviside(xt) -
##        xt^2 * beta[5] * Heaviside(-xt))/outer(sigma.t^2, pt^2)
## 
##     ## Combine everything into an array.
##     ##
##     ## TODO (Michael):Improve this part, it feels a little bit
##     ##                clumsy. Also try to avoid abind.
##     dldparams <- array(c(dldsigma * dsigmadalpha0/c(1, 2 * sigma.t[-1]),
##                          dldsigma * dsigmadalpha1/c(1, 2 * sigma.t[-1]),
##                          dldsigma * dsigmadbeta1/c(1, 2 * sigma.t[-1])),
##                          dim = c(T, lpt, lb - 1))
##     dl$db <- abind(dldparams, dldxi, along = 3)
##   }
##   if(which[3] == 1){
##     dl$dt = -1/matrix(rep(pt, each = T), nc = lpt) +
##       ((xt * beta[5])^2 * Heaviside(-xt) +
##        (xt / beta[5])^2 * Heaviside(xt))/outer(sigma.t^2, pt^3)
##   }
##   dl
## }

## Full logd for all initialisation
logd.mgarch <- function(xt, beta, pt, which){
  initMethod = attr(xt, "initialisation")
  if(initMethod == "Smooth"){
    sigma1 = beta[5]
    k = attr(xt, "smoothWindow")
  } else {
    sigma1 = attr(xt, "sigma1")
    k = 1
  }
  T <- length(xt)
  lpt <- length(pt)
  lb <- length(beta)
  dl <- vector("list", length = 4)
  names(dl) <- c("ld", "db", "dt")
  
  ## Calculate the conditional variance
  betaSum <- c(filter(xt[c(k:(T - 1))]^2, beta[2], "recursive"))
    sigma.t <-
        sqrt(c(rep(sigma1^2, k),
               beta[1] * (1 - cumprod(rep.int(beta[2], T - k)))/
               (1 - beta[2]) +
               cumprod(rep.int(beta[2], T - k)) * sigma1^2 +
               beta[3] * betaSum))

  ## Calculate the density
  if(which[1] == 1){
    dl$ld = log(2) - log(beta[4] + 1/beta[4]) -
      0.5 * log(2 * pi) - log(outer(sigma.t, pt)) -
      ((xt * beta[4])^2 * Heaviside(-xt) +
       (xt / beta[4])^2 * Heaviside(xt))/(2 * outer(sigma.t^2, pt^2))
  }

  ## Calculate the derivatives
  if(which[2] == 1){

    ## Piece wise Analytical derivative
    dldsigma = -1/sigma.t +
      ((xt * beta[4])^2 * Heaviside(-xt) +
       (xt / beta[4])^2 * Heaviside(xt))/outer(sigma.t^3, pt^2)

    convFilter <- 1:(T - k - 1) * beta[2]^(0:(T - k - 2))
    betaSum2 <- c(0, beta[3] * convolve(convFilter, rev(xt[c(k:(T - 1))]^2),
                                        type = "open")[1:(T - k - 1)])

    ## beta[1] - Mean of the conditional variance equation
    if(initMethod == "Smooth"){
      b11 = rep(0, k)
    } else if(initMethod == "BackFilter"){
      b11 = 0
    } else if(initMethod == "Unconditional"){
      b11 = (1 - beta[2] - beta[3])^-1
    }    
    dsigmadalpha0 <-
      c(b11, cumsum(beta[2]^(0:(T - k - 1))))


    ## beta[2] - Coefficient for lagged variance
    if(initMethod == "Smooth"){
      b21 = rep(0, k)
    } else if(initMethod == "BackFilter"){
      b21 = 0
    } else if(initMethod == "Unconditional"){
      b21 = (1 - beta[2] - beta[3])^-2
    }
    dsigmadalpha1 <-
      c(b21, ((beta[1] * cumsum(c(0:(T - k - 1) *
                                        beta[2]^(-1:(T - k - 2))))) +
        1:(T - k) * beta[2]^(0:(T - k - 1)) * sigma1^2 + betaSum2))


    ## beta[3] - Coefficient for lagged sqaured observation
    if(initMethod == "Smooth"){
      b31 = rep(0, k)
    } else if(initMethod == "BackFilter"){
      b31 = 0
    } else if(initMethod == "Unconditional"){
      b31 = (1 - beta[2] - beta[3])^-2
    }    
    dsigmadbeta1 <- c(b31, betaSum)

    ## beta[4] - Skewness parameter
    dldgamma = -((1 - 1/beta[4]^2)/(beta[4] + 1/beta[4])) +
      (xt^2 / beta[4]^3 * Heaviside(xt) -
       xt^2 * beta[4] * Heaviside(-xt))/outer(sigma.t^2, pt^2)

    ## beta[5] - Initial variance
    if(initMethod == "Smooth"){
      dsigmadsigmac <- c(rep(1, k), 2 * sigma1 * beta[2]^(1:(T - k)))
    }


    ## Combine everything
    if(initMethod == "Smooth"){
    dl$db <- array(c(dldsigma * c(dsigmadalpha0/
                         c(rep(1, k), 2 * sigma.t[-c(1:k)])),
                     dldsigma * c(dsigmadalpha1/
                         c(rep(1, k), 2 * sigma.t[-c(1:k)])),
                     dldsigma * c(dsigmadbeta1/
                         c(rep(1, k), 2 * sigma.t[-c(1:k)])),
                     dldgamma,
                     dldsigma * c(dsigmadsigmac/
                         c(rep(1, k), 2 * sigma.t[-c(1:k)]))),
                   dim = c(T, lpt, lb))
    } else {
    dl$db <- array(c(dldsigma * c(dsigmadalpha0/
                         c(rep(1, k), 2 * sigma.t[-c(1:k)])),
                     dldsigma * c(dsigmadalpha1/
                         c(rep(1, k), 2 * sigma.t[-c(1:k)])),
                     dldsigma * c(dsigmadbeta1/
                         c(rep(1, k), 2 * sigma.t[-c(1:k)])),
                     dldgamma), dim = c(T, lpt, lb))
    }
  }
  if(which[3] == 1){
    dl$dt = -1/matrix(rep(pt, each = T), nc = lpt) +
      ((xt * beta[4])^2 * Heaviside(-xt) +
       (xt / beta[4])^2 * Heaviside(xt))/outer(sigma.t^2, pt^3)
  }
  dl
}

## Compile the function to increase the speed
logd.mgarch <- cmpfun(logd.mgarch)

## Test whether the parameters are valid
valid.mgarch <- function(x, beta, mix){
    if(attr(x, "initialisation") == "Smooth"){
        beta[1] > 0 &&
        beta[2] >= 0 && beta[2] <= 1 &&
        beta[3] >= 0 && beta[3] <= 1 &&
        beta[4] > 0 &&
        beta[5] > 0
    } else {
        beta[1] > 0 &&
        beta[2] >= 0 && beta[2] <= 1 &&
        beta[3] >= 0 && beta[3] <= 1 &&
        beta[4] > 0
    }
}


initial.mgarch <- function(x, beta = NULL, mix = NULL, kmax = NULL){
  initMethod = attr(x, "initialisation")
  if(is.null(beta)){
    gf <- try(garchFit(data = as.numeric(x), trace = FALSE,
                       include.mean = FALSE, cond.dist = "snorm"))
    if(initMethod == "Smooth"){
      if(!(inherits(gf, "try-error"))){
        cgf <- coef(gf)
        beta <- c(cgf["omega"], cgf["beta1"], cgf["alpha1"], cgf["skew"],
                  gf@sigma.t[1])
      } else {
        beta <- c(1e-5, 0.8, 0.1, 1, sd(x))
      }
    } else {
      if(!(inherits(gf, "try-error"))){
        cgf <- coef(gf)
        beta <- c(cgf["omega"], cgf["beta1"], cgf["alpha1"], cgf["skew"])
      } else {
        beta <- c(1e-5, 0.8, 0.1, 1)
      }
    }
  } else {
    beta = beta
  }
  if(initMethod == "Smooth"){
    names(beta) <- c("omega", "beta1", "alpha1", "xi", "sigma1")
  } else {
    names(beta) <- c("omega", "beta1", "alpha1", "xi")
  }
  if(is.null(mix)){
    mix <- disc(1, 1)
  } else {
    mix = mix
  }
  list(beta = beta, mix = mix)
}


## The weights function, this is useful when the data is discrete and
## there are duplicated.
weights.mgarch <- function(x) 1

## Function to determine the grid
gridpoints.mgarch <- function(x, beta, grid){
  seq(1e-3, 20, length = grid)
}

## Function to restrict the space of the support points
suppspace.mgarch <- function(x, beta, mix){
  c(1e-3, 20)
}


## Function for converting different class of time series to mgarch
##
as.mgarch <- function(x, window = 5,
                      init = c("Smooth", "BackFilter", "Unconditional")){
    init = match.arg(init)
    mgts <- as.numeric(data.matrix(x))
    class(mgts) = "mgarch"
    attr(mgts, "initialisation") = init
    attr(mgts, "sigma1") = abs(x[1]) # For back estimate
    attr(mgts, "smoothWindow") = window # For smooth initial
    mgts
}


## summary result of the solution
summary.mgarch <- function(sol, digits = 4){
  cat(paste("Solution Standard Deviation: ",
            round(sqrt(sum(((sol$beta[4]^3 + 1/sol$beta[4]^3))/
                           (sol$beta[4] + 1/sol$beta[4]) *
                           sol$mix$pr * sol$mix$pt^2)), 4),
            "\n", sep = ""))
  cat(paste("Persistence: ", round(sum(sol$beta[2:3]), 4), "\n", sep = ""))
  cat(paste("Maximum gradient: ", round(max(abs(sol$grad)), 7), "\n", sep = ""))
  cat(paste("Log-likelihood: ", round(sol$ll, 4), "\n", sep = ""))
}

## Function for calculating the conditional standard deviation
cond.sd <- function(x, beta){
  T <- length(x)
  if(attr(x, "initialisation") == "Smooth"){
    sigma1 = beta[5]
    k = attr(x, "smoothWindow")
  } else {
    sigma1 = attr(x, "sigma1")
    k = 1
  }
  betaSum <- c(filter(x[c(k:(T - 1))]^2, beta[2], "recursive"))
  sigma.t <-
    sqrt(c(rep(sigma1^2, k),
           beta[1] * (1 - cumprod(rep.int(beta[2], T - k)))/
           (1 - beta[2]) +
           cumprod(rep.int(beta[2], T - k)) * sigma1^2 +
           beta[3] * betaSum))
  sigma.t
}


## NOTE (Michael): Sometimes the gradient are not zero because the
##                 likelihood does not change anymore. This function
##                 checks whether the log-likelihood is maximised by
##                 altering the coefficient and check the
##                 log-likelihood around its neighbourhood.
gradCheck <- function(x, sol){
  for(i in 1:length(sol$beta)){
    null <- rep(0, length(sol$beta))
    if(i == 1){
      null[i] = 1e-8
    } else {
      null[i] <- 1e-2
    }
    max.est <- all(logLik.snpmle(x, sol$beta, sol$mix)[1] >
                   logLik.snpmle(x, sol$beta + null, sol$mix)[1],
                   logLik.snpmle(x, sol$beta, sol$mix)[1] >
                   logLik.snpmle(x, sol$beta - null, sol$mix)[1])
    cat(paste("Estimation of ", names(sol$beta)[i], " is maximum: ",
              max.est, "\n", sep = ""))
  }
}


##' Prediction function
##'
##' The generic function for generating n-step forecast for a
##' SNM-GARCH(1, 1) model.
##'
##' @param sol The solution from the SNM-GARCH model
##' @param n.ahead The number of step of forecast.
##' @return The n-step forecast of the model
##' @export
##' @examples
##'
predict.mgarch <- function(sol, n.ahead = 1){
  h1 <- sqrt(sol$beta[1] + sol$beta[3] * sol$x[length(sol$x)]^2 +
             mybeta[2] * sol$sigma.t[length(x)]^2)
  if(n.ahead > 1){
    persist = sol$beta[2] + sol$beta[3]
    h = (sol$beta[1] * (1 - (persist)^(1:(n.ahead - 1))))/(1 - persist) +
      persist^(1:(n.ahead - 1)) * h1
  } else {
    h = NULL
  }
  c(h1, h)
}
