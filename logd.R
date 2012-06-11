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
## param[4] = sigma0
## param[5] = xi

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


## logd.mgarch <- function(xt, beta, pt, which){
##   ## xt <- as.numeric(xt)
##   ## Initialise variables and list
##   #k <- 20
##   ## sigma1 <- attr(xt, "sigma1")
##   sigma1 <- beta[4]
##   T <- length(xt)
##   lpt <- length(pt)
##   lb <- length(beta)
##   dl <- vector("list", length = 4)
##   names(dl) <- c("ld", "db", "dt")

##   ## Calculate the conditional variance
##   ## sigma.t <- cond.sd(xt, beta)
##   betaSum <- as.numeric(filter(xt[-T]^2, beta[2], "recursive"))
##     sigma.t <-
##         sqrt(c(sigma1^2,
##                beta[1] * (1 - cumprod(rep.int(beta[2], T - 1)))/
##                (1 - beta[2]) +
##                cumprod(rep.int(beta[2], T - 1)) * sigma1^2 +
##                beta[3] * betaSum))
  
##   ## Calculate the density
##   if(which[1] == 1){
##     dl$ld = log(2) - log(beta[5] + 1/beta[5]) -
##       0.5 * log(2 * pi) - log(outer(sigma.t, pt)) -
##       ((xt * beta[5])^2 * Heaviside(-xt) +
##        (xt / beta[5])^2 * Heaviside(xt))/(2 * outer(sigma.t^2, pt^2))
##   }

##   ## Calculate the derivatives
##   if(which[2] == 1){

##     ## Piece wise Analytical derivative
##     dldsigma = -1/sigma.t +
##       ((xt * beta[5])^2 * Heaviside(-xt) +
##        (xt / beta[5])^2 * Heaviside(xt))/outer(sigma.t^3, pt^2)

##     convFilter <- 1:(T - 2) * cumprod(c(1, rep(beta[2], T - 3)))
##     myconvolve <- beta[3] * convolve(convFilter, rev(xt^2),
##                                      type = "open")[1:(T - 2)]
##     cp.beta <- c(0, myconvolve)

##     ## beta[1] - Mean of the conditional variance equation    
##     dsigmadalpha0 <-
##       c(-0.5 * ((sigma.t[2]^2 - beta[1] -
##                  beta[2] * sigma.t[1]^2) * beta[3])^-0.5,
##         ((1 - cumprod(rep.int(beta[2], T - 1))))/
##         ((1 - beta[2])))
##     ##plot(dsigmadalpha0)

##     ## beta[2] - Coefficient for lagged variance
##     dsigmadalpha1 <-
##       c(-0.5 * sigma.t[1]^2 * ((sigma.t[2]^2 - beta[1] -
##                  beta[2] * sigma.t[1]^2) * beta[3])^-0.5,
##       ((beta[1] * cumsum(c(0, 1:(T - 2) * beta[2]^(0:(T - 3))))) +
##         1:(T - 1) * beta[2]^(0:(T - 2)) * beta[4]^2 +
##          cp.beta))

##     ## plot(dsigmadalpha1)
    
##     ## beta[3] - Coefficient for lagged sqaured observation
##     dsigmadbeta1 <- c(-0.5 * (sigma.t[2]^2 - beta[1] - beta[2] *
##                               sigma.t[1]^2)^0.5 * beta[3]^-1.5, betaSum)

##     ## beta[5] - The initial variance
##     dsigmadsigma <-
##       c(-beta[2] * sigma.t[1] * ((sigma.t[2]^2 - beta[1] -
##                  beta[2] * sigma.t[1]^2) * beta[3])^-0.5,
##         (2 * beta[4] * beta[2]^(1:(T - 1))))

##     ## beta[5] - Skewness parameter
##     dldxi = -((1 - 1/beta[5]^2)/(beta[5] + 1/beta[5])) + 
##       (xt^2 / beta[5]^3 * Heaviside(xt) -
##        xt^2 * beta[5] * Heaviside(-xt))/outer(sigma.t^2, pt^2)

##     ## Combine everything into an array.
##     ##
##     ## TODO (Michael):Improve this part, it feels a little bit
##     ##                clumsy. Also try to avoid abind.
##     dldparams <- array(c(dldsigma * dsigmadalpha0/c(1, 2 * sigma.t[-1]),
##                          dldsigma * dsigmadalpha1/c(1, 2 * sigma.t[-1]),
##                          dldsigma * dsigmadbeta1/c(1, 2 * sigma.t[-1]),
##                          dldsigma * dsigmadsigma/c(1, 2 * sigma.t[-1])),
##                            dim = c(T, lpt, lb - 1))
##     dl$db <- abind(dldparams, dldxi, along = 3)
##   }
##   if(which[3] == 1){
##     dl$dt = -1/matrix(rep(pt, each = T), nc = lpt) +
##       ((xt * beta[5])^2 * Heaviside(-xt) +
##        (xt / beta[5])^2 * Heaviside(xt))/outer(sigma.t^2, pt^3)
##   }
##   dl
## }

## NOTE (Michael): Then the problem is that the xt is not outputed

logd.mgarch <- function(xt, beta, pt, which){
  ## Initialise variables and list
  ## if(!is.null(attr(xt, "init.sigma")[2])){
  ##   attr(xt, "init.sigma")[1] = as.numeric(sqrt(abs((s[2]^2 - beta[1] -
  ##                           beta[2] * s[1]^2)/beta[3])))
  ## }
  ## sigma1 = attr(xt, "init.sigma")[1]
  sigma1 = attr(xt, "sigma1")
  T <- length(xt)
  lpt <- length(pt)
  lb <- length(beta)
  dl <- vector("list", length = 4)
  names(dl) <- c("ld", "db", "dt")

  ## Calculate the conditional variance
  ## sigma.t <- cond.sd(xt, beta)
  betaSum <- as.numeric(filter(xt[-T]^2, beta[2], "recursive"))
    sigma.t <-
        sqrt(c(sigma1^2,
               beta[1] * (1 - cumprod(rep.int(beta[2], T - 1)))/
               (1 - beta[2]) +
               cumprod(rep.int(beta[2], T - 1)) * sigma1^2 +
               beta[3] * betaSum))
  ## attr(xt, "init.sigma")[2] <- sigma.t[2]
  
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

    convFilter <- 1:(T - 2) * cumprod(c(1, rep(beta[2], T - 3)))
    myconvolve <- beta[3] * convolve(convFilter, rev(xt^2),
                                     type = "open")[1:(T - 2)]
    cp.beta <- c(0, myconvolve)

    ## beta[1] - Mean of the conditional variance equation
    ## d1 = -0.5 * ((sigma.t[2]^2 - beta[1] -
    ##               beta[2] * sigma.t[1]^2) * beta[3])^-0.5
    ##
    d1 = -0.5 * ((sigma.t[2]^2 - beta[1] - beta[3] * xt[1]^2) * beta[2])^-0.5
    d1 = 0
    dsigmadalpha0 <-
      c(d1, ((1 - cumprod(rep.int(beta[2], T - 1))))/((1 - beta[2])))

    ## beta[2] - Coefficient for lagged variance
    ## d2 = -0.5 * sigma.t[1]^2 * ((sigma.t[2]^2 - beta[1] -
    ##                              beta[2] * sigma.t[1]^2) * beta[3])^-0.5
    d2 = -0.5 * (sigma.t[2]^2 - beta[1] - beta[3] * xt[1]^2)^0.5 * beta[2]^-1.5
    d2 = 0
    dsigmadalpha1 <-
      c(d2, ((beta[1] * cumsum(c(0, 1:(T - 2) * beta[2]^(0:(T - 3))))) +
        1:(T - 1) * beta[2]^(0:(T - 2)) * sigma1^2 + cp.beta))
    
    ## beta[3] - Coefficient for lagged sqaured observation
    ## d3 = -0.5 * (sigma.t[2]^2 - beta[1] - beta[2] *
    ##   sigma.t[1]^2)^0.5 * beta[3]^-1.5
    d3 = -0.5 * xt[1]^2 * ((sigma.t[2]^2 - beta[1] - beta[3] * xt[1]^2) *
                           beta[2])^-0.5
    d3 = 0
    dsigmadbeta1 <- c(d3, betaSum)

    ## beta[4] - Skewness parameter
    dldxi = -((1 - 1/beta[4]^2)/(beta[4] + 1/beta[4])) + 
      (xt^2 / beta[4]^3 * Heaviside(xt) -
       xt^2 * beta[4] * Heaviside(-xt))/outer(sigma.t^2, pt^2)

    ## Combine everything into an array.
    ##
    ## TODO (Michael):Improve this part, it feels a little bit
    ##                clumsy. Also try to avoid abind.
    dldparams <- array(c(dldsigma * dsigmadalpha0/c(1, 2 * sigma.t[-1]),
                         dldsigma * dsigmadalpha1/c(1, 2 * sigma.t[-1]),
                         dldsigma * dsigmadbeta1/c(1, 2 * sigma.t[-1])),
                         dim = c(T, lpt, lb - 1))
    dl$db <- abind(dldparams, dldxi, along = 3)
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

## ## Test whether the parameters are valid
## valid.mgarch <- function(x, beta, mix){
##   beta[1] > 0 &&
##   beta[2] >= 0 && beta[2] <= 1 &&
##   beta[3] >= 0 && beta[3] <= 1 &&
##   beta[4] > 0 &&
##   beta[5] > 0
## }

## Test whether the parameters are valid
valid.mgarch <- function(x, beta, mix){
  beta[1] > 0 &&
  beta[2] >= 0 && beta[2] <= 1 &&
  beta[3] >= 0 && beta[3] <= 1 &&
  beta[4] > 0
}

## initial.mgarch <- function(x, beta = NULL, mix = NULL, kmax = NULL){
##   if(is.null(beta)){
##     gf <- try(garchFit(data = as.numeric(x), trace = FALSE,
##                        include.mean = FALSE, cond.dist = "snorm"))
##     if(!(inherits(gf, "try-error"))){
##       cgf <- coef(gf)
##       beta <- c(cgf["omega"], cgf["beta1"], cgf["alpha1"],
##                 gf@sigma.t[1], cgf["skew"])
##     } else {
##       beta <- c(1e-6, 0.8, 0.1, sd(x), 1)
##     }
##   } else {
##     beta = beta
##   }
##   names(beta) <- c("omega", "beta1", "alpha1", "sigma0", "xi")
##   if(is.null(mix)){
##     mix <- disc(1, 1)
##   } else {
##     mix = mix
##   }
##   list(beta = beta, mix = mix)
## }

initial.mgarch <- function(x, beta = NULL, mix = NULL, kmax = NULL){
  if(is.null(beta)){
    gf <- try(garchFit(data = as.numeric(x), trace = FALSE,
                       include.mean = FALSE, cond.dist = "snorm"))
    if(!(inherits(gf, "try-error"))){
      cgf <- coef(gf)
      beta <- c(cgf["omega"], cgf["beta1"], cgf["alpha1"], cgf["skew"])
    } else {
      beta <- c(1e-6, 0.8, 0.1, 1)
    }
  } else {
    beta = beta
  }
  names(beta) <- c("omega", "beta1", "alpha1", "xi")
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
  seq(0.1, 20, length = grid)
}

## Function to restrict the space of the support points
suppspace.mgarch <- function(x, beta, mix){
  c(0.1, 20)
}

## Function for converting different class of time series to mgarch
##
## TODO (Michael): Check John Chamber's book and all the time series
##                 classes in R. This is not a very nice solution, but
##                 just one thats convinient
as.mgarch <- function(x){
    mgts <- as.numeric(data.matrix(x))
    class(mgts) <- "mgarch"
    attr(mgts, "sd") = sd(x)
    attr(mgts, "sigma1") = abs(x[1])
    ## attr(mgts, "sigma1") = sd(x)
    mgts
}


## ## summary result of the solution
## summary.mgarch <- function(sol, digits = 4){
##   cat(paste("Solution Standard Deviation: ",
##             round(sqrt(sum(((sol$beta[5]^3 + 1/sol$beta[5]^3))/
##                            (sol$beta[5] + 1/sol$beta[5]) *
##                            sol$mix$pr * sol$mix$pt^2)), 4),
##             "\n", sep = ""))
##   cat(paste("Persistence: ", round(sum(sol$beta[2:3]), 4), "\n", sep = ""))
##   cat(paste("Maximum gradient: ", round(max(abs(sol$grad)), 7), "\n", sep = ""))
##   cat(paste("Log-likelihood: ", round(sol$ll, 4), "\n", sep = ""))
## }

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

## Extract coefficients
coef.mgarch <- function(sol){
  sol$beta
}

## Extract mixture, need to define mix as a S3 generic function
mix <- function(sol){
  sol$mix
}


## NOTE (Michael): Sometimes the gradient are not zero because the
##                 likelihood does not change anymore. This function
##                 checks whether the log-likelihood is maximised by
##                 altering the coefficient and check the
##                 log-likelihood around its neighbourhood.
gradCheck <- function(x, sol, tol = 1e-15){
  for(i in 1:length(sol$beta)){
    null <- rep(0, length(sol$beta))
    null[i] <- tol
    max.est <- all(logLik.snpmle(as.mgarch(x), sol$beta, sol$mix)[1] >
                   logLik.snpmle(as.mgarch(x), sol$beta + null, sol$mix)[1],
                   logLik.snpmle(as.mgarch(x), sol$beta, sol$mix)[1] >
                   logLik.snpmle(as.mgarch(x), sol$beta - null, sol$mix)[1])
    cat(paste("Estimation of ", names(sol$beta)[i], " is maximum: ",
              max.est, "\n", sep = ""))
  }
}

## This function constructs the condition standard deviation of a
## garch model.
##
## NOTES (Michael): Don't like the name of this function, will change it
##                  later.
##
## cond.sd <- function(x, beta){
##   T <- length(x)
##   betaSum <- as.numeric(filter(x[-T]^2, beta[2], "recursive"))
##   sigma.t <-
##     sqrt(c(beta[4]^2,
##            beta[1] * ((1 - cumprod(rep.int(beta[2], T - 1))))/
##            (1 - beta[2]) +
##            cumprod(rep.int(beta[2], T - 1)) * beta[4]^2 +
##            beta[3] * betaSum))
##   sigma.t
## }


cond.sd <- function(x, beta){
  sigma1 <- attr(x, "sigma1")
  T <- length(x)
  betaSum <- as.numeric(filter(x[-T]^2, beta[2], "recursive"))
  sigma.t <-
    sqrt(c(sigma1^2,
           beta[1] * ((1 - cumprod(rep.int(beta[2], T - 1))))/
           (1 - beta[2]) +
           cumprod(rep.int(beta[2], T - 1)) * sigma1^2 +
           beta[3] * betaSum))
  sigma.t
}
