########################################################################
## This is a rewrite of the log.mgarch function in which the mean
## estimated seperately and the density is reparametised
## (Numerical gradient version)
########################################################################

## Load library
library(nspmix)
library(compiler)
library(Rcpp)
library(rbenchmark)
library(abind)
library(fGarch)

## Parameters:
## param[1] = omega
## param[2] = beta
## param[3] = alpha
## param[4] = sigma0
## param[5] = xi

## TODO (Michael): Handle the class in a better way. Don't like coercing
##                 to numeric.
## TODO (Michael): Improve the efficiency of this code (e.g. remove outer)
## TODO (Michael): Using the unbiased estimate of sigma0 causes problem
##                 when carrying our the scaling. Since now sigma0 is a
##                 function of both beta[1], beta[2] and beta[3] and a
##                 constant scaling factor may not exist.

## NOTES (Michael): It appears that the analytical gradient are correct
##                  since if we reduce the size of increments, the
##                  difference between the analytical gradient and the
##                  numerical gradients becomes zero.

## NOTES (Michael): Use parallel computing to do dbeta if possible

dmdsnorm <- function(xt, beta, pt){
  T <- length(xt)
  lpt <- length(pt)
  lb <- length(beta)
  
  betaSum <- as.numeric(filter(xt[-T]^2, beta[2], "recursive"))
  sigma.t <-
    sqrt(c(beta[4]^2,
           beta[1] * ((1 - cumprod(rep.int(beta[2], T - 1))))/
           (1 - beta[2]) +
           cumprod(rep.int(beta[2], T - 1)) * beta[4]^2 +
           beta[3] * betaSum))

  sd_xi <- sqrt((1 - 2/pi) * (beta[5]^2 + beta[5]^-2) + (4/pi - 1))
  mu_xi <- sqrt(2/pi) * (beta[5] - beta[5]^-1)
  z_xi <- outer(sigma.t^-1, pt^-1) * sd_xi * xt + mu_xi

  ll <- 0.5 * log(2/pi) + 0.5 * log(sd_xi^2) - log(beta[5] + beta[5]^-1) -
    log(sigma.t) - log(rep(pt, each = T)) -
      0.5 * ((outer(sigma.t^-1, pt^-1) * sd_xi * xt + mu_xi) *
             (beta[5]^-1 * Heaviside(z_xi) + beta[5] * Heaviside(-z_xi)))^2
  ll
}



logd.mgarch <- function(xt, beta, pt, which){
  ## Initialise variables and list
  xt <- as.numeric(xt)
  T <- length(xt)
  lpt <- length(pt)
  lb <- length(beta)
  dl <- vector("list", length = 3)
  names(dl) <- c("ld", "db", "dt")
  incrmt <- 1e-12
  
  ## Calculate the density
  if(which[1] == 1){
    dl$ld <- dmdsnorm(xt = xt, beta = beta, pt = pt)
  }

  ## Calculate the derivatives

  if(which[2] == 1){
    dl$db <- array(0, dim = c(T, lpt, lb))
    index <- rep(0, lb)
    for(i in 1:lb){
      index0 <- index
      index0[i] <- 1
      beta1 <- beta + index0 * incrmt
      beta2 <- beta - index0 * incrmt
      lplus <- dmdsnorm(xt = xt, beta = beta1, pt = pt)
      lminus <- dmdsnorm(xt = xt, beta = beta2, pt = pt)
      dl$db[,, i] <- (lplus - lminus)/(2 * incrmt)
    }
  }
  
  if(which[3] == 1){
    pt1 = pt + incrmt
    pt2 = pt - incrmt
    dl$dt <-
      (dmdsnorm(xt = xt, beta = beta, pt = pt1) -
       dmdsnorm(xt = xt, beta = beta, pt = pt2))/(2 * incrmt)
  }
  dl
}

## Compile the function to increase the speed
logd.mgarch <- cmpfun(logd.mgarch)

## Test whether the parameters are valid
valid.mgarch <- function(x, beta, mix){
    beta[1] > 0 &&
    beta[2] >= 0 &&
    beta[3] >= 0 &&
    beta[4] > 0 &&
    beta[5] > 0 &&
    (beta[2] + beta[3]) < 1
}


## Function for initialising the parameters
## NOTES (Michael): Sometimes the coefficients given by fGarch voilates
## the stationarity assumption so we scale the coefficients to satisfy
## the assumption
initial.mgarch <- function(x, beta = NULL, mix = NULL, kmax = NULL){
    if(is.null(beta)){
        gf <- garchFit(data = as.numeric(x), trace = FALSE,
                             include.mean = FALSE)
        cgf <- coef(gf)
        if((sum(cgf[2:3]) >= 1))
          cgf[2:3] = cgf[2:3]/sum(cgf[2:3]) - 1e-3
        beta <- c(cgf[1], cgf[3], cgf[2], gf@sigma.t[1], 1)
        names(beta) <- c("omega", "beta1", "alpha1", "sigma0", "xi")
        mix <- disc(1, 1)
        list(beta = beta, mix = mix)
    }
}

## Function to determine the grid
## TODO (Michael): Write the min of the grid and support space function
##                 as a function of the data
gridpoints.mgarch <- function(x, beta, grid){
  seq(0.05, 20, length = grid)
}

## The weights function, this is useful when the data is discrete and
## there are duplicated.
weights.mgarch <- function(x) 1

## Function to restrict the space of the support points
## NOTES (Michael): The support space must be a subspace of the grid
##                  points
suppspace.mgarch <- function(x, beta, mix){
  c(0.05, 20)
}

## Function for converting different class of time series to mgarch
## TODO (Michael): Check John Chamber's book and survey all the time
##                 series classes in R
as.mgarch <- function(x){
  mgts <- as.numeric(data.matrix(x))
  class(mgts) <- "mgarch"
  mgts
}

## summary result of the solution
summary.mgarch <- function(sol, tol = 1e-5){
  cat(paste("Solution Standard Deviation: ",
              sqrt(sum(sol$mix$pr * sol$mix$pt^2)), "\n", sep = ""))
  cat(paste("Persistence :", sum(sol$beta[2:3])), "\n", sep = "")
  cat(paste("Maximum gradient equal to zero :", max(abs(sol$grad)) == 0,
            "\n", sep = ""))
}



######################################################################
## Modify the cnmms function just for scaling
######################################################################

cnmms <- function (x, init = NULL, maxit = 1000,
                   model = c("spmle", "npmle"), tol = 1e-6,
                   grid = 100, kmax = Inf,
                   plot = c("null", "gradient", "prob"), verb = 0){
    plot = match.arg(plot)
    model = match.arg(model)
    if (kmax == Inf) 
        init = initial.snpmle(x, init)
    else init = initial.snpmle(x, init, kmax = kmax)
    beta = init$beta
    if (is.null(beta) || is.na(beta)) 
        model = "npmle"
    nb = length(beta)
    mix = init$mix
    gradient = "Not computed"
    switch(plot, gradient = plotgrad(x, beta, mix, pch = 1), 
        prob = plot(x, mix, beta))
    ll1 = -Inf
    convergence = 1
    w = weights(x)
    wr = sqrt(w)
    for (i in 1:maxit) {
        l = logd(x, beta, mix$pt, which = c(1, 0, 0))$ld
        ma = apply(l, 1, max)
        dmix = drop(exp(l - ma) %*% mix$pr) + 1e-100
        if (length(mix$pt) < kmax) {
            gp = gridpoints(x, beta, grid)
            g = maxgrad(x, beta, dmix, ma, grid = gp, tol = -Inf)
            gradient = max(g$grad)
            kpt = min(kmax - length(mix$pt), length(g$pt))
            jpt = order(g$grad, decreasing = TRUE)
            mix = disc(c(mix$pt, g$pt[jpt][1:kpt]), c(mix$pr, 
                rep(0, kpt)))
        }
        lpt = logd(x, beta, mix$pt, which = c(1, 0, 0))$ld
        dpt = pmin(exp(lpt - ma), 1e+100)
        a = wr * (dpt/dmix - 2)
        grad.support = colSums(w * (dpt/dmix - 1))
        r = nnls(rbind(a, rep(1, length(mix$pt))), c(rep(0, nrow(a)), 
            1))
        sol = r$x/sum(r$x)
        r = lsch(mix, beta, disc(mix$pt, sol), beta, x, which = c(1, 
            0, 0))
        mix = collapse.snpmle(r$mix, beta, x)
        if (max(grad.support) < 1e+05) {
            r = switch(model, spmle = bfgs(mix, beta, x, which = c(1, 
                1, 1)), npmle = bfgs(mix, beta, x, which = c(1, 
                1, 0)))
            beta = r$beta
            mix = collapse.snpmle(r$mix, beta, x)
        }
        switch(plot, gradient = plotgrad(x, beta, mix, pch = 1), 
            prob = plot(x, mix, beta))

        ## NOTES (Michael): There is a problem here, some times the
        ##                  algorithm just breaks out here before the
        ##                  last step of scaling so that the final
        ##                  solution is actually not scaled.

        ## Start of scaling
        sc <- sqrt(sum(mix$pr * mix$pt^2))
        print("Likelihood Before Scaling:")
        print.snpmle(verb, x, mix, beta, gradient)
        if(abs(sc - 1) > tol & sc < 5){
          new.sc <- ifelse(valid(x, c(beta[1] * sc,
                                      beta[2],
                                      beta[3] * sc,
                                      beta[4] * sc,
                                      beta[5])),
                           sc, (1/(beta[2] + beta[3]) - 1e-7))
          print(paste("Scaled by :", new.sc, sep = ""))
          cat(paste("Close to boundary?: ", !(new.sc == sc), "\n", sep = ""))
          mix$pt <- mix$pt/new.sc
          beta[c(1, 3)] <- beta[c(1, 3)] * new.sc^2
          beta[4] <- beta[4] * new.sc
          
        }
        print("Likelihood After Scale:")
        print.snpmle(verb, x, mix, beta, gradient)

        print(beta)
        ## End of scaling
        
        if (r$ll >= ll1 - tol * abs(r$ll) && r$ll <= ll1 + tol * 
            abs(r$ll)) {
            convergence = 0
            break
        }
        ll1 = r$ll
        print.snpmle(verb, x, mix, beta, gradient)
    }
    m = length(mix$pt)
    if (m < length(r$mix$pt)) {
        d = dll.snpmle(x, mix, beta, which = c(0, 1, 1, 1))
        grad = c(d$dp, d$dt, d$db)
        names(grad) = c(paste("pr", 1:k, sep = "."), paste("pt", 
            1:k, sep = "."), paste("beta", 1:length(beta), sep = "."))
    }
    else grad = r$grad
    grad[1:m] = grad[1:m] - sum(rep(w, len = length(x)))
    list(mix = mix, beta = beta, num.iterations = i,
         ll = attr(mix,"ll")[1], grad = grad, convergence = convergence)
}


