########################################################################
## This is a rewrite of the log.mgarch function in which the mean
## estimated seperately and the density is reparametised
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

## TODO (Michael): Still not sure why a few analytical derivatives are
##                 significantly differet to the numerical gradient for
##                 the derivative of beta[1] and second order for theta.

logd.mgarch <- function(xt, beta, pt, which){
  ## Initialise variables and list
  xt <- as.numeric(xt)
  T <- length(xt)
  lpt <- length(pt)
  lb <- length(beta)
  dl <- vector("list", length = 4)  
  names(dl) <- c("ld", "db1", "dt1", "dt2")

  ## Calculate the conditional variance
  betaSum <- as.numeric(filter(xt[-T]^2, beta[2], "recursive"))
    sigma.t <-
        sqrt(c(beta[4]^2,
               beta[1] * ((1 - cumprod(rep.int(beta[2], T - 1))))/
               (1 - beta[2]) +
               cumprod(rep.int(beta[2], T - 1)) * beta[4]^2 +
               beta[3] * betaSum))

  ## Calculate the transformed moments of the skewed distribution
  sd_xi <- sqrt((1 - 2/pi) * (beta[5]^2 + beta[5]^-2) + (4/pi - 1))
  mu_xi <- sqrt(2/pi) * (beta[5] - beta[5]^-1)
  z_xi <- outer(sigma.t^-1, pt^-1) * sd_xi * xt + mu_xi

  ## Calculate the density
  if(which[1] == 1){
    dl$ld <- 0.5 * log(2/pi) + 0.5 * log(sd_xi^2) - log(beta[5] + beta[5]^-1) -
      log(sigma.t) - log(rep(pt, each = T)) -
        0.5 * ((outer(sigma.t^-1, pt^-1) * sd_xi * xt + mu_xi) *
               (beta[5]^-1 * Heaviside(z_xi) + beta[5] * Heaviside(-z_xi)))^2
  }

  ## Calculate the derivatives
  if(which[2] == 1){
    ## Use product rule to compute the first order derivative of the
    ## conditional variance parameters (omega, alpha, beta)
     
    dldsigma <- -1/sigma.t +
      (outer(sigma.t^-1, pt^-1) * sd_xi * xt + mu_xi) *
        (outer(sigma.t^-2, pt^-1) * sd_xi * xt) *
          (beta[5]^-2 * Heaviside(z_xi) + beta[5]^2 * Heaviside(-z_xi))

    sig.vec <- 2 * sigma.t
    convFilter <- 1:(T - 2) * cumprod(c(1, rep(beta[2], T - 3)))
    myconvolve <- beta[3] * convolve(convFilter, rev(xt^2),
                                     type = "open")[1:(T - 2)]
    cp.beta <- c(0, 0, myconvolve)

    ## beta[1] - Mean of the conditional variance equation
    dsigmadalpha0 <-
      c(0, ((1 - cumprod(rep.int(beta[2], T - 1))))/(1 - beta[2]))
    
    ## beta[2] - Coefficient for lagged variance
    dsigmadalpha1 <-
      (beta[1] * cumsum(c(0, 0, 1:(T - 2) * beta[2]^(0:(T - 3))))) +
        c(0, (1:(T - 1)) * beta[2]^(0:(T - 2)) * beta[4]^2) + 
         cp.beta
    
    ## beta[3] - Coefficient for lagged sqaured observation
    dsigmadbeta1 <- c(0, betaSum)

    ## beta[4] - The initial variance
    dsigmadsigma <- 2 * beta[4] * beta[2]^(0:(T - 1))
    
    ## beta[5] - Skewness parameter
    dldxi <- ((2 - 4/pi) * (beta[5] - beta[5]^-3))/
      (2 * ((1 - 2/pi) * (beta[5]^2 + beta[5]^-2) + (4/pi - 1))) +
        (beta[5]^-2 - 1)/(beta[5] + beta[5]^-1) -
          ((outer(sigma.t^-1, pt^-1) * sd_xi * xt + mu_xi) *
            (beta[5]^-1 * Heaviside(z_xi) + beta[5] * Heaviside(-z_xi))) *
              ((outer(sigma.t^-1, pt^-1) * sd_xi^-1 * (1 - 2/pi) *
                (beta[5] - beta[5]^-3) * xt + sqrt(2/pi) * (1 + beta[5]^-2)) *
               (beta[5]^-1 * Heaviside(z_xi) + beta[5] * Heaviside(-z_xi)) +
               (outer(sigma.t^-1, pt^-1) * sd_xi * xt + mu_xi) *
               (Heaviside(-z_xi) - beta[5]^-2 * Heaviside(z_xi)))

    ## Combine everything into an array.
    ## TODO (Michael):Improve this part, it feels a little bit
    ## clumsy. Also try to avoid abind.
    dldparams <- array(c(dldsigma * dsigmadalpha0/sig.vec,
                         dldsigma * dsigmadalpha1/sig.vec,
                         dldsigma * dsigmadbeta1/sig.vec,
                         dldsigma * dsigmadsigma/sig.vec),
                           dim = c(T, lpt, lb - 1))
    dl$db1 <- abind(dldparams, dldxi, along = 3)
  }
  if(which[3] == 1){
    dl$dt <- -1/rep(pt, each = T) +
      (outer(sigma.t^-1, pt^-1) * sd_xi * xt + mu_xi) *
        (outer(sigma.t^-1, pt^-2) * sd_xi * xt) *
         (beta[5]^-2 * Heaviside(z_xi) + beta[5]^2 * Heaviside(-z_xi))    
  }
  ## ## TODO(Michael): Check why this derivative is different to the
  ## ## numerical solution
  ## ## NOTES (Michael): This is not required in the new cnmms algorithm
  ## if(which[4] == 1){
  ##   dl$dt2 <- 1/rep(pt^2, each = T) -
  ##     (outer(sigma.t^-1, pt^-3) * sd_xi * xt) *
  ##       (3 * outer(sigma.t^-1, pt^-1) * sd_xi * xt + 2 * mu_xi) *          
  ##         (beta[5]^-2 * Heaviside(z_xi) + beta[5]^2 * Heaviside(-z_xi))
  ## }
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
gridpoints.mgarch <- function(x, beta, grid){
  seq(1e-10, 20, length = grid)
}

## NOTES (Michael): What is this weights function?
weights.mgarch <- function(x){
  return(1)
}

## Function to restrict the space of the support points
suppspace.mgarch <- function(x, beta, mix){
  c(0, Inf)
}


######################################################################
## Modify the cnmms function just for scaling
######################################################################

cnmms <- function (x, init = NULL, maxit = 1000,
                   model = c("spmle", "npmle"), tol = 1e-15,
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

        ## Scale the error distribution and the beta
        sc <- sqrt(sum(r$mix$pr * r$mix$pt^2))
        print("Likelihood Before Scaling:")
        print.snpmle(verb, x, r$mix, r$beta, gradient)
        if(abs(sc - 1) > tol & sc < 5){
            new.sc <- ifelse(valid(x, c(r$beta[1] * sc,
                                        r$beta[2],
                                        r$beta[3] * sc,
                                        r$beta[4] * sc,
                                        r$beta[5])),
                             sc, (1/(r$beta[2] + r$beta[3]) - 1e-7))
            print(paste("Scaled by :", new.sc, sep = ""))
            r$mix$pt <- r$mix$pt/new.sc
            r$beta[c(1, 3)] <- r$beta[c(1, 3)] * new.sc^2
            r$beta[4] <- r$beta[4] * new.sc
            r$conv <- 4
        }
        print("Likelihood After Scale:")
        print.snpmle(verb, x, r$mix, r$beta, gradient)
        
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
    list(mix = mix, beta = beta, num.iterations = i, ll = attr(mix, 
        "ll")[1], grad = grad, convergence = convergence)
}
