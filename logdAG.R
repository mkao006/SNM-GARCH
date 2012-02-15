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
## TODO (Michael): Still not sure why a few analytical derivatives are
##                 significantly different to the numerical gradient for
##                 the derivative of beta[1] and second order for theta.

## NOTES (Michael): I think the biggest problem is actually numerical
##                  stability of the logd function
## NOTES (Michael): It appears that the analytical gradient are correct
##                  since if we reduce the size of increments, the
##                  difference between the analytical gradient and the
##                  numerical gradients becomes zero.

## NOTES (Michael): The use of heaviside function is so that there are
##                  not undefined points as would have happen if we used
##                  ifelse or something similar. However, will need to
##                  use something else as the heaviside gives equal
##                  weight to the left and the right hand side which is
##                  appropriate consider our non-linear function.


logd.mgarch <- function(xt, beta, pt, which){
  ## Initialise variables and list
  xt <- as.numeric(xt)
  T <- length(xt)
  lpt <- length(pt)
  lb <- length(beta)
  dl <- vector("list", length = 3)
  names(dl) <- c("ld", "db", "dt")

  ## Calculate the conditional variance
  betaSum <- as.numeric(filter(xt[-T]^2, beta[2], "recursive"))
    sigma.t <-
        sqrt(c(beta[4]^2,
               beta[1] * ((1 - cumprod(rep.int(beta[2], T - 1))))/
               (1 - beta[2]) +
               cumprod(rep.int(beta[2], T - 1)) * beta[4]^2 +
               beta[3] * betaSum))
  ##     sigma.t2 <-
  ##       c(beta[4]^2,
  ##              beta[1] * ((1 - cumprod(rep.int(beta[2], T - 1))))/
  ##              (1 - beta[2]) +
  ##              cumprod(rep.int(beta[2], T - 1)) * beta[4]^2 +
  ##              beta[3] * betaSum)
  ## print(head(sigma.t2, 10))
  ## dl$sigma1 <- sigma.t
  
  ## sigma2 <- double(T)
  ## sigma2[1] <- beta[4]^2
  ## for(i in 2:T){
  ##   sigma2[i] <- beta[1] + beta[2] * sigma2[i - 1] + beta[3] * xt[i - 1]^2
  ## }
  ## sigma.t2 <- sqrt(sigma2)
  ## dl$sigma2 <- sigma.t2

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

    ## Piece wise Analytical derivative
    dldsigma <- -1/sigma.t +
      (outer(sigma.t^-1, pt^-1) * sd_xi * xt + mu_xi) *
        (outer(sigma.t^-2, pt^-1) * sd_xi * xt) *
          (beta[5]^-2 * Heaviside(z_xi) + beta[5]^2 * Heaviside(-z_xi))

    ## NOTES (Michael): Check again by solving with Heaviside
    ##                  function. Unfortunately this solution is wrong
    ##                  because the sign function may be a reasonable
    ##                  representation of the Heaviside function but the
    ##                  property is quite different when taking into
    ##                  account of the skewness parameter
    ## dldsigma <-
    ##   -1/sigma.t -
    ##     (beta[5]^(sign(-z_xi)) * (z_xi)) *
    ##       (beta[5]^(-(outer(sigma.t^-2, pt^-1) * sd_xi * xt) * Delta(-z_xi)) *
    ##        (z_xi) +
    ##      beta[5]^(sign(-z_xi)) *
    ##       (-(outer(sigma.t^-2, pt^-1) * sd_xi * xt))
    
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
    dsigmadsigma <-
      2 * beta[4] * beta[2]^(0:(T - 1))
      ## c(2 * beta[4], rep(0, (T - 1)))
      
    
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
    dl$db <- abind(dldparams, dldxi, along = 3)
  }
  if(which[3] == 1){
    dl$dt <- -1/rep(pt, each = T) +
      (outer(sigma.t^-1, pt^-1) * sd_xi * xt + mu_xi) *
        (outer(sigma.t^-1, pt^-2) * sd_xi * xt) *
         (beta[5]^-2 * Heaviside(z_xi) + beta[5]^2 * Heaviside(-z_xi))    
  }
  dl
}

## Compile the function to increase the speed
logd.mgarch <- cmpfun(logd.mgarch)

## Test whether the parameters are valid
## NOTES (Michael): Removed the stationarity restriction
valid.mgarch <- function(x, beta, mix){
    beta[1] > 0 &&
    beta[1] < max(x)^2 &&
    beta[2] >= 0 &&
    beta[2] < 1 &&
    beta[3] >= 0 &&
    beta[3] < 1 &&
    beta[4] > 0 &&
    beta[4] < max(x) &&
    beta[5] > 0
}


## Function for initialising the parameters
##
## NOTES (Michael): Sometimes the coefficients given by fGarch voilates
##                  the stationarity assumption so we scale the
##                  coefficients to satisfy the assumption. Should we
##                  remove this restriction as well?
##
## TODO (Michael): Write a catch function if the garchFit fails,
##                 sometimes the GARCH can not be initialised when the
##                 inverse of the hessian does not exists
##
## NOTES (Michael): Use the standard deviation of the data to initiate,
##                   this avoids the problem of having initial value too
##                   different to the data in which somehow causes
##                   differential in the derivatives

initial.mgarch <- function(x, beta = NULL, mix = NULL, kmax = NULL){
    if(is.null(beta)){
        gf <- garchFit(data = as.numeric(x), trace = FALSE,
                             include.mean = FALSE)
        cgf <- coef(gf)
        if((sum(cgf[2:3]) >= 1))
          cgf[2:3] = cgf[2:3]/sum(cgf[2:3]) - 1e-3
        ## beta <- c(cgf[1], cgf[3], cgf[2], gf@sigma.t[1], 1)
        beta <- c(cgf[1], cgf[3], cgf[2], gf@sigma.t[1], 1)
        names(beta) <- c("omega", "beta1", "alpha1", "sigma0", "xi")
        mix <- disc(1, 1)
        list(beta = beta, mix = mix)
    }
}

## Function to determine the grid
gridpoints.mgarch <- function(x, beta, grid){
  seq(0.005, 20, length = grid)
}

## The weights function, this is useful when the data is discrete and
## there are duplicated.
weights.mgarch <- function(x) 1

## Function to restrict the space of the support points
suppspace.mgarch <- function(x, beta, mix){
  c(0.005, 20)
}

## Function for converting different class of time series to mgarch
## TODO (Michael): Check John Chamber's book and all the time series
##                 classes in R
as.mgarch <- function(x){
  mgts <- as.numeric(data.matrix(x))
  class(mgts) <- "mgarch"
  mgts
}  


## summary result of the solution
summary.mgarch <- function(sol, digits = 4){
  cat(paste("Solution Standard Deviation: ",
              round(sqrt(sum(sol$mix$pr * sol$mix$pt^2)), 4), "\n", sep = ""))
  cat(paste("Persistence: ", round(sum(sol$beta[2:3]), 4), "\n", sep = ""))
  cat(paste("Maximum gradient equal to zero: ",
            max(abs(sol$grad)) <=
            eval(parse(text = paste("1e-", digits, sep = ""))),
            "\n", sep = ""))
  cat(paste("Log-likelihood: ", round(sol$ll, 4), "\n", sep = ""))
}

## Extract coefficients
coef.mgarch <- function(sol){
  round(sol$beta, 7)
}

## TODO (Michael): This takes too long and inflexible, take out the
##                 fGarch component, and just write a pure plotting
##                 function.

## Plot the fit

plot.mgarch <- function(xt, sol, myylim = c(0, 0.75), bin = 10){
  ## Calculate the fit by different distribution
  norm.fit <- garchFit(data = xt, cond.dist = "snorm",
                       include.mean = FALSE, trace = FALSE)
  ged.fit <- garchFit(data = xt, cond.dist = "sged",
                      include.mean = FALSE, trace = FALSE)
  t.fit <- garchFit(data = xt, cond.dist = "sstd",
                    include.mean = FALSE, trace = FALSE)

  par(mfrow = c(2, 2), mar = c(2.1, 4.1, 4.1, 1.1))
  ## Plot the standard normal fit
  hist(xt/norm.fit@sigma.t, freq = FALSE,
       breaks = length(xt)/bin, xlim = c(-4, 4), ylim = myylim,
       main = "Standard Normal")
  curve(dsnorm(x, 0, 1, xi = coef(norm.fit)["skew"]),
        add = TRUE, col = "blue", lwd = 3)
  lines(density(xt/norm.fit@sigma.t), col = "red", lwd = 3)
  legend("topleft", legend = c("Density", "Fitted"),
       col = c("red", "blue"), bty = "n", lty = 1, lwd = 3)
  box()

  ## Plot the t-distribution fit
  par(mar = c(2.1, 3.1, 4.1, 2.1))
  hist(xt/t.fit@sigma.t, freq = FALSE,
       breaks = length(xt)/bin, xlim = c(-4, 4),
       ylim = c(0, 0.6),
       main = paste("t (", round(coef(t.fit)[5], 2), ")", sep = ""))
  curve(dsstd(x, 0, 1, nu = coef(t.fit)["shape"], xi = coef(t.fit)["skew"]),
        add = TRUE, col = "blue", lwd = 3)
  lines(density(xt/t.fit@sigma.t), col = "red", lwd = 3)
  box()

  ## Plot the ged fit
  par(mar = c(5.1, 4.1, 1.1, 1.1))
  hist(xt/ged.fit@sigma.t, freq = FALSE,
       breaks = length(xt)/bin, xlim = c(-4, 4),
       ylim = myylim,
       main = paste("Generalised Error Distribution (",
         round(coef(ged.fit)["shape"], 2), ")", sep = ""),
       ylab = "", xlab = "")
  curve(dsged(x, 0, 1, nu = coef(ged.fit)["shape"], xi = coef(ged.fit)["skew"]),
        add = TRUE, col = "blue", lwd = 3)
  lines(density(xt/ged.fit@sigma.t), col = "red", lwd = 3)
  box()

  ## Plot the mixture fit
  mix.sd <- cond.sd(xt, sol$beta)
  hist(xt/mix.sd, freq = FALSE,
       breaks = length(xt)/bin, ylim = myylim, xlim = c(-4, 4), main = "",
     xlab = "Error Distribution")
  lines(density(xt/mix.sd), col = "red", lwd = 3)
  curve(dmsnorm(x, sd = 1, varmix = sol$mix, xi = sol$beta["xi"]),
        add = TRUE, col = "blue", lwd = 3)
  box()
  for(i in 1:length(sol$mix$pt)){
    curve(sol$mix$pr[i] * dsnorm(x, sd = sol$mix$pt[i], xi = sol$beta["xi"]),
          add = TRUE, col = "light blue", lwd = 2)
  }
}

## x <- seq(-4, 4, length = 1000)
## den <- double(length(x))
## for(i in 1:length(x)){
##   den[i] <- dg.mg$mix$pr %*% dsnorm(x[i], sd = dg.mg$mix$pt, xi = 0.865)
## }
## lines(x, den, col = "orange", lwd = 3, lty = 3)

## lines(x, dmsnorm(x, sd = 1, varmix = dg.mg$mix, xi = dg.mg$beta["xi"]),
##       col = "green")



dmsnorm <- function(x, mean = 0, sd = 1, varmix = disc(1, 1), xi){
  n = length(x)
  n.mix = length(varmix$pt)

  ## Account for different cases of sigma
  ## (1) If the length of sd == 1 then it's constant variance
  ## (2) If the length of sd == n then it's conditional variance
  ## Otherwise there is no interpretation/meaning using the recyclying rule
  if(length(sd) == 1){
    sd <- rep(sd, n)
  }else if(length(sd) != n){
    stop("length of the standard deviation is not correct")
  }

  ## Calculate the transformed moments
  sd_xi <- sqrt((1 - 2/pi) * (xi^2 + xi^-2) + (4/pi - 1))
  mu_xi <- sqrt(2/pi) * (xi - xi^-1)
  Sigma <- outer(sd^-1, varmix$pt^-1)
  z_xi <- Sigma * sd_xi * x + mu_xi

  ## Calculate the log-density then the weighted density
  ll <- 0.5 * log(2/pi) + 0.5 * log(sd_xi^2) - log(xi + xi^-1) -
                  log(Sigma) -
        0.5 * ((Sigma * sd_xi * x + mu_xi) *
               (xi^-1 * Heaviside(z_xi) + xi * Heaviside(-z_xi)))^2
  exp(ll) %*% varmix$pr
  ## ((2 * Sigma)/(xi + 1/xi) * dnorm(z_xi)) %*% varmix$pr
}

## NOTES (Michael): Don't like the name of this function, will change it
##                  later.
cond.sd <- function(x, beta){
  T <- length(x)
  betaSum <- as.numeric(filter(x[-T]^2, beta[2], "recursive"))
  sigma.t <-
    sqrt(c(beta[4]^2,
           beta[1] * ((1 - cumprod(rep.int(beta[2], T - 1))))/
           (1 - beta[2]) +
           cumprod(rep.int(beta[2], T - 1)) * beta[4]^2 +
           beta[3] * betaSum))
  sigma.t
}

## curve(dmsnorm(x, mean = 0, sd = 1, dg.mg$mix, xi = dg.mg$beta[5]), -3, 3,
##       ylim = c(0, 0.5), col = "blue")
## curve(dmsnorm(x, mean = 0, sd = 1, xi = 1), lty = 2, add = TRUE)
## curve(dnorm(x), add = TRUE, col = "red", lty = 2)
## curve(dmsnorm(x, mean = 0, sd = 1, xi = 2), lty = 3, add = TRUE,
##       col = "orange", lwd = 2)
## curve(dsnorm(x, xi = 2), lty = 3, col = "purple", add = TRUE)
## dmsnorm(1:5, mean = 0, sd = 1:5, dg.mg$mix, xi = dg.mg$beta[5])

######################################################################
## Modify the cnmms function for scaling and also adding the class to
## the returned object
######################################################################

cnmms <- function (x, init = NULL, maxit = 1000,
                   model = c("spmle", "npmle"), tol = 1e-5,
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
        ## NOTES (Michael): The scaling restriction has been removed.

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
          print(paste("Scaled by :", sc, sep = ""))
          cat(paste("Parameters violating the stationarity boundary?: ",
                    !(new.sc == sc), "\n", sep = ""))
          mix$pt <- mix$pt/sc
          beta[c(1, 3)] <- beta[c(1, 3)] * sc^2
          beta[4] <- beta[4] * sc          
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
    result <- list(mix = mix, beta = beta, num.iterations = i,
                   ll = attr(mix,"ll")[1], grad = grad,
                   convergence = convergence)
    attr(result, "class") <- "mgarch"
    result
}
