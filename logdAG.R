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
  dl <- vector("list", length = 4)
  names(dl) <- c("ld", "db", "dt", "sigma.t")

  ## Calculate the conditional variance
  betaSum <- as.numeric(filter(xt[-T]^2, beta[2], "recursive"))
    sigma.t <-
        sqrt(c(beta[4]^2,
               beta[1] * ((1 - cumprod(rep.int(beta[2], T - 1))))/
               (1 - beta[2]) +
               cumprod(rep.int(beta[2], T - 1)) * beta[4]^2 +
               beta[3] * betaSum))
##   dl$sigma.t <- sigma.t
  
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
valid.mgarch <- function(x, beta){
    beta[1] > 0 &&
    beta[2] >= 0 &&
    beta[3] >= 0 &&
    beta[4] > 0 &&
    beta[5] > 0 
##    mix$pt >= 0.1
}

## valid.snpmle <- function(x, beta, mix){
##   bs = suppspace(x)
##   valid.mgarch(x, beta, mix) && all(mix$pr >= 0, mix$pt >= bs[1], mix$pt <= 
##         bs[2])
## }


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
        gf <- try(garchFit(data = as.numeric(x), trace = FALSE,
                           include.mean = FALSE))
        if(!(inherits(gf, "try-error"))){
          cgf <- coef(gf)
        } else {
          cgf <- c(1e-6, 0.1, 0.8, sd(x), 1)
        }
        ## if((sum(cgf[2:3]) >= 1))
        ##   cgf[2:3] = cgf[2:3]/sum(cgf[2:3]) - 1e-3
        ## beta <- c(cgf[1], cgf[3], cgf[2], gf@sigma.t[1], 1)
        beta <- c(cgf[1], cgf[3], cgf[2], gf@sigma.t[1], 1)
        names(beta) <- c("omega", "beta1", "alpha1", "sigma0", "xi")
        mix <- disc(1, 1)
        list(beta = beta, mix = mix)
    }
}

## Function to determine the grid
gridpoints.mgarch <- function(x, beta, grid){
  seq(0.1, 20, length = grid)
}

## The weights function, this is useful when the data is discrete and
## there are duplicated.
weights.mgarch <- function(x) 1

## Function to restrict the space of the support points
suppspace.mgarch <- function(x, beta, mix){
  c(0.1, 20)
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

## Error function taken from the VGAM package
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

## Plot the fit, also return the Komogorov-smirnov statistic
plot.mgarch <- function(xt, sol, plot = TRUE){
  T <- length(xt)
  betaSum <- as.numeric(filter(xt[-T]^2, sol$beta[2], "recursive"))
  sigma.t <-
    sqrt(c(sol$beta[4],
           sol$beta[1] * ((1 - cumprod(rep.int(sol$beta[2], T - 1))))/
           (1 - sol$beta[2]) +
           cumprod(rep.int(sol$beta[2], T - 1)) * sol$beta[4]^2 +
           sol$beta[3] * betaSum))
  res <- sort(as.numeric(xt/sigma.t))
  res.ecdf <- ecdf(res)
  if(plot){
    plot(res.ecdf,
         main = "Empirical CDF of the innovation\n with fitted mixture",
         xlim = c(-6, 6), lwd = 2)
    curve(pmsnorm(x, varmix = sol$mix, xi = sol$beta["xi"]), add = TRUE,
          col = "red", lwd = 1.5)
    curve(pnorm(x), add = TRUE, col = "blue")
  }
    ks <- 1:length(res)/length(res) -
      pmsnorm(res, varmix = sol$mix, xi = sol$beta["xi"])
  list(res = res, ks = ks)
}


disStudy <- function(n.samp = 100, n.iter = 50,
                     param = list(omega = 1e-6, alpha = 0.15, beta = 0.75,
                       skew = 1),
                     cond.dist = "snorm", mix = disc(c(1, 1)), plot = FALSE){
  ## Initialisation
  n <- n.samp
  iter <- n.iter
  norm.res <- matrix(0, nc = iter, nr = n)
  norm.ks <- matrix(0, nc = iter, nr = n)
  t.res <- matrix(0, nc = iter, nr = n)
  t.ks <- matrix(0, nc = iter, nr = n)
  ged.res <- matrix(0, nc = iter, nr = n)
  ged.ks <- matrix(0, nc = iter, nr = n)
  mnorm.res <- matrix(0, nc = iter, nr = n)
  mnorm.ks <- matrix(0, nc = iter, nr = n)
  for(i in 1:iter){

    cat("\t\t\t--------------\n")
    cat(paste("\t\t\tSimulation: ", i, "\n", sep = ""))
    cat("\t\t\t--------------\n")
    
    if(!(cond.dist %in% c("msnorm", "stable"))){
      myspec <- garchSpec(model = param, cond.dist = cond.dist, rseed = i)
      n.sim <- garchSim(spec = myspec, n = n)
    } else {
      msnorm.beta <- c(param$omega, param$beta, param$alpha, 0.01, param$skew)
      n.sim <- mgarchSim(n = n, beta = msnorm.beta, mix = mix, seed = i,
                         cond.dist = cond.dist)$x
    }

    norm.fit <- try(garchFit(data = n.sim, trace = FALSE, cond.dist = "snorm",
                           include.mean = FALSE))
    if(inherits(norm.fit, "try-error")){
      norm.res[, i] <- rep(NA, n)
      norm.ks[, i] <- rep(NA, n)
    } else {
      norm.res[, i] <- sort(as.numeric(n.sim/norm.fit@sigma.t))
      norm.ks[, i] <- abs(1:n/n - psnorm(norm.res[, i],
                                         xi = coef(norm.fit)["skew"]))
    }

    t.fit <- try(garchFit(data = n.sim, trace = FALSE, cond.dist = "sstd",
                          include.mean = FALSE))
    if(inherits(t.fit, "try-error")){
      t.res[, i] <- rep(NA, n)
      t.ks[, i] <- rep(NA, n)
    } else {
      t.res[, i] <- sort(as.numeric(n.sim/t.fit@sigma.t))
      t.ks[, i] <- abs(1:n/n - psstd(t.res[, i], xi = coef(t.fit)["skew"],
                                     nu = coef(t.fit)["shape"]))
    }

    ged.fit <- try(garchFit(data = n.sim, trace = FALSE, cond.dist = "sged",
                            include.mean = FALSE))
    if(inherits(ged.fit, "try-error")){
      ged.res[, i] <- rep(NA, n)
      ged.ks[, i] <- rep(NA, n)
    } else {
      ged.res[, i] <- sort(as.numeric(n.sim/ged.fit@sigma.t))
      ged.ks[, i] <- abs(1:n/n - psged(ged.res[, i], xi = coef(ged.fit)["skew"],
                                       nu = coef(ged.fit)["shape"]))
    }
    mnorm.fit <- try(cnmms(as.mgarch(n.sim), plot = "null", grid = 1000,
                           verb = 0, tol = 1e-5))
    if(inherits(mnorm.fit, "try-error")){
      mnorm.res[, i] <- rep(NA, n)
      mnorm.ks[, i] <- rep(NA, n)
    } else {
      tmp <- plot(xt = as.mgarch(n.sim), sol = mnorm.fit, plot = plot)
      mnorm.res[, i] <- tmp$res
      mnorm.ks[, i] <- abs(tmp$ks)
    }
  }

  plot.new()
  plot.window(xlim = c(1, n),
              ylim = c(0, max(cbind(norm.ks, t.ks, ged.ks, mnorm.ks),
                na.rm = TRUE)))
  for(i in 1:iter){
    lines(norm.ks[, i])
    lines(t.ks[, i], col = "green")
    lines(ged.ks[, i], col = "red")
    lines(mnorm.ks[, i], col = "orange")
  }
  legend("topleft", col = c("black", "green", "red", "orange"), lty = 1,
         lwd = 3, legend = c("Normal", "t", "ged", "mnorm"), bty = "n")
  text(n * 0.8,
       c(max(cbind(norm.ks, t.ks, ged.ks, mnorm.ks), na.rm = TRUE) *
         rev(seq(8, 9.5, by = 0.5)/10)),
       labels = c("norm", "t", "ged", "mnorm"))
  text(n,
       c(max(cbind(norm.ks, t.ks, ged.ks, mnorm.ks), na.rm = TRUE) *
         rev(seq(8, 10, by = 0.5)/10)),
       labels = c("max",
         round(c(max(norm.ks, na.rm = TRUE), max(t.ks, na.rm = TRUE),
               max(ged.ks, na.rm = TRUE), max(mnorm.ks, na.rm = TRUE)), 8)))
  text(n * 0.9,
       c(max(cbind(norm.ks, t.ks, ged.ks, mnorm.ks), na.rm = TRUE) *
         rev(seq(8, 10, by = 0.5)/10)),
       labels = c("mean",
         round(c(mean(norm.ks, na.rm = TRUE), mean(t.ks, na.rm = TRUE),
                 mean(ged.ks, na.rm = TRUE), mean(mnorm.ks, na.rm = TRUE)), 8)))
  axis(1)
  axis(2)
  box()
  lines(rowMeans(norm.ks, na.rm = TRUE), lwd = 5)
  lines(rowMeans(t.ks, na.rm = TRUE), col = "green", lwd = 5)
  lines(rowMeans(ged.ks, na.rm = TRUE), col = "red", lwd = 5)
  lines(rowMeans(mnorm.ks, na.rm = TRUE), col = "orange", lwd = 5)
  print(sum(colSums(norm.res, na.rm = TRUE) == 0))
  print(sum(colSums(t.res, na.rm = TRUE) == 0))
  print(sum(colSums(ged.res, na.rm = TRUE) == 0))
  print(sum(colSums(mnorm.res, na.rm = TRUE) == 0))
  list(norm.res, norm.ks, t.res, t.ks, ged.res, ged.ks, mnorm.res, mnorm.ks)
}


prConstruct <- function(pr, pt) sqrt((1 - pt^2 * pr)/(1 - pr))


## ## This is old and probably not ideal (Takes too long to run)
## plot.mgarch <- function(xt, sol, myylim = c(0, 0.75), bin = 10){
##   ## Calculate the fit by different distribution
##   norm.fit <- garchFit(data = xt, cond.dist = "snorm",
##                        include.mean = FALSE, trace = FALSE)
##   ged.fit <- garchFit(data = xt, cond.dist = "sged",
##                       include.mean = FALSE, trace = FALSE)
##   t.fit <- garchFit(data = xt, cond.dist = "sstd",
##                     include.mean = FALSE, trace = FALSE)
##
##   par(mfrow = c(2, 2), mar = c(2.1, 4.1, 4.1, 1.1))
##   ## Plot the standard normal fit
##   hist(xt/norm.fit@sigma.t, freq = FALSE,
##        breaks = length(xt)/bin, xlim = c(-4, 4), ylim = myylim,
##        main = "Standard Normal")
##   curve(dsnorm(x, 0, 1, xi = coef(norm.fit)["skew"]),
##         add = TRUE, col = "blue", lwd = 3)
##   lines(density(xt/norm.fit@sigma.t), col = "red", lwd = 3)
##   legend("topleft", legend = c("Density", "Fitted"),
##        col = c("red", "blue"), bty = "n", lty = 1, lwd = 3)
##   box()
##
##   ## Plot the t-distribution fit
##   par(mar = c(2.1, 3.1, 4.1, 2.1))
##   hist(xt/t.fit@sigma.t, freq = FALSE,
##        breaks = length(xt)/bin, xlim = c(-4, 4),
##        ylim = c(0, 0.6),
##        main = paste("t (", round(coef(t.fit)[5], 2), ")", sep = ""))
##   curve(dsstd(x, 0, 1, nu = coef(t.fit)["shape"], xi = coef(t.fit)["skew"]),
##         add = TRUE, col = "blue", lwd = 3)
##   lines(density(xt/t.fit@sigma.t), col = "red", lwd = 3)
##   box()
##
##   ## Plot the ged fit
##   par(mar = c(5.1, 4.1, 1.1, 1.1))
##   hist(xt/ged.fit@sigma.t, freq = FALSE,
##        breaks = length(xt)/bin, xlim = c(-4, 4),
##        ylim = myylim,
##        main = paste("Generalised Error Distribution (",
##          round(coef(ged.fit)["shape"], 2), ")", sep = ""),
##        ylab = "", xlab = "")
##   curve(dsged(x, 0, 1, nu = coef(ged.fit)["shape"], xi = coef(ged.fit)["skew"]),
##         add = TRUE, col = "blue", lwd = 3)
##   lines(density(xt/ged.fit@sigma.t), col = "red", lwd = 3)
##   box()
##
##   ## Plot the mixture fit
##   mix.sd <- cond.sd(xt, sol$beta)
##   hist(xt/mix.sd, freq = FALSE,
##        breaks = length(xt)/bin, ylim = myylim, xlim = c(-4, 4), main = "",
##      xlab = "Error Distribution")
##   lines(density(xt/mix.sd), col = "red", lwd = 3)
##   curve(dmsnorm(x, varmix = sol$mix, xi = sol$beta["xi"]),
##         add = TRUE, col = "blue", lwd = 3)
##   box()
##   for(i in 1:length(sol$mix$pt)){
##     curve(sol$mix$pr[i] * dsnorm(x, sd = sol$mix$pt[i], xi = sol$beta["xi"]),
##           add = TRUE, col = "light blue", lwd = 2)
##   }
## }
##
## x <- seq(-4, 4, length = 1000)
## den <- double(length(x))
## for(i in 1:length(x)){
##   den[i] <- dg.mg$mix$pr %*% dsnorm(x[i], sd = dg.mg$mix$pt, xi = 0.865)
## }
## lines(x, den, col = "orange", lwd = 3, lty = 3)
##
## lines(x, dmsnorm(x, sd = 1, varmix = dg.mg$mix, xi = dg.mg$beta["xi"]),
##       col = "green")


dmsnorm <- function(x, mean = 0, varmix = disc(1, 1), xi = 1){
  ## Function for calculating the PDF of the scaled skewed normal
  ## mixture distribution.
  ##
  ## Args:
  ##  x:      The observations
  ##  mean:   Mean of the mixture (can only be a scalar)
  ##  varmix: The variance mixture of the distribution
  ##  xi:     The skewness of the mixture distribution
  ##
  ## Returns:
  ##  The density of the scale skewed normal mixture distribution
  ##

  n = length(x)
  n.mix = length(varmix$pt)
  sx = (x - mean)/matrix(rep(varmix$pt, each = n), nc = n.mix)

  ## Calculate the transformed moments
  sd_xi <- sqrt((1 - 2/pi) * (xi^2 + xi^-2) + (4/pi - 1))
  mu_xi <- sqrt(2/pi) * (xi - xi^-1)
  z_xi <- sx * sd_xi + mu_xi
  Xi = xi^sign(z_xi)
  g = 2/(xi + 1/xi)
  den = g * dnorm(x = z_xi/Xi)
  Density = (den * sd_xi)/matrix(rep(varmix$pt, each = n), nc = n.mix)  
  Density %*% varmix$pr
}


pmsnorm <- function(q, mean = 0, varmix = disc(1, 1), xi = 1){
  ## Function for calculating the CDF of the mixture skewed normal
  ## distribution.
  ##
  ## Args:
  ##  x:      The observations
  ##  mean:   Mean of the mixture (can only be a scalar)
  ##  sd:     The standard deviation over time (only for time series)
  ##          otherwise treated as constant 1.
  ##  varmix: The mixture of the distribution
  ##  xi:     The skewness of the mixture distribution
  ##
  ## Returns:
  ##  The CDF of the scale mixture skewed normal distribution.
  ##
  
  n = length(q)
  n.mix = length(varmix$pt)
  sq = (q - mean)/matrix(rep(varmix$pt, each = n), nc = n.mix)
  

  ## Calculate the transformed moments
  sd_xi <- sqrt((1 - 2/pi) * (xi^2 + xi^-2) + (4/pi - 1))
  mu_xi <- sqrt(2/pi) * (xi - xi^-1)
  z_xi <- sq * sd_xi  + mu_xi
  Xi = xi^sign(z_xi)
  g = 2/(xi + 1/xi)
  prob = Heaviside(z_xi) - sign(z_xi) * g * Xi * pnorm(q = -abs(z_xi)/Xi)
  prob %*% varmix$pr
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
        if(verb != 0){
          print("Likelihood Before Scaling:")
          print.snpmle(verb, x, mix, beta, gradient)
        }
        if(abs(sc - 1) > tol & sc < 5){
          new.sc <- ifelse(valid.mgarch(x, c(beta[1] * sc,
                                      beta[2],
                                      beta[3] * sc,
                                      beta[4] * sc,
                                      beta[5])),
                           sc, (1/(beta[2] + beta[3]) - 1e-7))
          if(verb != 0){
            print(paste("Scaled by :", sc, sep = ""))
            cat(paste("Parameters violating the stationarity boundary?: ",
                      !(new.sc == sc), "\n", sep = ""))
          }
          mix$pt <- mix$pt/sc
          beta[c(1, 3)] <- beta[c(1, 3)] * sc^2
          beta[4] <- beta[4] * sc          
        }
        if(verb != 0){
          print("Likelihood After Scale:")
          print.snpmle(verb, x, mix, beta, gradient)
          print(beta)
        }

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
    ## sigma.t <- logd.mgarch(x, beta, mix$pt, which = c(0, 0, 0))$sigma.t
    ## result <- list(mix = mix, beta = beta, num.iterations = i,
    ##                ll = attr(mix,"ll")[1], grad = grad,
    ##                convergence = convergence,
    ##                sigma.t = sigma.t)
    result <- list(mix = mix, beta = beta, num.iterations = i,
                   ll = attr(mix,"ll")[1], grad = grad,
                   convergence = convergence)
    attr(result, "class") <- "mgarch"
    result
}



## NOTES (Michael):The current number of iteration in line search is
## insufficient in our case, thus I have increased the number of
## iterations so the algorithm does not break pre-maturely
lsch <- function (mix1, beta1, mix2, beta2, x, maxit = 100,
                  which = c(1, 1, 1), brkt = FALSE){
    k = length(mix1$pt)
    convergence = 1
    dl1 = dll.snpmle(x, mix1, beta1, which = c(1, which))
    lla = ll1 = dl1$ll
    names.grad = c(if (which[1]) paste("pr", 1:k, sep = ".") else NULL, 
        if (which[2]) paste("pt", 1:k, sep = ".") else NULL, 
        if (which[3]) paste("beta", 1:length(beta1), sep = ".") else NULL)
    grad1 = c(if (which[1]) dl1$dp else NULL, if (which[2]) dl1$dt else NULL, 
        if (which[3]) dl1$db else NULL)
    names(grad1) = names.grad
    d1 = c(if (which[1]) mix2$pr - mix1$pr else NULL, if (which[2]) mix2$pt - 
        mix1$pt else NULL, if (which[3]) beta2 - beta1 else NULL)
    d1.norm = sqrt(sum(d1 * d1))
    s = d1/d1.norm
    g1d1 = sum(grad1 * d1)
    dla = g1s = g1d1/d1.norm
    if (d1.norm == 0 || g1s <= 0) {
        return(list(mix = mix1, beta = beta1, grad = grad1, ll = ll1, 
            convergence = 3))
    }
    a = 0
    b = 1
    if (which[1] && any(mix2$pr == 0)) 
        brkt = FALSE
    for (i in 1:maxit) {
        for (j in 1:1000) {
            m = disc((1 - b) * mix1$pt + b * mix2$pt, (1 - b) * 
                mix1$pr + b * mix2$pr)
            beta = if (is.null(beta1)) 
                NULL
            else (1 - b) * beta1 + b * beta2
            if (valid.snpmle(x, beta, m)) 
                break
            brkt = FALSE
            b = 0.5 * a + 0.5 * b
        }
        if (j == 1000) 
            warning("Can not produce valid interior point in lsch()")
        dl = dll.snpmle(x, m, beta, which = c(1, which))
        ll = dl$ll
        grad = c(if (which[1]) dl$dp else NULL, if (which[2]) dl$dt else NULL, 
            if (which[3]) dl$db else NULL)
        gs = sum(grad * s)
        if (brkt && gs > g1s * 0.5 && ll >= ll1 + g1d1 * b * 
            0.33) {
            a = b
            b = 2 * b
            lla = ll
            dla = gs
        }
        else break
    }
    if (i == maxit) 
        brkt = FALSE
    alpha = b
    llb = ll
    dlb = gs
    for (i in 1:maxit) {
        g1d = g1d1 * alpha
        if (ll >= ll1 - 1e-15 * abs(ll1) && g1d <= 1e-15 * abs(ll)) {
            convergence = 2
            break
        }
        if (brkt) {
            if (ll >= ll1 + g1d * 0.33 && abs(gs) <= g1s * 0.5) {
                convergence = 0
                break
            }
            if (ll >= ll1 + g1d * 0.33 && gs > g1s * 0.5) {
                a = alpha
                lla = ll
                dla = gs
            }
            else {
                b = alpha
                llb = ll
                dlb = gs
            }
        }
        else {
            if (ll >= ll1 + g1d * 0.33) {
                convergence = 0
                break
            }
            else {
                b = alpha
                llb = ll
                dlb = gs
            }
        }
        alpha = (a + b) * 0.5
        m = disc((1 - alpha) * mix1$pt + alpha * mix2$pt, (1 - 
            alpha) * mix1$pr + alpha * mix2$pr)
        beta = if (is.null(beta1)) 
            NULL
        else (1 - alpha) * beta1 + alpha * beta2
        dl = dll.snpmle(x, m, beta, which = c(1, which))
        ll = dl$ll
        grad = c(if (which[1]) dl$dp else NULL, if (which[2]) dl$dt else NULL, 
            if (which[3]) dl$db else NULL)
        gs = sum(grad * s)
    }
    names(grad) = names.grad
    beta = if (is.null(beta1)) 
        NULL
    else (1 - alpha) * beta1 + alpha * beta2
    list(mix = disc((1 - alpha) * mix1$pt + alpha * mix2$pt, 
        (1 - alpha) * mix1$pr + alpha * mix2$pr), beta = beta, 
        grad = grad, ll = ll, convergence = convergence, num.iterations = i)
}


sd.disc <- function(dis){
  sqrt(sum(dis$pt^2 * dis$pr))
}


## NOTES (Michael): Now the only error left is the following:
##
## Error in if (brkt && gs > g1s * 0.5 && ll >= ll1 + g1d1 * b * 0.33) {
## : missing value where TRUE/FALSE needed


## setGeneric("valid", function(x, beta, mix) standardGeneric("valid"))

