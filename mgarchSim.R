######################################################################
## Functions to support the simulation study
######################################################################

## TODO (Michael): Rewrite this section so that you can generate a Garch
##                 process with different distributions


## Random number generation for the skewed scale normal mixture
## distribution
rssnormm <- function(n, mix = disc(1, 1), xi = 1){
    ## A function to generate samples from a scale normal mixture
    ## distribution.
    ##
    ## Args:
    ##   n:   The number of observation to be generated as used in other
    ##        "rdist" functions
    ##   mix: The mixture function generated using the disc function
    ##   xi: The skewness parameter of the distribution
    ##
    ## Output:
    ##   m:   The index of which component distribution the random
    ##        number was drawn from. (e.g. 1 equal the first component
    ##        density etc...)
    ##   sim: The list of random numbers generated.

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
    sim <- fGarch::rsnorm(n, mean = 0, sd = mix$pt[m], xi = xi)
    list(m = m, sim = sim)
}

## check <- rssnormm(1000, mix = disc(c(0, 5, 10), c(0.5, 0.25, 0.25)), xi = 1)
## table(check$m)
## hist(check$sim, breaks = 50)


mgarchSim <- function(n, beta, mix = NA, cond.dist = "msnorm", seed = NULL){
    ## Function to generate GARCH time series with scale normal
    ## mixture distribution
    ##
    ## Args:
    ##   beta: [1] omega: The mean of the conditional variance equation.
    ##         [2] beta1: The coefficient for the variance in the
    ##                    conditional variance equation
    ##         [3] alpha1: The coefficient of the squared observation
    ##                     in the conditional variance equation.
    ##         [4] sigma_0: The estimated initial value of the
    ##                      conditional variance
    ##         [5] xi: The skewness parameter of the distribution.
    ##   n: The length of the simulated series to be generated.
    ##   mix: The mixing distribution to be used.
    ##   seed: The initial seed for the random number.
    ##
    ## Output:
    ##   x: The GARCH time series with the desired mixing error
    ##      distribution.
    ##   sigma.t: The underlying volatility that was used to generate
    ##            the time series.

    if(!is.null(seed)) set.seed(seed)

    ## Initialise the vector and their first observations.
    sim <- double(n)
    sigma2 <- double(n)
    
    sigma2[1] <- beta[4]^2
    if(cond.dist == "msnorm"){
      rand <- rssnormm(1, mix, beta[5])$sim
    } else if(cond.dist == "stable"){
      rand <- rstable(1, alpha = 1.6, beta = 0, gamma = 0.7)
    }
    sim[1] <- sqrt(sigma2[1]) * rand

    ## Loop to generate the underlying volatility first and then the
    ## observations.

    for(i in 2:n){
        sigma2[i] = beta[1] + beta[2] * sigma2[i - 1] +
            beta[3] * sim[i - 1]^2
        if(cond.dist == "msnorm"){
          error <- rssnormm(1, mix, beta[5])$sim
        } else if(cond.dist == "stable"){
          error <- rstable(1, alpha = 1.6, beta = 0, gamma = 0.7)
        }
        sim[i] = sqrt(sigma2[i]) * error
    }
    ## class(sim) <- "mgarch"
    list(x = sim, sigma2 = sigma2)
}

mysim <- mgarchSim(beta = c(1e-3, 0.8, 0.1, 0.001, 2), n = 1000,
                   cond.dist = "msnorm",
                   mix = disc(c(1, 1.2, 1.5), rep(1/3, 3)))
plot(mysim$x, type = "l")


mysim <- mgarchSim(beta = c(1e-3, 0.8, 0.1, 0.001, 2), n = 1000,
                   cond.dist = "stable")
plot(mysim$x, type = "l")


## Testing codes for stable distribution
## Alpha
pSpace <- seq(1e-3, 2, length = 20)
curve(dnorm(x), -5, 5, col = "red", lwd = 2)
for(i in 1:length(pSpace)){
  curve(dstable(x, alpha = pSpace[i], beta = 0, pm = 0), add = TRUE,
        col = rgb(0, 0, i/length(pSpace)))
}


## Beta
pSpace <- seq(-1, 1, length = 20)
curve(dnorm(x), -5, 5, col = "red", lwd = 2)
for(i in 1:length(pSpace)){
  curve(dstable(x, alpha = 1, beta = pSpace[i], pm = 0), add = TRUE,
        col = rgb(0, 0, i/length(pSpace)))
}

## Gamma
pSpace <- seq(0.1, 5, length = 50)
curve(dnorm(x), -5, 5, col = "red", lwd = 2)
for(i in 1:length(pSpace)){
  curve(dstable(x, alpha = 1, beta = 1, gamma = pSpace[i], pm = 0), add = TRUE,
        col = rgb(0, 0, i/length(pSpace)))
}


## Delta
pSpace <- seq(1e-5, 1, length = 20)
curve(dnorm(x), -5, 5, col = "red", lwd = 2)
for(i in 1:length(pSpace)){
  curve(dstable(x, alpha = 1, beta = 1, delta = pSpace[i], pm = 0), add = TRUE,
        col = rgb(0, 0, i/length(pSpace)))
}

## This is the distribution we can generate from
curve(dnorm(x), -5, 5, col = "red", lwd = 2)
curve(dstd(x, nu = 4), add = TRUE, col = "green")
curve(dstable(x, alpha = 1.6, beta = 0, gamma = 0.7), add = TRUE, col = "blue")
curve(dstable(x, alpha = 1.57, beta = 0.159, gamma = 6.76 * 10^-3,
              delta = 3.5*10^-3), col = "purple", add = TRUE)




## Old mgarchSim just for backup
## mgarchSim <- function(n, beta, mix, seed = NULL){
##     ## Function to generate GARCH time series with scale normal
##     ## mixture distribution
##     ##
##     ## Args:
##     ##   beta: [1] omega: The mean of the conditional variance equation.
##     ##         [2] beta1: The coefficient for the variance in the
##     ##                    conditional variance equation
##     ##         [3] alpha1: The coefficient of the squared observation
##     ##                     in the conditional variance equation.
##     ##         [4] sigma_0: The estimated initial value of the
##     ##                      conditional variance
##     ##         [5] xi: The skewness parameter of the distribution.
##     ##   n: The length of the simulated series to be generated.
##     ##   mix: The mixing distribution to be used.
##     ##   seed: The initial seed for the random number.
##     ##
##     ## Output:
##     ##   x: The GARCH time series with the desired mixing error
##     ##      distribution.
##     ##   sigma.t: The underlying volatility that was used to generate
##     ##            the time series.

##     if(!is.null(seed)) set.seed(seed)

##     ## Initialise the vector and their first observations.
##     sim <- double(n)
##     sigma2 <- double(n)
    
##     sigma2[1] <- beta[4]^2
##     sim[1] <- sqrt(sigma2[1]) * rssnormm(1, mix, beta[5])$sim

##     ## Loop to generate the underlying volatility first and then the
##     ## observations.
##     for(i in 2:n){
##         sigma2[i] = beta[1] + beta[2] * sigma2[i - 1] +
##             beta[3] * sim[i - 1]^2
##         sim[i] = sqrt(sigma2[i]) * rssnormm(1, mix, beta[5])$sim
##     }
##     ## class(sim) <- "mgarch"
##     list(x = sim, sigma2 = sigma2)
## }
