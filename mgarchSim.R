######################################################################
## Functions to support the simulation study
######################################################################


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
    sim <- rsnorm(n, mean = 0, sd = mix$pt[m], xi = xi)
    list(m = m, sim = sim)
}

## check <- rssnormm(1000, mix = disc(c(0, 5, 10), c(0.5, 0.25, 0.25)), xi = 1)
## table(check$m)
## hist(check$sim, breaks = 50)


## TODO (Michael): Rewrite this function so that the model works the
## same way as what is adopted in ARIMA


mgarchSim <- function(beta, n, mix, seed = NULL){
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
    sim[1] <- sqrt(sigma2[1]) * rssnormm(1, mix, beta[5])$sim

    ## Loop to generate the underlying volatility first and then the
    ## observations.
    for(i in 2:n){
        sigma2[i] = beta[1] + beta[2] * sigma2[i - 1] +
            beta[3] * sim[i - 1]^2
        sim[i] = sqrt(sigma2[i]) * rssnormm(1, mix, beta[5])$sim
    }
    ## class(sim) <- "mgarch"
    list(x = sim, sigma2 = sigma2)
}

mysim <- mgarchSim(beta = c(1e-3, 0.8, 0.1, 0.001, 2), n = 1000,
                   mix = disc(c(1, 1.2, 1.5), rep(1/3, 3)))
plot(mysim$x, type = "l")



## NOTES (Michael): This section is deprecated due to the fact that the
## beta have extended from alpha1, beta1, mu, omega to include the
## initial variance(sigma_1^2) and the skewness (xi)

## mgarchSim <- function(beta, n, mix, seed = 587){
##     ## Function to generate GARCH time series with scale normal
##     ## mixture distribution
##     ##
##     ## Args:
##     ##   beta: A vector containing the parameters of the beta.
##     ##          Currently the mu, alpha, beta are supported.
##     ##   n: The length of the simulated series to be generated.
##     ##   mix: The mixing distribution to be used.
##     ##   seed: The initial seed for the random number.
##     ##
##     ## Output:
##     ##   x: The GARCH time series with the desired mixing error
##     ##      distribution.
##     ##   sigma.t: The underlying volatility that was used to generate
##     ##            the time series.
##     set.seed(seed)
##     ## Initialise parameters
##     sim <- double(n)
##     sigma.t <- double(n)
##     ## Use the theoretical asymptotic variance as the initial value of
##     ## the conditional standard deviation
##     sigma.t[1] <- beta[1]/(1 - beta[2])
##     ## Initialise the first observation
##     sim[1] <- beta[4] + sqrt(sigma.t[1]) * rsnorm(1, mix)$sim
##     ## Loop to generate the underlying volatility and then the
##     ## observations.
##     for(i in 2:n){
##         sigma.t[i] = beta[1] + beta[2] * sigma.t[i - 1] +
##             beta[3] * sim[i - 1]^2
##         sim[i] = beta[4] + sqrt(sigma.t[i]) * rsnorm(1, mix)$sim
##     }
##     class(sim) <- "mgarch"
##     list(x = sim, sigma.t = sigma.t)
## }
## mysim <- mgarchSim(beta = c(1e-3, 0.8, 0.1), n = 1000,
##                    mix = disc(c(1, 1.2, 1.5), rep(1/3, 3)))
## plot(mysim$x, type = "l")
