######################################################################
## Functions to support the simulation study
######################################################################

rsnorm <- function(n, mix){
    ## A function to generate samples from a scale normal mixture
    ## distribution.
    ##
    ## Args:
    ##   n: The number of observation to be generated as used in other
    ##   "rdist" functions
    ##   mix: The mixture function generated using the dden function
    ##
    ## Output:
    ##   m: The index of which component distribution the random
    ##      number was drawn from. (e.g. 1 equal the first component
    ##      density etc...)
    ##   sim: The list of random numbers generated.

    ## To determin which distribution to draw from
    w <- runif(n)

    ## Select the right component distribution and generate the random
    ## sample
    breaks = cumsum(mix$pr)

    ## TODO: Try vectorize this
    m <- double(n)
    for(i in 1:n){
        m[i] = sum(w[i] > breaks) + 1
    }
    sim <- rnorm(n, 0, mix$pt[m])
    list(m = m, sim = sim)

}


# check <- rmnorm(1000, dden(c(0, 5, 10), c(0.5, 0.25, 0.25)))
# table(check$m)
# hist(check$sim, breaks = 100)


## Rewrite this function so that the model works the same way as what
## is adopted in ARIMA
mgarchSim <- function(model, n, mix, seed = 587){
    ## Function to generate GARCH time series with scale normal
    ## mixture distribution
    ##
    ## Args:
    ##   model: A vector containing the parameters of the model.
    ##          Currently the mu, alpha, beta are supported.
    ##   n: The length of the simulated series to be generated.
    ##   mix: The mixing distribution to be used.
    ##   seed: The initial seed for the random number.
    ##
    ## Output:
    ##   x: The GARCH time series with the desired mixing error
    ##      distribution.
    ##   sigma.t: The underlying volatility that was used to generate
    ##            the time series.

    set.seed(seed)

    ## Initialise parameters
    sim <- double(n)
    sigma.t <- double(n)


    ## use the theoretical asymptotic variance as the initial value
    sigma.t[1] <- model[1]/(1 - model[2])

    ## just simulate the data
    sim[1] <- model[4] + sqrt(sigma.t[1]) * rsnorm(1, mix)$sim

    ## Loop to generate the underlying volatility and then the
    ## observations.
    for(i in 2:n){
        sigma.t[i] = model[1] + model[2] * sigma.t[i - 1] +
            model[3] * sim[i - 1]^2
        sim[i] = model[4] + sqrt(sigma.t[i]) * rsnorm(1, mix)$sim
    }
    class(sim) <- "mgarch"
    list(x = sim, sigma.t = sigma.t)
}

#mysim <- mgarchSim(model = c(1e-3, 0.8, 0.1), n = 1000,
#                   mix = dden(c(1, 1.2, 1.5), rep(1/3, 3)))
#plot(mysim$x, type = "l")
