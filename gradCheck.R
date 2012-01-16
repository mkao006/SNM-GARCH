######################################################################
## Check the derivatives
######################################################################

## These functions serves as support to check whether the derivatives of
## the logd functions are calculated correctly

set.seed(587)
Data <- rnorm(100, 0, 1)
## Data <- rnorm(100, 0, 1e-3)
## Data <- rnorm(100, 0, 1e-5)
## Data <- rnorm(100, 0, 1e-10)
class(Data) <- "mgarch"

## Data the derivatives w.r.t the betas
beta.diff <- function(x = Data,
                      beta = c(0.001, 0.8, 0.1, sd(x), 2), incremt = 1e-6,
                      pt = c(1, 2)) {
    index0 = rep(0, 5)
    for(i in 1:5) {
        index = index0
        index[i] = 1
        d1 = (logd.mgarch(x, pt = pt, beta = beta + index * incremt,
                          which=c(1, 0, 0, 0))$ld -
              logd.mgarch(x, pt = pt, beta = beta - index * incremt,
                          which = c(1, 0, 0, 0))$ld)/(2 * incremt)
    d2 = logd.mgarch(x, pt = pt, beta = beta, which = c(0, 1, 0, 0))$db1[,,i]
    print(max(abs(d1 - d2)))
    }
#    list(nd = d1, ad = d2)
}

beta.diff()

## Data the derivatives w.r.t the support points (theta)
theta1.diff <- function(x = Data,
                      beta = c(0.001, 0.8, 0.1, sd(x), 2), incremt = 1e-6,
                      pt = c(1, 2)) {
  d1 = (logd.mgarch(x, pt = pt + incremt, beta = beta,
                    which = c(1, 0, 0, 0))$ld -
        logd.mgarch(x, pt = pt - incremt, beta = beta,
                    which = c(1, 0, 0, 0))$ld)/(2 * incremt)
  d2 = logd.mgarch(x, pt = pt, beta = beta, which = c(0, 0, 1, 0))$dt
  print(max(abs(d1 - d2)))
#  list(nd = d1, ad = d2)
}

theta1.diff()


## Deprecated, since second order derivative is not required in the new
## implementation of the cnmms algorithm.
## ## Data the 2nd order derivatives w.r.t the support points (theta)
## theta2.diff <- function(x = Data,
##                       beta = c(0.001, 0.8, 0.1, sd(x), 2), incremt = 1e-6,
##                       pt = c(1, 2)) {
##   d1 = (logd.mgarch(x, pt = pt + incremt, beta = beta,
##                     which = c(1, 0, 0, 0))$ld +  ## f(x + h)
##         logd.mgarch(x, pt = pt - incremt, beta = beta,
##                     which = c(1, 0, 0, 0))$ld -  ## f(x - h)
##                       (2 * logd.mgarch(x, pt = pt, beta = beta,
##                     which = c(1, 0, 0, 0))$ld))/ ## f(x)
##                       (incremt^2)
##   d2 = logd.mgarch(x, pt = pt, beta = beta, which = c(0, 0, 0, 1))$dt2
##   print(max(abs(d1 - d2)))
## #  list(nd = d1, ad = d2)
## }

## theta2.diff()
