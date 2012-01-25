######################################################################
## This script checks whether the logd.mgarch function is implemented
## correctly
######################################################################

set.seed(587)
## Data <- rnorm(100, 0, 1)
Data <- rnorm(10, 0, 1e-10)
class(Data) <- "mgarch"

## Data the derivatives w.r.t the betas
beta.diff <- function(x = Data,
                      beta = c(1e-4, 0.2, 0.8, 1e-8, 2),
                      incremt = 1e-10,
                      pt = c(1, 2), which.beta = 1:length(beta),
                      out = c("print", "list")) {
    out = match.arg(out)
    index0 = rep(0, 5)
    for(i in which.beta) {
        index = index0
        index[i] = 1
        d1 = (logd.mgarch(x, pt = pt, beta = beta + index * incremt,
                          which=c(1, 0, 0))$ld -
              logd.mgarch(x, pt = pt, beta = beta - index * incremt,
                          which = c(1, 0, 0))$ld)/(2 * incremt)
    d2 = logd.mgarch(x, pt = pt, beta = beta, which = c(0, 1, 0))$db[,,i]
    if(out == "print") print(max(abs(d1 - d2)))
    }
    if(out == "list") list(nd = d1, ad = d2)
}
cat("Checking derivatives of Beta (Should all be zero):\n")
beta.diff()


## Data the derivatives w.r.t the support points (theta)
theta1.diff <- function(x = Data,
                      beta = c(1e-4, 0.2, 0.8, 1e-8, 2),
                      incremt = 1e-10,
                      pt = c(1, 2)) {
  d1 = (logd.mgarch(x, pt = pt + incremt, beta = beta,
                    which = c(1, 0, 0, 0))$ld -
        logd.mgarch(x, pt = pt - incremt, beta = beta,
                    which = c(1, 0, 0, 0))$ld)/(2 * incremt)
  d2 = logd.mgarch(x, pt = pt, beta = beta, which = c(0, 0, 1, 0))$dt
  print(max(abs(d1 - d2)))
#  list(nd = d1, ad = d2)
}

cat("Checking derivatives of Theta (Should all be zero):\n")
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


## Check whether the log-likelihood matches that of the fGarch package

#### dem2gbp data
## The vanilla normal case
check <- garchFit(data = dem2gbp, trace = FALSE, include.mean = FALSE)
checkData <- c(unlist(dem2gbp))
inits <- initial.mgarch(checkData)
cat("Log-likelihoods (Both value should be the same):\n")
print(sum(logd.mgarch(checkData, inits$beta, inits$mix$pt,
                      which = c(1, 0, 0))$ld))
print(-check@fit$ll)

## The skewed normal case
check <- garchFit(data = dem2gbp, trace = FALSE,
                  include.mean = FALSE, cond.dist = "snorm")
inits <- list(beta = c(coef(check)[c(1, 3, 2)],
                check@sigma.t[1], coef(check)[4]),
              mix = disc(1, 1))
cat("Log-likelihoods (Both value should be the same):\n")
print(sum(logd.mgarch(checkData, inits$beta, inits$mix$pt,
                      which = c(1, 0, 0))$ld))
print(-check@fit$ll)

## #### S&P data
## ## The vanilla normal case
## check <- garchFit(data = sp, trace = FALSE, include.mean = FALSE)
## checkData <- c(unlist(sp))
## inits <- initial.mgarch(checkData)
## print(sum(logd.mgarch(checkData, inits$beta, inits$mix$pt,
##                       which = c(1, 0, 0))$ld))
## print(check@fit$ll)

## ## The skewed normal case
## check <- garchFit(data = sp, trace = FALSE,
##                   include.mean = FALSE, cond.dist = "snorm")
## inits <- list(beta = c(coef(check)[c(1, 3, 2)],
##                 check@sigma.t[1], coef(check)[4]),
##               mix = disc(1, 1))
## print(sum(logd.mgarch(checkData, inits$beta, inits$mix$pt,
##                       which = c(1, 0, 0))$ld))
## print(check@fit$ll)
