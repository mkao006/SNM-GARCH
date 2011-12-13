source("functions.R")

######################################################################
## DEM2GBP
######################################################################

## Example for the model

tdg <- as.numeric(data.matrix(dem2gbp))[1:200]
tmdg <- tdg
class(tmdg) <- "mgarch"
(dg.mg <- cnmms(tmdg - coef(dg.t0)[1],
                plot = "gradient", grid = 100, verb = 4))

dg.t0 <- garchFit(data = tdg, cond.dist = "norm")
dg.t1 <- garchFit(data = tdg, cond.dist = "ged")
dg.t2 <- garchFit(data = tdg, cond.dist = "std")

dg.t0@fit$ll
dg.t1@fit$ll
dg.t2@fit$ll

check.ll <- function(beta){
    -sum(logd.mgarch(tdg - coef(dg.t0)[1], beta = beta,
                     pt = 1, which = c(1, 0, 0, 0))$ld)
}

optim(rep(0.4, 4), check.ll)


######################################################################
## Use semi-parametric model to estimate linear regression to test
## what is going wrong
######################################################################

logd.tlm <- function(x, beta, pt, which){
    y <- x$y
    x <- x$x
    lpt <- length(pt)
    dl <- vector("list", length = 4)
    names(dl) <- c("ld", "db1", "dt1", "dt2")
    if(which[1] == 1){
        dl$ld <- matrix(-1/2 * log(2 * pi) -
            log(pt) - ((y - x * beta)^2)/
                        (2 * rep(pt, each = length(y))^2), nc = lpt)
    }
    if(which[2] == 1){
        dl$db1 <- array((-(y - x * beta) * x)/
                         rep(pt, each = length(y))^2,
                        dim = c(length(y), lpt, 1))
    }
    if(which[3] == 1){
        dl$dt1 <- matrix(-1/pt + ((y - x * beta)^2)/
                         rep(pt, each = length(y))^3, nc = lpt)
    }
    if(which[4] == 1){
        dl$dt2 <- matrix(1/pt^2 - (3 * (y - x * beta)^2)/
                         rep(pt, each = length(y))^4, nc = lpt)
    }
    dl
}



valid.tlm <- function(x, beta, mix) TRUE
initial.tlm <- function(x, beta = NULL, mix = NULL, kmax = NULL){
    if(is.null(beta)){
        beta = c(coef(lm(x$y ~ x$x)))[2]
        mix <- dden(1, 1)
        list(beta = beta, mix = mix)
    }
}

x <- 1:100
y <- c(3 * x + rnorm(100, 0, 20))
plot(x, y)

mydf <- data.frame(y = y, x = x)
class(mydf) <- "tlm"
cnmms(mydf, grid = 100, verb = 4)

