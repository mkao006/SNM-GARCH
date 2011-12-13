######################################################################
## A class for mixture linear state-space is developed
######################################################################

## Will use the nyse data and the model described in Shumway & Stoffer
## as benchmark

## At Home
load("c:/Users/user/Dropbox/Uni Work/Masters Thesis/Codes/State Space Model/tsa3.rda")

## At Work
#load("~/Dropbox/Uni Work/Masters Thesis/Codes/State Space Model/tsa3.rda")

## Plot the original and the transformed time series
par(mfrow = c(2, 1))
plot(nyse)
y <- log(nyse^2)
plot(y)



## This is a class of functions that are used to implement the mixture
## CNM method in a state-space frame work

## Need an initial function, will use the original estimatin function
## provided by Shumway and Stoffer.
initial <- function(){
}

## A valid function to be implemented within lsch, it determines
## whether the parameter are valid. (e.g. variance cannot be negative)
valid <- function(){
}

## The most important function, the log density function. This
## function will be based on Kalman filter for the case of linear
## mixture state-space models, we will then extended to the non-linear
## case further down the track

## sigma assumed to be 1?
logd.msv <- function(x, beta, pt, which = c(1, 0, 0, 0)){
    dl = vector("list", 4)
    names(dl) = c("ld", "db1", "dt1", "dt2")
    if(which[1] == 1)
        eps <- x - c(kf(x, beta[1], beta[2], beta[3],
                      beta[4], beta[5], beta[6], beta[7])$at[-1])
        dl$ld = -log(2 * pi) - (outer(eps, pt, "-")^2)/2
    if(which[2] == 1)
        dl$db1 = b
    if(which[3] == 1)
        dl$dt1 = outer(eps, pt, "-")
    if(which[4] == 1)
        ## equal to the variance, can be estimated from the Kalman
        ## Filteras well
        dl$dt2 = 1
    dl
}


### NOTES:
## The kalman filter will be used as part of the logd function, maybe
## we can test which of the kalman filter is better? We will use kf
## first as it is implemented in Fotran and it returns all the desired
## statistics

## Test functions
test <- kf(y, 1, 1, 1, 1, 1, 1, 1)
plot(y)
lines(c(test$at), col = "red")

plot(c(y) - c(test$at), type = "l")
hist(c(y) - c(test$at), breaks = 100)

hist(diff(y), freq = FALSE, breaks = 500)
curve(dnorm(x, mean(diff(y)), sd(diff(y))), add = TRUE, col = "red")
plot(y)plot(y)plot(y)plot(y)

lines(c(test$at), col = "red")

lines(c(test$at), col = "red")

args(kf)
function (yt, Zt, Tt, Rt, Ht, Qt, a1, P1, P1inf = 0, optcal = c(TRUE,
    TRUE, TRUE, TRUE), tol = 1e-07)

apply(logd.msv(y, beta = rep(1, 7), pt = c(-3:3))[[1]], 2, sum)

apply(logd.msv(y, beta = rep(1, 7),
               pt = c(-3:3), which = c(1, 0, 1, 0))[[3]], 2, sum)



kf(y[-(1:10)], 1, 1, 1, 1, 1, 1, 1)$lik +
kf(y[-c(1:10)], 1, 1, 1, 1, 1, kf(y[1:10], 1, 1, 1, 1, 1, 1, 1)$at[, 11],
   tail(kf(y[1:10], 1, 1, 1, 1, 1, 1, 1)$Pt, 1))$lik

kf(y, 1, 1, 1, 1, 1, 1, 1)$lik


######################################################################
## Books to read
######################################################################

## (1)Financial risk management with Bayesian estimation of GARCH
## models theory and applications

