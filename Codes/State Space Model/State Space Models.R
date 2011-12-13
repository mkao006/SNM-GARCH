library(sspir)
library(KFAS)

######################################################################
## This is a demonstration of how the Kalman Filter recovers the state
## series most of the time, however the algorithm fails sometimes and
## therefore it would be important to ues methods like bootstrap to
## verify the estimates.
######################################################################

######################################################################
## This part demonstrates that the same script does not work on some
## parameters
######################################################################

## (a) Fail to recover state
## Initialise the parameters
set.seed(123)
N = 300
Phi = 0.3
A = 2
sim.x <- arima.sim(N, model = list(ar = Phi, sd = 1))
sim.y <- A * sim.x + rnorm(N, 0, 3)

## the plot of the simulated data
plot(sim.y, col = "red")
lines(sim.x)

real.par <- c(Phi, A, 1, 3)

## Estimate the par through Kalman filter MLE
logLik.kf <- function(theta){
    lik <- kf(yt = sim.y, Zt = theta[3], Tt = theta[1],
              Rt = 0, Ht = theta[2], Qt = 1,
              a1 = 0, P1 = 1e7)
    -lik$lik
}

(sim.fit <- optim(par = c(1, 1, 1), fn = logLik.kf,
                  lower = rep(1e-4, 4), method = "L-BFGS-B"))


## Plot the fit and probability limits
filtersim <- kf(yt = sim.y, Zt = sim.fit$par[3],
                Tt = sim.fit$par[1], Rt = 0, Ht = sim.fit$par[2],
                Qt = 1, a1 = 0, P1 = 1e7)
smoothsim <- ks(filtersim)
attach(smoothsim)
hwidth <- qnorm(0.05, lower = FALSE) * sqrt(drop(Vt))
sm <- cbind(drop(ahat), as.vector(ahat) + hwidth %o% c(-1, 1))
sm <- ts(sm, start = start(sim.y))
plot(sm, plot.type = "s", type = "l", lty = c(1, 5, 5),
     ylab = "Level", xlab = "", ylim = range(sim.y))
lines(sim.y, type = "o", col = "darkgrey")
legend("bottomleft", col = c("darkgrey", rep("black", 2), "red"),
       lty = c(1, 1, 5, 1), pch = c(1, NA, NA, NA), bty = "n", legend =
       c("data", "smoothed level", "90% probability limits", "state"))
detach(smoothsim)
lines(sim.x, col = "red")

## (b) Success in recovering the state
## Initialise the parameters
set.seed(587)
N = 300
Phi = 0.3
A = 2
sim.x <- arima.sim(N, model = list(ar = Phi, sd = 1))
sim.y <- A * sim.x + rnorm(N, 0, 3)

## the plot of the simulated data
plot(sim.y, col = "red")
lines(sim.x)

real.par <- c(Phi, A, 1, 3)

## Estimate the par through Kalman filter MLE
logLik.kf <- function(theta){
    lik <- kf(yt = sim.y, Zt = theta[3], Tt = theta[1],
              Rt = 0, Ht = theta[2], Qt = 1,
              a1 = 0, P1 = 1e7)
    -lik$lik
}

(sim.fit <- optim(par = c(1, 1, 1), fn = logLik.kf,
                  lower = rep(1e-4, 4), method = "L-BFGS-B"))


## Plot the fit and probability limits
filtersim <- kf(yt = sim.y, Zt = sim.fit$par[3],
                Tt = sim.fit$par[1], Rt = 0, Ht = sim.fit$par[2],
                Qt = 1, a1 = 0, P1 = 1e7)
smoothsim <- ks(filtersim)
attach(smoothsim)
hwidth <- qnorm(0.05, lower = FALSE) * sqrt(drop(Vt))
sm <- cbind(drop(ahat), as.vector(ahat) + hwidth %o% c(-1, 1))
sm <- ts(sm, start = start(sim.y))
plot(sm, plot.type = "s", type = "l", lty = c(1, 5, 5),
     ylab = "Level", xlab = "", ylim = range(sim.y))
lines(sim.y, type = "o", col = "darkgrey")
legend("bottomleft", col = c("darkgrey", rep("black", 2), "red"),
       lty = c(1, 1, 5, 1), pch = c(1, NA, NA, NA), bty = "n", legend =
       c("data", "smoothed level", "90% probability limits", "state"))
detach(smoothsim)
lines(sim.x, col = "red")




## This is employing the sspir package
## Ignore this part
#sim.yts <- ts(matrix(sim.y), freq = 12)
#y.ss <- SS(sim.yts)
#phi(y.ss) <- c(1, 0)
#m0(y.ss) <- matrix(5)
#C0(y.ss) <- matrix(5)

#y.kfit <- kfilter(y.ss)
#plot(y.kfit$y)
#lines(y.kfit$m, lty = 2, col = "red")


######################################################################
## This part is just some test that I ran to see the effect of changin
## the number of parameters estimated, different error distribution
## etc....
######################################################################

## Currently it shows that if I tried to estimate the Qt (state error
## variance) on top of the Zt (Observation transition matrix) and Tt
## (state transition matrix), the estimation collapse, reason unknow.

N = 300
Phi = 0.3
A = 0.5
sim.x <- arima.sim(N, model = list(ar = Phi, sd = 1))
sim.y <- A * sim.x + rnorm(N, 0, 3)

## Lets try adjusted gamma error
#gammaerror <- rgamma(N, 2)
#gammaerror <- gammaerror - mean(gammaerror)
#sim.y <- A * sim.x + gammaerror


## the plot of the simulated data
plot(sim.y, col = "red")
lines(sim.x)


logLik.kf <- function(theta){
    ll <- kf(yt = sim.y, Zt = theta[1], Tt = theta[2], Rt = 1,
             Ht = 3, Qt = theta[3], a1 = 0, P1 = 1e7)
    -ll$lik
}

sim.fit <- optim(par = c(1, 1, 1), fn = logLik.kf)

filtersim <- kf(yt = sim.y, Zt = sim.fit$par[1],
                Tt = sim.fit$par[2], Rt = 0, Ht = 3,
                Qt = sim.fit$par[3], a1 = 0, P1 = 1e7)
smoothsim <- ks(filtersim)
attach(smoothsim)
hwidth <- qnorm(0.05, lower = FALSE) * sqrt(drop(Vt))
sm <- cbind(drop(ahat), as.vector(ahat) + hwidth %o% c(-1, 1))
sm <- ts(sm, start = start(sim.y))
plot(sm, plot.type = "s", type = "l", lty = c(1, 5, 5),
     ylab = "Level", xlab = "", ylim = range(sim.y))
lines(sim.y, type = "o", col = "darkgrey")
legend("bottomleft", col = c("darkgrey", rep("black", 2), "red"),
       lty = c(1, 1, 5, 1), pch = c(1, NA, NA, NA), bty = "n", legend =
       c("data", "smoothed level", "90% probability limits", "state"))
detach(smoothsim)
lines(sim.x, col = "red")





######################################################################
## Trying to improve the stability by taking 300 sample generated from
## the same ARIMA process using same coefficients
######################################################################

## Repeated estimation with 300 sample generated using the same
## parameters
simmean <- matrix(nc = 2, nr = 300)
for(i in 1:300){
    sim.x <- arima.sim(N, model = list(ar = Phi, sd = 1))
    sim.y <- A * sim.x + rnorm(N, 0, 3)
    simmean[i, ] <- optim(par = c(1, 1), fn = logLik.kf)$par
}
simmean.par <- apply(simmean, 2, mean)

## The corresponding plot
filtersim <- kf(yt = sim.y, Zt = simmean.par$par[1],
                Tt = simmean.par$par[2], Rt = 0, Ht = 3,
                Qt = 1, a1 = 0, P1 = 1e7)
smoothsim <- ks(filtersim)
attach(smoothsim)
hwidth <- qnorm(0.05, lower = FALSE) * sqrt(drop(Vt))
sm <- cbind(drop(ahat), as.vector(ahat) + hwidth %o% c(-1, 1))
sm <- ts(sm, start = start(sim.y))
plot(sm, plot.type = "s", type = "l", lty = c(1, 5, 5),
     ylab = "Level", xlab = "", ylim = range(sim.y))
lines(sim.y, type = "o", col = "darkgrey")
legend("bottomleft", col = c("darkgrey", rep("black", 2), "red"),
       lty = c(1, 1, 5, 1), pch = c(1, NA, NA, NA), bty = "n", legend =
       c("data", "smoothed level", "90% probability limits", "state"))
detach(smoothsim)
lines(sim.x, col = "red")


## The histogram of the fit
est <- c("Rt", "Tt", "Ht", "Qt")
par(mfrow = c(2, 1))
for(i in 1:2){
    hist(simmean[, i], breaks = 50, main = est[i],
         xlab = est[i])
}


