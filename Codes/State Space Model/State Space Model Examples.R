## Example 6.5
set.seed(1)
num = 50
w = rnorm(num + 1, 0, 1)
v = rnorm(num, 0, 1)
mu = cumsum(w) # states: mu[0], ..., mu[50]
y = mu[-1] + 1 # obs: y[1], ..., y[50]

# filter and smooth (Ksmooth0 does both)
mu0 = 0
sigma0 = 1
phi = 1
cQ= 1
cR = 1
ks = Ksmooth0(num, y, 1, mu0, sigma0, phi, cQ, cR)

## Plots
Time = 1:num
par(mfrow = c(3, 1))

plot(Time, mu[-1], main = "Prediction", ylim = c(-5, 10))
lines(ks$xp)
lines(ks$xp + 2 * sqrt(ks$Pp), lty = "dashed", col = "blue")
lines(ks$xp - 2 * sqrt(ks$Pp), lty = "dashed", col = "blue")

plot(Time, mu[-1], main = "Filter", ylim = c(-5, 10))
lines(ks$xf)
lines(ks$xf + 2 * sqrt(ks$Pf), lty = "dashed", col = "blue")
lines(ks$xf - 2 * sqrt(ks$Pf), lty = "dashed", col = "blue")

plot(Time, mu[-1], main = "Smoother", ylim = c(-5, 10))
lines(ks$xs)
lines(ks$xs + 2 * sqrt(ks$Ps), lty = "dashed", col = "blue")
lines(ks$xs - 2 * sqrt(ks$Ps), lty = "dashed", col = "blue")

## Parameters
mu[1]
ks$x0n
sqrt(ks$P0n)

dev.new()
plot(Time, mu[-1], ylim = c(-2, 7))
lines(ks$xp, col = 4)
lines(ks$xf, col = 3)
lines(ks$xs, col = 2)
names = c("Predictor", "Filter", "Smoother")
legend("bottomright", names, col = 4:2, lty = 1, bty = "n")

## Example 6.6
set.seed(999)
num = 100
N = num + 1
x = arima.sim(n = N, list(ar = .8, sd = 1))
y = ts(x[-1] + rnorm(num, 0, 1))

## Initial estimates
u = ts.intersect(y, lag(y, -1), lag(y, -2))
varu <- var(u)
coru <- cor(u)
phi <- coru[1, 3]/coru[1, 2]
q <- (1 - phi^2) * varu[1, 2]/phi
r <- varu[1, 1] - q/(1 - phi^2)
(init.par = c(phi, sqrt(q), sqrt(r)))

## Function to evaluate the likelihood
Linn <- function(para){
    phi <- para[1]
    sigw <- para[2]
    sigv <- para[3]
    Sigma0 <- (sigw^2)/(1 - phi^2)
    Sigma0[Sigma0 < 0] = 0
    kf <- Kfilter0(num, y, 1, mu0 = 0, Sigma0,
                   phi, sigw, sigv)
    kf$like
}

## Estimation
(est <- optim(init.par, Linn, gr = NULL,
              method = "BFGS", hessian = TRUE,
              control = list(trace = 1, REPORT = 1)))
SE = sqrt(diag(solve(est$hessian)))
cbind(estimate = c(phi = est$par[1], sigw = est$par[2], sigv = est$par[3]), SE)


########################################################################
## The stats space model for the log of the NYSE (6.9)
########################################################################

y <- log(nyse^2)
num <- length(y)

## Initial parameters
phi0 = 0
phi1 = 0.95
sQ = 0.2
alpha = mean(y)
sR0 = 1
mu1 = -3
sR1 = 2
init.par = c(phi0, phi1, sQ, alpha, sR0, mu1, sR1)

## Innovations Likelihood
Linn <- function(para){
    phi0 = para[1]
    phi1 = para[2]
    sQ = para[3]
    alpha = para[4]
    sR0 = para[5]
    mu1 = para[6]
    sR1 = para[7]
    sv = SVfilter(num, y, phi0, phi1, sQ, alpha, sR0, mu1, sR1)
    sv$like
}

## Estimation
est <- optim(init.par, Linn, NULL, method = "BFGS", hessian = TRUE,
             control = list(trace = 1, REPORT = 1))
est
SE = sqrt(diag(solve(est$hessian)))
u <- cbind(estimates = est$par, SE)
rownames(u) <- c("phi0", "phi1", "sQ", "alpha", "sigv0", "mu1", "sigv1")

## Graphics

phi0 <- est$par[1]
phi1 <- est$par[2]
sQ <- est$par[3]
alpha <- est$par[4]
sR0 <- est$par[5]
mu1 <- est$par[6]
sR1 <- est$par[7]

sv <- SVfilter(num, y, phi0, phi1, sQ, alpha, sR0, mu1, sR1)

Time <- 801:1000
## Density plot
x <- seq(-15, 6, by = 0.1)
f <- exp(-0.5 * (exp(x) - x))/(sqrt(2 * pi))
f0 <- exp(-0.5 * (x^2)/sR0^2)/(sR0 * sqrt(2 * pi))
f1 <- exp(-0.5 * (x - mu1)^2/sR1^2)/(sR1 * sqrt(2 * pi))
fm = (f0 + f1)/2
plot(x, f, type = "l")
lines(x, fm, lty = 2, lwd = 2)
dev.new()
par(mfrow = c(2, 1))
plot(Time, y[Time], type = "l", main = "log(Squared NYSE Returns)")
plot(Time, sv$xp[Time], type = "l", main = "Predicted log-Volatility",
     ylim = c(-1.5, 1.8), ylab = "", xlab = "")
lines(Time, sv$xp[Time] + 2 * sqrt(sv$Pp[Time]), lty = 2)
lines(Time, sv$xp[Time] - 2 * sqrt(sv$Pp[Time]), lty = 2)


########################################################################
## Example of the Kalman filter (6.2)
########################################################################

## Generate the data
set.seed(1)

num <- 50
w = rnorm(num + 1)
v = rnorm(num)
mu <- cumsum(w)
y <- 0.8 *  mu[-1] + v
plot(y, type = "b")

## Filter and Smooth
mu0 = 0
sigma0 = 1
phi = 0.1
cQ = 1
cR = 1
ks <- Ksmooth0(num, y, 1, mu0, sigma0, phi, cQ, cR)

## Start figure
par(mfrow = c(3, 1))
Time = 1:num
plot(Time, mu[-1], main = "Prediction", ylim = c(-5, 10))
lines(ks$xp)
lines(ks$xp + 2 * sqrt(ks$Pp), lty = 2, col = "blue")
lines(ks$xp - 2 * sqrt(ks$Pp), lty = 2, col = "blue")
lines(y, col = "red")

plot(Time, mu[-1], main = "Filter", ylim = c(-5, 10))
lines(ks$xf)
lines(ks$xf + 2 * sqrt(ks$Pp), lty = 2, col = "blue")
lines(ks$xf - 2 * sqrt(ks$Pp), lty = 2, col = "blue")
lines(y, col = "red")

plot(Time, mu[-1], main = "Smoother", ylim = c(-5, 10))
lines(ks$xs)
lines(ks$xs + 2 * sqrt(ks$Pp), lty = 2, col = "blue")
lines(ks$xs - 2 * sqrt(ks$Pp), lty = 2, col = "blue")
lines(y, col = "red")


## MLE of state state space model
# Generate Data
set.seed(999)

num = 200
N = num + 1
x = arima.sim(n = N, list(ar = .3, sd = 0.1))
y = ts(x[-1] + rnorm(num,0,1))

# Initial Estimates
u = ts.intersect(y, lag(y,-1), lag(y,-2))
varu = var(u)
coru = cor(u)
phi = coru[1,3]/coru[1,2]
q = (1-phi^2)*varu[1,2]/phi
r = varu[1,1] - q/(1-phi^2)
(init.par = c(phi, sqrt(q), sqrt(r))) # = .91, .51, 1.03

# Function to evaluate the likelihood
Linn <- function(para){
    phi = para[1]
    sigw = para[2]
    sigv = para[3]
    #Ahat = para[4]
    Sigma0 = (sigw^2)/(1-phi^2)
    Sigma0[Sigma0<0]=0
    kf = Kfilter0(num, y, A = 1, mu0=0, Sigma0, phi, sigw, sigv)
    return(kf$like)
}

# Estimation (partial output shown)
est = optim(init.par, Linn, gr=NULL, method="BFGS", hessian=TRUE)
SE = sqrt(diag(solve(est$hessian)))
c(0.3, 1, 1)
est$par

fit <- Kfilter0(num, y, A = 1, mu0 = 0, Sigma0 = est$par[2]^2/(1 - est$par[1]^2),
                est$par[1], est$par[2], est$par[3])
KFstate <- c(matrix(fit$xf))

## Plot the fit and the original
#par(mfrow = c(4, 1))
plot(y)
lines(x, col = "red", lty = 2)
lines(KFstate, col = "blue", lty = 2)
legend("topleft", legend = c("Observed", "State", "Kalman-Filtered State"),
       col = c("black", "red", "blue"), lty = c(1, 2, 2), bty = "n")

plot(x - KFstate)
abline(h = 0, col = "red", lty = 2)

acf(x - KFstate)
pacf(x - KFstate)





## Second example of MLE
# Setup

y = cbind(gtemp,gtemp2)

N = 100
dft <- cumsum(rep(0.2, N))
state <- arima.sim(n = N, list(ar = .5, sd = .1)) + dft
ob1 <- state + rnorm(N, 0, 1)
ob2 <- state + rnorm(N, 0, 1)

plot(state)
lines(ob1, col = "red")
lines(ob2, col = "blue")

y <- cbind(ob1, ob2)
num = nrow(y)
input = rep(1,num)
A = array(rep(1,2), dim=c(2,1,num))
mu0 = 0
Sigma0 = 1
Phi = 1

# Function to Calculate Likelihood
Linn <- function(para){
    cQ = para[1] # sigma_w
    cR1 = para[2] # 11 element of chol(R)
    cR2 = para[3] # 22 element of chol(R)
    cR12 = para[4] # 12 element of chol(R)
    cR = matrix(c(cR1,0,cR12,cR2),2) # put the matrix together
    drift = para[5]
    Phihat = para[6]
    kf = Kfilter1(num,y,A,mu0,Sigma0,Phi = Phihat,drift,0,cQ,cR,input)
    return(kf$like)
}

# Estimation
init.par = c(.1,.1,.1,0,.05, 1) # initial values of parameters
(est <- optim(init.par, Linn, NULL, method="BFGS", hessian=TRUE,
              control=list(trace=1,REPORT=1, maxit = 1000)))

#SE = sqrt(diag(solve(est$hessian)))

# display estimates
#u = cbind(estimate=est$par, SE)
#rownames(u)=c("sigw","cR11", "cR22", "cR12", "drift")
#u

# Smooth (first set parameters to their final estimates)
cQ=est$par[1]
cR1=est$par[2]
cR2=est$par[3]
cR12=est$par[4]
cR = matrix(c(cR1,0,cR12,cR2), 2)
R = t(cR)%*%cR # to view the estimated R matrix
drift = est$par[5]
ks = Ksmooth1(num,y,A,mu0,Sigma0,Phi,drift,0,cQ,cR,input)

#Plot
xsmooth = ts(as.vector(ks$xf), start=1)
plot(xsmooth, lwd=2, ylim = c(range(cbind(ob1, ob2))), ylab = "Temperature Deviations")
lines(ob1, col="blue", lty=1) # color helps here
lines(ob2, col="red", lty=2)
lines(lowess((ob1 + ob2)/2, f = 1/100000), lty = 2)



lines(gtemp, col="blue", lty=1) # color helps here
lines(gtemp2, col="red", lty=2)
lines(lowess((gtemp + gtemp2)/2, f = 1/100000), lty = 2)


######################################################################
## Use the dlm package
######################################################################
library(dlm)

rw <- dlm(m0 = 0, C0 = 10, FF = 1, V = 1.4, GG = 1, W = 0.2)
unlist(rw)
lg <- dlm(FF = matrix(c(1, 0), nr = 1),
          V = 1.4,
          GG = matrix(c(1, 0, 1, 1), nr = 2),
          W = diag(c(0, 0.2)),
          m0 = rep(0, 2),
          C0 = 10 * diag(2))
lg

x <- rnorm(100) # covariates
dlr <- dlm(FF = matrix(c(1, 0), nr = 1),
           V = 1.3,
           GG = diag(2),
           W = diag(c(0.4, 0.2)),
           m0 = rep(0, 2), C0 = 10 * diag(2),
           JFF = matrix(c(0, 1), nr = 1),
           X = x)
dlr



######################################################################
## GARCH and data extraction
######################################################################

library(tseries)
library(fGarch)

goog.ts <- get.hist.quote(instrument = "GOOG", start = c(2004-08-19))

GOpen.ts <-goog.ts[, 1]
plot(GOpen.ts)

DLGO.ts <- ts(diff(log(GOpen.ts)))

g.garch <- garch(DLGO.ts, order = c(1, 1))
summary(g.garch)
plot(g.garch)



######################################################################
## Examples in "State space models in R"
######################################################################

library(dlm)
library(KFAS)
library(forecast)
library(numDeriv)

fitNile <- StructTS(Nile, "level")
fitNile

plot(Nile, type = "o")
lines(fitted(fitNile), lty = 2, lwd = 2)
lines(tsSmooth(fitNile), lty = 3, lwd = 2)

plot(forecast(fitNile, level = c(seq(50, 90, by = 10)), h = 10),
     xlim = c(1950, 1980))

mod <- dlmModPoly(1, dV = 0.3, dW = 0.01)

buildNile <- function(theta){
    dlmModPoly(order = 1, dV = theta[1], dW = theta[2])
}

## Similar to fitNile from StructTS
fit <- dlmMLE(Nile, parm = c(100, 2), buildNile, lower = rep(1e-4, 2))

modNile <- buildNile(fit$par)
drop(V(modNile))
drop(W(modNile))

## Hessian of the estimates at the optim
hs <- hessian(function(x) dlmLL(Nile, buildNile(x)), fit$par)
all(eigen(hs, only.values = TRUE)$values > 0)

## Variance and standard deviation from the hessain
aVar <- solve(hs)
sqrt(diag(aVar))

## The smooth
smoothNile <- dlmSmooth(Nile, modNile)
hwidth <- qnorm(0.05, lower = FALSE) * sqrt(unlist(dlmSvd2var(smoothNile$U.S, smoothNile$D.S)))
sm <- cbind(smoothNile$s, bound = as.vector(smoothNile$s) + hwidth %o% c(-1, 1))

## Now use Kalman filter on the with the estimates obtained previously
filterNile <- dlmFilter(Nile, modNile)
plot(residuals(filterNile, sd = FALSE), type = "o",
     ylab = "Statndardized prediction error")
abline(h = 0)

## Producing the forecasts
foreNile <- dlmForecast(filterNile, nAhead = 10)
attach(foreNile)
hwidth <- qnorm(0.25, lower = FALSE) * sqrt(unlist(Q))
fore <- cbind(f, as.vector(f) + hwidth %o% c(-1, 1))
rg <- range(c(fore, window(Nile, start = c(1951, 1))))
plot(fore, type = "o", pch = 16, plot.type = "s", lty = c(1, 3, 3),
     ylab = "Nile level", xlab = "", xlim = c(1951, 1980), ylim = rg)
lines(window(Nile, start = c(1951, 1)), type = 'o')
R> lines(window(smoothNile$s, start = c(1951,1)), lty = 5)
abline(v = mean(c(time(f)[1], tail(time(Nile), 1))),
       lty = "dashed", col = "darkgrey")
legend("topleft", lty = c(1, 5, 1, 3), pch = c(1, NA, 16, 16), bty = "n",
       legend = c("observed level", "smoothed level", "forecasted level",
       "50% probability limits"))
detach(foreNile)


## Use the KFAS package
## Define the loglikelihood function obtained from kalman filter
logLik.kf <- function(theta){
    lik <- kf(yt = Nile, Zt = 1, Tt = 1, Rt = 1, Ht = theta[1],
              Qt = theta[2], a1 = 0, P1 = 1e7)
    -lik$lik
}

## Maximize the loglikelihood
fit <- optim(par = c(100, 2), fn = logLik.kf, lower = rep(1e-4, 2))

## Plot the fit and probability limits
filterNile <- kf(yt = Nile, Zt = 1, Tt = 1, Rt = 1, Ht = fit$par[1],
                 Qt = fit$par[2], a1 = 0, P1 = 1e7)
smoothNile <- ks(filterNile)
attach(smoothNile)
hwidth <- qnorm(0.05, lower = FALSE) * sqrt(drop(Vt))
sm <- cbind(drop(ahat), as.vector(ahat) + hwidth %o% c(-1, 1))
sm <- ts(sm, start = start(Nile))
plot(sm, plot.type = "s", type = "l", lty = c(1, 5, 5),
     ylab = "Level", xlab = "", ylim = range(Nile))
lines(Nile, type = "o", col = "darkgrey")
legend("bottomleft", col = c("darkgrey", rep("black", 2)),
       lty = c(1, 1, 5), pch = c(1, NA, NA), bty = "n", legend =
       c("data", "smoothed level", "90% probability limits"))
detach(smoothNile)

## residuals
residNile <- drop(filterNile$vtuni / sqrt(filterNile$Ftuni))


foreNile <- forecast(filterNile, fc = 9)
attach(foreNile)
hwidth <- qnorm(0.25, lower = FALSE) * sqrt(drop(Pt.fc))
fore <- ts(cbind(drop(at.fc), drop(at.fc) + hwidth %o%
                 c(-1, 1)), start = 1 + end(Nile)[1])
rg <- range(c(fore, window(Nile, start = c(1951, 1))))
detach(foreNile)



######################################################################
## A test on simulate some data, it appears that if you estimate more
## than two parameters using the Kalman MLE method, the state becomes
## flat (In some circumstances). In general, this has worked
## sufficiently well.
######################################################################

N = 1000
Phi = 0.3
A = 2
sim.x <- arima.sim(N, model = list(ar = Phi, sd = 1))
sim.y <- A * sim.x + rnorm(N, 0, 5)

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

######################################################################
## Bayesian approach of dlm with MCMC
######################################################################

set.seed(123)
gibbsOut <- dlmGibbsDIG(Nile, mod = dlmModPoly(1), shape.y = 0.1,
                        rate.y = 0.1, shape.theta = 0.1,
                        rate.theta = 0.1, n.sample = 10000, thin = 9)

burn <- 1:1000
attach(gibbsOut)
ts.plot(ergMean(dV[-burn]), ylab = "sample mean", xlab = "iterations",
        main = "obs variance")
ts.plot(ergMean(dW[-burn]), ylab = "sample mean", xlab = "iterations",
        main = "evolution variance")
acf(dV[-burn])
acf(dW[-burn])



plot(density(dV[-burn]), xlim = c(2000, 34000), ylab = "", main = "")
hist(dV[-burn], prob = TRUE, add = TRUE)
curve(dgamma(1/x, shape = 0.1, rate = 0.1) / x^2, lty = "dashed",
      add = TRUE)
plot(density(dW[-burn]), ylab = "", xlim = c(0, 16000), main = "")
hist(dW[-burn], prob = TRUE, add = TRUE)
curve(dgamma(1/x, shape = 0.1, rate = 0.1) / x^2, lty = "dashed",
      add = TRUE)
plot(dV[-burn], dW[-burn], pch = ".", cex = 1.5, ylab = "")



######################################################################
## Additional univariate and multivariate examples
######################################################################

x <- matrix(c(rep(0, 27), rep(1, length(Nile) - 27)), ncol = 1)
modNileReg <- dlmModReg(x, dW = c(1, 0))
buildFun <- function(theta){
    V(modNileReg) <- exp(theta[1])
    diag(W(modNileReg))[1] <- exp(theta[2])
    modNileReg
}

fit <- dlmMLE(Nile, parm = rep(0, 2), build = buildFun)
modNileReg <- buildFun(fit$par)

modSmooth <- dlmSmooth(Nile, mod = modNileReg)
plot(Nile, type = "o")
lines(ts(modSmooth$s[-1, 1] + modSmooth$s[-1, 2] * x, start = 1871),
      lty = 2, col = "red")


lGas <- log(UKgas)
dlmGas <- dlmModPoly() + dlmModSeas(4)
buildFun <- function(x) {
    diag(W(dlmGas))[2:3] <- exp(x[1:2])
    V(dlmGas) <- exp(x[3])
    dlmGas
}
fit <- dlmMLE(lGas, parm = rep(0, 3), build = buildFun)
dlmGas <- buildFun(fit$par)

gasSmooth <- dlmSmooth(lGas, mod = dlmGas)
x <- cbind(lGas, dropFirst(gasSmooth$s[, c(1, 3)]))
colnames(x) <- c("Gas", "Trend", "seasonal")
plot(x, type = "o", main = "UK Gas Consumption")

######################################################################
## Formulating State Space Models in R with Focus on Logitudinal
## Regression Models
######################################################################

