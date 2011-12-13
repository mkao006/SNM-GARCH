######################################################################
## Reproduceable example
######################################################################

library(compiler)
library(tseries)
library(fGarch)
source("~/dropbox/Uni Work/Masters Thesis/Codes/Constraint Newton Method/cnm.R")
source("~/dropbox/Uni Work/Masters Thesis/Codes/Constraint Newton Method/dden.R")

#source("functions.R")
#source("functionMeanSkew.R")
souce("functionMean.R")
source("mgarchsim.R")




######################################################################
## DEM2GBP
######################################################################

## Example for the model
tdg <- as.numeric(data.matrix(dem2gbp))
tmdg <- tdg
class(tmdg) <- "mgarch"
dg.t0 <- garchFit(data = tdg, cond.dist = "norm")
dg.t1 <- garchFit(data = tdg, cond.dist = "ged")
dg.t2 <- garchFit(data = tdg, cond.dist = "std")
(dg.mg <- cnmms(tmdg, plot = "gradient", grid = 100, verb = 4))


dg.t0@fit$llh
dg.t1@fit$llh
dg.t2@fit$llh
dg.mg$ll

## Now compute the sigma.t
dg.sigmat <- c(length(tdg))
dg.sigmat[1] <- dg.mg$beta[4]
for(i in 2:length(tdg)){
    dg.sigmat[i] <- dg.mg$beta[1] +
        dg.mg$beta[2] * dg.sigmat[i - 1] +
            dg.mg$beta[3] * tdg[i - 1]^2
}
dg.sigmat <- sqrt(dg.sigmat)

par(mfrow = c(2, 2), mar = c(2.1, 4.1, 4.1, 1.1))
hist((tdg - coef(dg.t0)[1])/dg.t0@sigma.t, freq = FALSE,
     breaks = 200, xlim = c(-4, 4), ylim = c(0, 0.6), main = "")
curve(dnorm(x, 0, 1), add = TRUE, col = "blue", lwd = 3)
lines(density((tdg - coef(dg.t0)[1])/dg.t0@sigma.t), col = "red", lwd = 3)
text(-4, 0.6, "Standard Normal", adj = 0, cex = 1.5)
box()
legend("topright", legend = c("Density", "Fitted"),
       col = c("red", "blue"), bty = "n", lty = 1, lwd = 3)
par(mar = c(2.1, 3.1, 4.1, 2.1))
hist((tdg - coef(dg.t2)[1])/dg.t2@sigma.t, freq = FALSE,
     breaks = 200, xlim = c(-4, 4),
     ylim = c(0, 0.6), main = "", xlab = "Error Distribution")
curve(dstd(x, 0, 1, coef(dg.t2)[5]), add = TRUE, col = "blue", lwd = 3)
lines(density((tdg - coef(dg.t2)[1])/dg.t2@sigma.t), col = "red", lwd = 3)
text(-4, 0.6, paste("t (", round(coef(dg.t2)[5], 2), ")", sep = ""),
     adj = 0, cex = 1.5)
box()
par(mar = c(5.1, 4.1, 1.1, 1.1))
hist((tdg - coef(dg.t1)[1])/dg.t1@sigma.t, freq = FALSE,
     breaks = 200, xlim = c(-4, 4),
     ylim = c(0, 0.6), main = "", ylab = "", xlab = "")
curve(dged(x, 0, 1, coef(dg.t1)[5]), add = TRUE, col = "blue", lwd = 3)
lines(density((tdg - coef(dg.t1)[1])/dg.t1@sigma.t), col = "red", lwd = 3)
text(-4, 0.6, paste("Generalised Error Distribution (", round(coef(dg.t1)[5], 2), ")", sep = ""),
     adj = 0, cex = 1.5)
box()
par(mar = c(5.1, 3.1, 1.1, 2.1))
hist((tdg - coef(dg.t0)[1])/dg.sigmat, freq = FALSE,
     breaks = 200, xlim = c(-4, 4), ylim = c(0, 0.6), main = "",
     xlab = "Error Distribution")
#curve(sum(dg.mg$mix$pr * dnorm(x, 0, dg.mg$mix$pt)), add = TRUE, col = "blue", lwd = 3)
curve(dg.mg$mix$pr[1] * dnorm(x, 0, dg.mg$mix$pt[1]) +
      dg.mg$mix$pr[2] * dnorm(x, 0, dg.mg$mix$pt[2]) +
      dg.mg$mix$pr[3] * dnorm(x, 0, dg.mg$mix$pt[3]) +
      dg.mg$mix$pr[4] * dnorm(x, 0, dg.mg$mix$pt[4]) +
      dg.mg$mix$pr[5] * dnorm(x, 0, dg.mg$mix$pt[5]) +
      dg.mg$mix$pr[6] * dnorm(x, 0, dg.mg$mix$pt[6]),
      add = TRUE, col = "blue", lwd = 3)
curve(dg.mg$mix$pr[1] * dnorm(x, 0, dg.mg$mix$pt[1]),
      add = TRUE, col = "light blue", lwd = 2)
curve(dg.mg$mix$pr[2] * dnorm(x, 0, dg.mg$mix$pt[2]),
      add = TRUE, col = "light blue", lwd = 2)
curve(dg.mg$mix$pr[3] * dnorm(x, 0, dg.mg$mix$pt[3]),
      add = TRUE, col = "light blue", lwd = 2)
curve(dg.mg$mix$pr[4] * dnorm(x, 0, dg.mg$mix$pt[4]),
      add = TRUE, col = "light blue", lwd = 2)
curve(dg.mg$mix$pr[5] * dnorm(x, 0, dg.mg$mix$pt[5]),
      add = TRUE, col = "light blue", lwd = 2)
curve(dg.mg$mix$pr[6] * dnorm(x, 0, dg.mg$mix$pt[6]),
      add = TRUE, col = "light blue", lwd = 2)
lines(density((tdg - coef(dg.t0)[1])/dg.sigmat), col = "red", lwd = 3)
text(-4, 0.6, "Scale Normal Mixture", adj = 0, cex = 1.5)
box()

######################################################################
## S&P example
######################################################################

sp <- get.hist.quote("SPY", quote = "Open", start = "2005-01-01",
                      provider = "yahoo")


tsp <- diff(log(as.numeric(sp)))[1:700]
tmsp <- tsp
class(tmsp) <- "mgarch"
sp.t0 <- garchFit(data = tsp, cond.dist = "norm")
sp.t1 <- garchFit(data = tsp, cond.dist = "ged")
sp.t2 <- garchFit(data = tsp, cond.dist = "std")
(sp.mg <- cnmms(tmsp, plot = "gradient", grid = 1000))


sp.t0@fit$llh
sp.t1@fit$llh
sp.t2@fit$llh
sp.mg$ll

## Now compute the sigma.t
sp.sigmat <- c(length(tsp))
sp.sigmat[1] <- sp.mg$beta[4]
for(i in 2:length(tsp)){
    sp.sigmat[i] <- sp.mg$beta[1] +
        sp.mg$beta[2] * sp.sigmat[i - 1] +
            sp.mg$beta[3] * tsp[i - 1]^2
}
sp.sigmat <- sqrt(sp.sigmat)

par(mfrow = c(2, 2), mar = c(2.1, 4.1, 4.1, 1.1))
hist((tsp - coef(sp.t0)[1])/sp.t0@sigma.t, freq = FALSE,
     breaks = 200, xlim = c(-4, 4), ylim = c(0, 0.6), main = "")
curve(dnorm(x, 0, 1), add = TRUE, col = "blue", lwd = 3)
lines(density((tsp - coef(sp.t0)[1])/sp.t0@sigma.t), col = "red", lwd = 3)
text(-4, 0.6, "Standard Normal", adj = 0, cex = 1.5)
box()
legend("topright", legend = c("Density", "Fitted"),
       col = c("red", "blue"), bty = "n", lty = 1, lwd = 3)
par(mar = c(2.1, 3.1, 4.1, 2.1))
hist((tsp - coef(sp.t2)[1])/sp.t2@sigma.t, freq = FALSE,
     breaks = 200, xlim = c(-4, 4),
     ylim = c(0, 0.6), main = "", xlab = "Error Distribution")
curve(dstd(x, 0, 1, coef(sp.t2)[5]), add = TRUE, col = "blue", lwd = 3)
lines(density((tsp - coef(sp.t2)[1])/sp.t2@sigma.t), col = "red", lwd = 3)
text(-4, 0.6, paste("t (", round(coef(sp.t2)[5], 2), ")", sep = ""),
     adj = 0, cex = 1.5)
box()
par(mar = c(5.1, 4.1, 1.1, 1.1))
hist((tsp - coef(sp.t1)[1])/sp.t1@sigma.t, freq = FALSE,
     breaks = 200, xlim = c(-4, 4),
     ylim = c(0, 0.6), main = "", ylab = "", xlab = "")
curve(dged(x, 0, 1, coef(sp.t1)[5]), add = TRUE, col = "blue", lwd = 3)
lines(density((tsp - coef(sp.t1)[1])/sp.t1@sigma.t), col = "red", lwd = 3)
text(-4, 0.6, paste("Generalised Error Distribution (", round(coef(sp.t1)[5], 2), ")", sep = ""),
     adj = 0, cex = 1.5)
box()
par(mar = c(5.1, 3.1, 1.1, 2.1))
hist((tsp - coef(sp.t0)[1])/sp.sigmat, freq = FALSE,
     breaks = 200, xlim = c(-4, 4), ylim = c(0, 0.6), main = "",
     xlab = "Error Distribution")
#curve(sum(sp.mg$mix$pr * dnorm(x, 0, sp.mg$mix$pt)), add = TRUE, col = "blue", lwd = 3)
curve(sp.mg$mix$pr[1] * dnorm(x, 0, sp.mg$mix$pt[1]) +
      sp.mg$mix$pr[2] * dnorm(x, 0, sp.mg$mix$pt[2]),
      add = TRUE, col = "blue", lwd = 3)
curve(sp.mg$mix$pr[1] * dnorm(x, 0, sp.mg$mix$pt[1]),
      add = TRUE, col = "light blue", lwd = 2)
curve(sp.mg$mix$pr[2] * dnorm(x, 0, sp.mg$mix$pt[2]),
      add = TRUE, col = "light blue", lwd = 2)
lines(density((tsp - coef(sp.t0)[1])/sp.sigmat), col = "red", lwd = 3)
text(-4, 0.6, "Scale Normal Mixture", adj = 0, cex = 1.5)
box()

######################################################################
## NYSE example
######################################################################

load("~/Dropbox/Uni Work/Masters Thesis/Codes/State Space Model/tsa3.rda")

tnyse <- as.numeric(nyse)[1:1000]
tmnyse <- tnyse
class(tmnyse) <- "mgarch"
(nyse.mg <- cnmms(tmnyse, plot = "gradient", grid = 300, verb = 4))
nyse.t0 <- garchFit(data = tnyse, cond.dist = "norm")
nyse.t1 <- garchFit(data = tnyse, cond.dist = "ged")
nyse.t2 <- garchFit(data = tnyse, cond.dist = "std")


nyse.t0@fit$llh
nyse.t1@fit$llh
nyse.t2@fit$llh
nyse.mg$ll

## Now compute the sigma.t
nyse.sigmat <- c(length(tnyse))
nyse.sigmat[1] <- nyse.mg$beta[4]
for(i in 2:length(tnyse)){
    nyse.sigmat[i] <- nyse.mg$beta[1] +
        nyse.mg$beta[2] * nyse.sigmat[i - 1] +
            nyse.mg$beta[3] * tnyse[i - 1]^2
}
nyse.sigmat <- sqrt(nyse.sigmat)

par(mfrow = c(2, 2), mar = c(2.1, 4.1, 4.1, 1.1))
hist((tnyse - coef(nyse.t0)[1])/nyse.t0@sigma.t, freq = FALSE,
     breaks = 200, xlim = c(-4, 4), ylim = c(0, 0.6), main = "")
curve(dnorm(x, 0, 1), add = TRUE, col = "blue", lwd = 3)
lines(density((tnyse - coef(nyse.t0)[1])/nyse.t0@sigma.t), col = "red", lwd = 3)
text(-4, 0.6, "Standard Normal", adj = 0, cex = 1.5)
box()
legend("topright", legend = c("Density", "Fitted"),
       col = c("red", "blue"), bty = "n", lty = 1, lwd = 3)
par(mar = c(2.1, 3.1, 4.1, 2.1))
hist((tnyse - coef(nyse.t2)[1])/nyse.t2@sigma.t, freq = FALSE,
     breaks = 200, xlim = c(-4, 4),
     ylim = c(0, 0.6), main = "", xlab = "Error Distribution")
curve(dstd(x, 0, 1, coef(nyse.t2)[5]), add = TRUE, col = "blue", lwd = 3)
lines(density((tnyse - coef(nyse.t2)[1])/nyse.t2@sigma.t), col = "red", lwd = 3)
text(-4, 0.6, paste("t (", round(coef(nyse.t2)[5], 2), ")", sep = ""),
     adj = 0, cex = 1.5)
box()
par(mar = c(5.1, 4.1, 1.1, 1.1))
hist((tnyse - coef(nyse.t1)[1])/nyse.t1@sigma.t, freq = FALSE,
     breaks = 200, xlim = c(-4, 4),
     ylim = c(0, 0.6), main = "", ylab = "", xlab = "")
curve(dged(x, 0, 1, coef(nyse.t1)[5]), add = TRUE, col = "blue", lwd = 3)
lines(density((tnyse - coef(nyse.t1)[1])/nyse.t1@sigma.t), col = "red", lwd = 3)
text(-4, 0.6, paste("Generalised Error Distribution (", round(coef(nyse.t1)[5], 2), ")", sep = ""),
     adj = 0, cex = 1.5)
box()
par(mar = c(5.1, 3.1, 1.1, 2.1))
hist((tnyse - coef(nyse.t0)[1])/nyse.sigmat, freq = FALSE,
     breaks = 200, xlim = c(-4, 4), ylim = c(0, 0.6), main = "",
     xlab = "Error Distribution")
#curve(sum(nyse.mg$mix$pr * dnorm(x, 0, nyse.mg$mix$pt)), add = TRUE, col = "blue", lwd = 3)
curve(nyse.mg$mix$pr[1] * dnorm(x, 0, nyse.mg$mix$pt[1]) +
      nyse.mg$mix$pr[2] * dnorm(x, 0, nyse.mg$mix$pt[2]) +
      nyse.mg$mix$pr[3] * dnorm(x, 0, nyse.mg$mix$pt[3]) +
      nyse.mg$mix$pr[4] * dnorm(x, 0, nyse.mg$mix$pt[4]) +
      nyse.mg$mix$pr[5] * dnorm(x, 0, nyse.mg$mix$pt[5]),
      add = TRUE, col = "blue", lwd = 3)
curve(nyse.mg$mix$pr[1] * dnorm(x, 0, nyse.mg$mix$pt[1]),
      add = TRUE, col = "light blue", lwd = 2)
curve(nyse.mg$mix$pr[2] * dnorm(x, 0, nyse.mg$mix$pt[2]),
      add = TRUE, col = "light blue", lwd = 2)
curve(nyse.mg$mix$pr[3] * dnorm(x, 0, nyse.mg$mix$pt[3]),
      add = TRUE, col = "light blue", lwd = 2)
curve(nyse.mg$mix$pr[4] * dnorm(x, 0, nyse.mg$mix$pt[4]),
      add = TRUE, col = "light blue", lwd = 2)
curve(nyse.mg$mix$pr[5] * dnorm(x, 0, nyse.mg$mix$pt[5]),
      add = TRUE, col = "light blue", lwd = 2)
lines(density((tnyse - coef(nyse.t0)[1])/nyse.sigmat), col = "red", lwd = 3)
text(-4, 0.6, "Scale Normal Mixture", adj = 0, cex = 1.5)
box()

