######################################################################
## Garch modeling from Sumhway and Stoffer
######################################################################

## At Home
load("c:/Users/user/Dropbox/Uni Work/Masters Thesis/Codes/State Space Model/tsa3.rda")

load("c:/Users/user/Dropbox/Uni Work/Masters THesis/Codes/Constraint Newton Method/cnm.R")

## At Work
#load("~/Dropbox/Uni Work/Masters Thesis/Codes/State Space Model/tsa3.rda")

## Have a look at the log gnp data and also fit an AR(1) process to
## remove the mean
gnpgr <- diff(log(gnp))
sarima(gnpgr, 1, 0, 0)
acf2(innov^2, 24)

## Fit the AR(1)-ARCH(1) model to the GNP data
library(fGarch)
gnp.garch <- garchFit(~arma(1, 0) + garch(1, 0), gnpgr)
summary(gnp.garch)

## Fit the GARCH(1, 1) model to the NYSE data
nyse.garch <- garchFit(~garch(1, 1), nyse)
summary(nyse.garch)

u <- nyse.garch@sigma.t
plot(window(nyse, start = 900, end = 1000),
     ylim = c(-.22, .2), ylab = "NYSE Returns")
lines(window(nyse - 2 * u, start = 900, end = 1000), lty = 2, col = "blue")
lines(window(nyse + 2 * u, start = 900, end = 1000), lty = 2, col = "blue")

######################################################################
## fGarch
######################################################################

## garchSpec -
    #Use default parameters beside alpha:
spec <- garchSpec(model = list(alpha = c(0.05, 0.05)))
spec
coef(spec)

## garchSim -
   # Simulate an univariate "timeSeries" series
x <- garchSim(spec, n = 200)

## garchFit -
fit <- garchFit( ~ garch(1, 1), data = x)

#####
## Swiss Pension
x <- as.timeSeries(data(LPP2005REC))

(fit <- garchFit(LPP40 ~ garch(1, 1),
                data = x, trace = FALSE, cond.dist = "norm"))

(fit <- garchFit(LPP40 ~ garch(1, 1),
                data = x, trace = FALSE, cond.dist = "snorm"))



(fit <- garchFit(LPP40 ~ garch(1, 1),
                data = x[-377, ], trace = FALSE))



gll <- double()
for(i in 10:nrow(x)){
    gll[i - 1] <- garchFit(LPP40 ~ garch(1, 1),
                       data = x[1:i, ], trace = FALSE)@fit$llh/
                           garchFit(LPP40 ~ garch(1, 1),
                                    data = x[1:(i - 1), ],
                                    trace = FALSE)@fit$llh
}






plot(x[, 8])
lines(as.timeSeries(fitted(fit)), col = "red", lty = 2)

plot(fit, which = c(1:13))

## sged -
par(mfrow = c(2, 2))
set.seed(1953)
r = rsged(n = 1000)
plot(r, type = "l", main = "sged", col = "steelblue")

# Plot empirical density and compare with true density:
hist(r, n = 25, probability = TRUE, border = "white", col = "steelblue")
box()
x = seq(min(r), max(r), length = 201)
lines(x, dsged(x), lwd = 2)

# Plot df and compare with true df:
plot(sort(r), (1:1000/1000), main = "Probability", col = "steelblue",
ylab = "Probability")
lines(x, psged(x), lwd = 2)

# Compute quantiles:
round(qsged(psged(q = seq(-1, 5, by = 1))), digits = 6)
## sgedFit -
sgedFit(r)

## Not run:
## sgedSlider -
if (require(tcltk)) {
sgedSlider("dist")
sgedSlider("rand")
}
## End(Not run)


