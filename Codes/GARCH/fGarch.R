######################################################################
## This script is used to understand the working of the garch function
## and hor the variables are calculated and initialised.
######################################################################

library(fGarch)


## garchFit

formula = ~garch(1, 1)
data = dem2gbp
init.rec = "mci"
delta = 2
skew = 1
shape = 4
cond.dist = "norm"
include.mean = TRUE
include.delta = NULL
include.skew = NULL
include.shape = NULL
leverage = NULL
trace = TRUE
algorithm = "nlminb"
hessian = "ropt"
control = list()
title = NULL
description = NULL
DEBUG = FALSE

init.rec = match.arg(init.rec)
cond.dist = match.arg(cond.dist)
hessian = match.arg(hessian)
algorithm = match.arg(algorithm)
CALL = match.call()
Name = capture.output(substitute(data))

if (is.character(data)) {
        eval(parse(text = paste("data(", data, ")")))
        data = eval(parse(text = data))
    }


data <- as.data.frame(data)

if (isUnivariate(data)) {
    colnames(data) <- "data"
} else {
        uniqueNames = unique(sort(colnames(data)))
        if (is.null(colnames(data))) {
            stop("Column names of data are missing.")
        }
        if (length(colnames(data)) != length(uniqueNames)) {
            stop("Column names of data are not unique.")
        }
    }

if (length(formula) == 3 && isUnivariate(data))
        formula[2] <- NULL


if (length(formula) == 2) {
        if (isUnivariate(data)) {
            formula = as.formula(paste("data", paste(formula,
                collapse = " ")))
        }
        else {
            stop("Multivariate data inputs require lhs for the formula.")
        }
    }

robust.cvar <- (cond.dist == "QMLE")

args = .garchArgsParser(formula = formula, data = data, trace = FALSE)

ans <- .garchFit(formula.mean = args$formula.mean,
        formula.var = args$formula.var,
        series = args$series, init.rec, delta, skew, shape, cond.dist,
        include.mean, include.delta, include.skew, include.shape,
        leverage, trace, algorithm, hessian, robust.cvar, control,
        title, description)

    ans@call = CALL
    attr(formula, "data") <- paste("data = ", Name, sep = "")
    ans@formula = formula
    ans



## .garchFit

formula.mean = args$formula.mean
formula.var = args$formula.var
series = args$series
DEBUG <- FALSE
.StartFit <- Sys.time()

con <- .garchOptimizerControl(algorithm, cond.dist)
con[(namc <- names(control))] <- control
data <- series
scale <- if (con$xscale) sd(series) else 1
series <- series/scale


.series <- .garchInitSeries(formula.mean = formula.mean,
        formula.var = formula.var, cond.dist = cond.dist[1],
        series = series, scale = scale, init.rec = init.rec[1],
        h.start = NULL, llh.start = NULL, trace = trace)

.setfGarchEnv(.series = .series)


.params <- .garchInitParameters(formula.mean = formula.mean,
        formula.var = formula.var, delta = delta, skew = skew,
        shape = shape, cond.dist = cond.dist[1], include.mean = include.mean,
        include.delta = include.delta, include.skew = include.skew,
        include.shape = include.shape, leverage = leverage,
        algorithm = algorithm[1],
        control = con, trace = trace)

.setfGarchEnv(.params = .params)
.setfGarchEnv(.garchDist = .garchSetCondDist(cond.dist[1]))

.setfGarchEnv(.llh = 1e+99)
.llh <- .getfGarchEnv(".llh")

## This part changes the .series
fit = .garchOptimizeLLH(hessian, robust.cvar, trace)

.series <- .getfGarchEnv(".series")
.params <- .getfGarchEnv(".params")

names(.series$h) <- NULL

fit$series = .series
fit$params = .params


residuals = .series$z
fitted.values = .series$x - residuals
h.t = .series$h

if (.params$includes["delta"])
        deltainv = 1/fit$par["delta"] else deltainv = 1/fit$params$delta

sigma.t = (.series$h)^deltainv

fit$cvar <- if (robust.cvar){
               (solve(fit$hessian) %*%
                (t(fit$gradient) %*% fit$gradient) %*%
                solve(fit$hessian))} else { -solve(fit$hessian)}


fit$se.coef = sqrt(diag(fit$cvar))

fit$tval = fit$coef/fit$se.coef

## The summary
fit$matcoef = cbind(fit$coef, fit$se.coef, fit$tval,
                    2 *(1 - pnorm(abs(fit$tval))))

dimnames(fit$matcoef) = list(names(fit$tval), c(" Estimate",
        " Std. Error", " t value", "Pr(>|t|)"))


if (is.null(title)) title = "GARCH Modelling"
if (is.null(description)) description = description()

Time = Sys.time() - .StartFit

## .garchInitSeries
h.start = NULL
llh.start = NULL
mm = length(formula.mean)

end = regexpr("\\(", as.character(formula.mean[mm])) - 1

model.mean = substr(as.character(formula.mean[mm]), 1, end)

mv = length(formula.var)
if (mv != 2)
        stop("Variance Formula misspecified")
end = regexpr("\\(", as.character(formula.var[mv])) - 1
model.var = substr(as.character(formula.var[mv]), 1, end)
if (!any(c("garch", "aparch") == model.var))
        stop("formula.var must be one of: garch, aparch")


model.order = as.numeric(strsplit(strsplit(strsplit(as.character(formula.mean),
"\\(")[[2]][2], "\\)")[[1]], ",")[[1]])

u = model.order[1]
v = 0

if (length(model.order) == 2) v = model.order[2]

maxuv = max(u, v)

model.order = as.numeric(strsplit(strsplit(strsplit(as.character(formula.var),
        "\\(")[[2]][2], "\\)")[[1]], ",")[[1]])

p = model.order[1]
q = 0
if (length(model.order) == 2) q = model.order[2]
maxpq = max(p, q)
max.order = max(maxuv, maxpq)

if (is.null(h.start))
        h.start = max.order + 1

if (is.null(llh.start))
        llh.start = 1

if (init.rec != "mci" & model.var != "garch") {
        stop("Algorithm only supported for mci Recursion")
}

## .garchInitParameters

.DEBUG = FALSE
.series <- .getfGarchEnv(".series")

model.order = as.numeric(strsplit(strsplit(strsplit(as.character(formula.mean),
        "\\(")[[2]][2], "\\)")[[1]], ",")[[1]])
u = model.order[1]
v = 0

if (length(model.order) == 2)
    v = model.order[2]

model.order = as.numeric(strsplit(strsplit(strsplit(as.character(formula.var),
        "\\(")[[2]][2], "\\)")[[1]], ",")[[1]])
p = model.order[1]
q = 0

if (length(model.order) == 2) q = model.order[2]

model.var = .series$model[2]

if (is.null(include.delta)) {
if (model.var == "garch") {
    include.delta = FALSE
} else {
            include.delta = TRUE
        }
}

if (is.null(leverage)) {
        if (model.var == "garch") {
            leverage = FALSE
        }
        else {
            leverage = TRUE
        }
    }

skewed.dists = c("snorm", "sged", "sstd", "snig")
if (is.null(include.skew)) {
        if (any(skewed.dists == cond.dist)) {
            include.skew = TRUE
        }
        else {
            include.skew = FALSE
        }
}

shaped.dists = c("ged", "sged", "std", "sstd", "snig")
    if (is.null(include.shape)) {
        if (any(shaped.dists == cond.dist)) {
            include.shape = TRUE
        }
        else {
            include.shape = FALSE
        }
    }

Names = c("mu", if (u > 0) paste("ar", 1:u, sep = ""), if (v >
        0) paste("ma", 1:v, sep = ""), "omega", if (p > 0) paste("alpha",
        1:p, sep = ""), if (p > 0) paste("gamma", 1:p, sep = ""),
        if (q > 0) paste("beta", 1:q, sep = ""), "delta", "skew",
        "shape")

## The mean were fitted via the ARIMA model first
fit.mean = arima(.series$x, order = c(u, 0, v),
                 include.mean = include.mean)$coef


alpha.start = 0.1
beta.start = 0.8
params = c(if (include.mean) fit.mean[length(fit.mean)] else 0,
        if (u > 0) fit.mean[1:u], if (v > 0) fit.mean[(u + 1):(length(fit.mean) -
            as.integer(include.mean))], var(.series$x, na.rm = TRUE) *
            (1 - alpha.start - beta.start), if (p > 0) rep(alpha.start/p,
            times = p), if (p > 0) rep(0.1, times = p), if (q >
            0) rep(beta.start/q, times = q), delta, skew, shape)
names(params) = Names


TINY = 1e-08
USKEW = 1/10
USHAPE = 1

if (cond.dist == "snig")
        USKEW = -0.99
    U = c(-10 * abs(mean(.series$x)), if (u > 0) rep(-1 + TINY,
        times = u), if (v > 0) rep(-1 + TINY, times = v), 1e-06 *
        var(.series$x), if (p > 0) rep(0 + TINY, times = p),
        if (p > 0) rep(-1 + TINY, times = p), if (q > 0) rep(0 +
            TINY, times = q), 0, USKEW, USHAPE)
names(U) = Names


VSKEW = 10
VSHAPE = 10

if (cond.dist == "snig")
        VSKEW = 0.99
    V = c(10 * abs(mean(.series$x)), if (u > 0) rep(1 - TINY,
        times = u), if (v > 0) rep(1 - TINY, times = v), 100 *
        var(.series$x), if (p > 0) rep(1 - TINY, times = p),
        if (p > 0) rep(1 - TINY, times = p), if (q > 0) rep(1 -
            TINY, times = q), 2, VSKEW, VSHAPE)
names(V) = Names



includes = c(include.mean, if (u > 0) rep(TRUE, times = u),
        if (v > 0) rep(TRUE, times = v), TRUE, if (p > 0) rep(TRUE,
            times = p), if (p > 0) rep(leverage, times = p),
        if (q > 0) rep(TRUE, times = q), include.delta, include.skew,
        include.shape)
names(includes) = Names




index = (1:length(params))[includes == TRUE]
names(index) = names(params)[includes == TRUE]

alpha <- beta <- NULL

if (p > 0)
        alpha = params[substr(Names, 1, 5) == "alpha"]

if (p > 0 & leverage)
        gamma = params[substr(Names, 1, 5) == "gamma"]

if (p > 0 & !leverage)
        gamma = rep(0, times = p)

if (q > 0)
        beta = params[substr(Names, 1, 4) == "beta"]

if (.series$model[2] == "garch") {
        persistence = sum(alpha) + sum(beta)
} else if (.series$model[2] == "aparch") {
        persistence = sum(beta)
        for (i in 1:p) persistence = persistence + alpha[i] *
            garchKappa(cond.dist, gamma[i], params["delta"],
                params["skew"], params["shape"])
    }
names(persistence) = "persistence"

ans = data.frame(U, V, params, includes)
rownames(ans) = paste("   ", names(params))


## Not sure what U and V really are. Will need to check this.


## .garchOptimizeLLH
DEBUG = FALSE
.series <- .getfGarchEnv(".series")
.params <- .getfGarchEnv(".params")
INDEX = .params$index
algorithm = .params$control$algorithm[1]
TOL1 = .params$control$tol1
TOL2 = .params$control$tol2

if (algorithm == "nlminb" | algorithm == "nlminb+nm") {
        fit <- .garchRnlminb(.params, .series, .garchLLH, trace)
        .params$llh = fit$llh
        .params$params[INDEX] = fit$par
        .setfGarchEnv(.params = .params)
    }

.params$llh = fit$llh
.params$params[INDEX] = fit$par
.setfGarchEnv(.params = .params)

if (hessian == "ropt") {
        fit$hessian <- -.garchRoptimhess(par = fit$par, .params = .params,
            .series = .series)
        titleHessian = "R-optimhess"
    }

if (.params$control$xscale) {
    .series$x <- .series$x * .series$scale
        if (.params$include["mu"])
            fit$coef["mu"] <- fit$par["mu"] <- .params$params["mu"] <- .params$params["mu"] * .series$scale
        if (.params$include["omega"])
            fit$coef["omega"] <- fit$par["omega"] <- .params$params["omega"] <- .params$params["omega"] * .series$scale^(.params$params["delta"])
        .setfGarchEnv(.params = .params)
        .setfGarchEnv(.series = .series)
    }

if (.params$control$xscale) {
        if (.params$include["mu"]) {
            fit$hessian[, "mu"] <- fit$hessian[, "mu"]/.series$scale
            fit$hessian["mu", ] <- fit$hessian["mu", ]/.series$scale
        }
        if (.params$include["omega"]) {
            fit$hessian[, "omega"] <- fit$hessian[, "omega"]/.series$scale^(.params$params["delta"])
            fit$hessian["omega", ] <- fit$hessian["omega", ]/.series$scale^(.params$params["delta"])
        }
    }

.llh <- fit$llh <- fit$value <- .garchLLH(fit$par, trace = FALSE,
        fGarchEnv = TRUE)

.series <- .getfGarchEnv(".series")

if (robust.cvar)
        fit$gradient <- -.garchRCDAGradient(par = fit$par, .params = .params,
            .series = .series)

N = length(.series$x)
NPAR = length(fit$par)

fit$ics = c(AIC = c((2 * fit$value)/N + 2 * NPAR/N), BIC = (2 *
        fit$value)/N + NPAR * log(N)/N, SIC = (2 * fit$value)/N +
        log((N + 2 * NPAR)/N), HQIC = (2 * fit$value)/N + (2 *
        NPAR * log(log(N)))/N)

names(fit$ics) <- c("AIC", "BIC", "SIC", "HQIC")
fit
