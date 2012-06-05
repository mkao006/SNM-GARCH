########################################################################
## This scripts extends the functions in the fGarch package to
## incorporate our study
########################################################################


garchSpec <- function(model = list(), presample = NULL,
                      cond.dist = c("norm", "ged", "std", "cauchy", "snorm",
                          "sged", "smnorm", "sstd", "ssmnorm"),
                      rseed = NULL){
  cond.dist = match.arg(cond.dist)
  skew = list(norm = NULL, ged = NULL, std = NULL, smnorm = NULL,
    snorm = 0.9, sged = 0.9, sstd = 0.9, ssmnorm = 0.9)
  shape = list(norm = NULL, ged = 2, std = 4, smnorm = disc(1, 1),
    sged = 2, sstd = 4, ssmnorm = disc(1, 1))
  control = list(omega = 1e-06, alpha = 0.1, gamma = NULL,
    beta = 0.8, mu = NULL, ar = NULL, ma = NULL, delta = 2,
    skew = skew[[cond.dist]], shape = shape[[cond.dist]])
  control[names(model)] <- model
  model <- control
  if (sum(c(model$alpha, model$beta)) > 1)
    warnings("sum(alpha)+sum(beta)>1")
  order.ar = length(model$ar)
  order.ma = length(model$ma)
  order.alpha = length(model$alpha)
  if (sum(model$beta) == 0) {
    order.beta = 0
  }
  else {
    order.beta = length(model$beta)
  }
  if (order.ar == 0 && order.ma == 0) {
    formula.mean = ""
  }
  if (order.ar > 0 && order.ma == 0) {
    formula.mean = paste("ar(", as.character(order.ar), ")",
      sep = "")
  }
  if (order.ar == 0 && order.ma > 0) {
    formula.mean = paste("ma(", as.character(order.ma), ")",
      sep = "")
  }
  if (order.ar > 0 && order.ma > 0) {
    formula.mean = paste("arma(", as.character(order.ar),
      ", ", as.character(order.ma), ")", sep = "")
  }
  formula.var = "garch"
  if (order.beta == 0)
    formula.var = "arch"
  if (!is.null(model$gamma) != 0)
    formula.var = "aparch"
  if (model$delta != 2)
    formula.var = "aparch"
  if (order.beta == 0) {
    formula.var = paste(formula.var, "(", as.character(order.alpha),
      ")", sep = "")
  }
  else {
    formula.var = paste(formula.var, "(", as.character(order.alpha),
      ", ", as.character(order.beta), ")", sep = "")
  }
  if (formula.mean == "") {
    formula = as.formula(paste("~", formula.var))
  }
  else {
    formula = as.formula(paste("~", formula.mean, "+", formula.var))
  }
  if (is.null(model$mu))
    model$mu = 0
  if (is.null(model$ar))
    model$ar = 0
  if (is.null(model$ma))
    model$ma = 0
  if (is.null(model$gamma))
    model$gamma = rep(0, times = order.alpha)
  if (is.null(rseed)) {
    rseed = 0
  }
  else {
    set.seed(rseed)
  }
  order.max = max(order.ar, order.ma, order.alpha, order.beta)
  iterate = TRUE
  if (!is.matrix(presample)) {
    if (is.null(presample)) {
      iterate = FALSE
      n.start = order.max
    }
    else {
      n.start = presample
    }
    z = rnorm(n = n.start)
    h = rep(model$omega/(1 - sum(model$alpha) - sum(model$beta)),
      times = n.start)
    y = rep(model$mu/(1 - sum(model$ar)), times = n.start)
  }
  else {
    z = presample[, 1]
    h = presample[, 2]
    y = presample[, 3]
  }
  presample = cbind(z, h, y)
  if (iterate) {
    n.iterate = length(z) - order.max
    deltainv = 1/model$delta
    for (i in n.iterate:1) {
      h[i] = model$omega + sum(model$alpha * (abs(abs(y[i +
         (1:order.alpha)]) - model$gamma *
         y[i + (1:order.alpha)])^model$delta)) +
           sum(model$beta * h[i + (1:order.beta)])
      y[i] = model$mu + sum(model$ar * y[i + (1:order.ar)]) +
        sum(model$ma * (h[i + (1:order.ma)]^deltainv)) +
          h[i]^deltainv * z[i]
    }
  }
  new("fGARCHSPEC", call = match.call(), formula = formula,
      model = list(omega = model$omega, alpha = model$alpha,
        gamma = model$gamma, beta = model$beta, mu = model$mu,
        ar = model$ar, ma = model$ma, delta = model$delta,
        skew = model$skew, shape = model$shape),
      presample = as.matrix(presample),
      distribution = as.character(cond.dist), rseed = as.numeric(rseed))
}

garchSim <- function (spec = garchSpec(), n = 100, n.start = 100,
                      extended = FALSE){
  stopifnot(class(spec) == "fGARCHSPEC")
  model = spec@model
  if (spec@rseed != 0)
    set.seed(spec@rseed)
  n = n + n.start
  if (spec@distribution == "norm")
    z = rnorm(n)
  if (spec@distribution == "ged")
    z = rged(n, nu = model$shape)
  if (spec@distribution == "std")
    z = rstd(n, nu = model$shape)
  if (spec@distribution == "smnorm")
    z = rsmnorm(n, mix = model$shape)
  if (spec@distribution == "snorm")
    z = rsnorm(n, xi = model$skew)
  if (spec@distribution == "sged")
    z = rsged(n, nu = model$shape, xi = model$skew)
  if (spec@distribution == "sstd")
    z = rsstd(n, nu = model$shape, xi = model$skew)
  if (spec@distribution == "ssmnorm")
    z = rssmnorm(n, mix = model$shape, xi = model$skew)
  delta = model$delta
  z = c(rev(spec@presample[, 1]), z)
  h = c(rev(spec@presample[, 2]), rep(NA, times = n))
  y = c(rev(spec@presample[, 3]), rep(NA, times = n))
  m = length(spec@presample[, 1])
  names(z) = names(h) = names(y) = NULL
  mu = model$mu
  ar = model$ar
  ma = model$ma
  omega = model$omega
  alpha = model$alpha
  gamma = model$gamma
  beta = model$beta
  deltainv = 1/delta
  order.ar = length(ar)
  order.ma = length(ma)
  order.alpha = length(alpha)
  order.beta = length(beta)
  eps = h^deltainv * z
  for (i in (m + 1):(n + m)) {
    h[i] = omega + sum(alpha * (abs(eps[i - (1:order.alpha)]) -
       gamma * (eps[i - (1:order.alpha)]))^delta) +
         sum(beta * h[i - (1:order.beta)])
    eps[i] = h[i]^deltainv * z[i]
    y[i] = mu + sum(ar * y[i - (1:order.ar)]) + sum(ma *
       eps[i - (1:order.ma)]) + eps[i]
  }
  data = cbind(z = z[(m + 1):(n + m)], sigma = h[(m + 1):(n +
                                         m)]^deltainv, y = y[(m + 1):(n + m)])
  rownames(data) = as.character(1:n)
  data = data[-(1:n.start), ]
  from <- timeDate(format(Sys.time(), format = "%Y-%m-%d")) -
    NROW(data) * 24 * 3600
  charvec <- timeSequence(from = from, length.out = NROW(data))
  ans <- timeSeries(data = data[, c(3, 2, 1)], charvec = charvec)
  colnames(ans) <- c("garch", "sigma", "eps")
  ans <- if (extended)
    ans
  else ans[, "garch"]
  attr(ans, "control") <- list(garchSpec = spec)
  ans
}
