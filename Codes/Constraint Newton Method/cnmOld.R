# =================================================================== #
# This file contains the constrained Newton method for computing the  #
# nonparametric maximum likelihood estimate of a mixing distribution, #
# where the component distributions are normal or Poisson             #
#                                                                     #
# Author: Yong Wang (yongwang@stat.auckland.ac.nz)                    #
#         Department of Statistics, University of Auckland            #
#         New Zealand                                                 #
# =================================================================== #

# ============================= #
# Functions for normal mixtures #
# ============================= #

# ------------------------------- #
# Constrained Newton method (CNM) #
# ------------------------------- #

# x         Observations
# verb      Verbosity level (0-4)
# prec      Precision level to combine components
# maxit     Maximum number of iterations
# check.points   Number of grid points for computing the derivatives of
#           the gradient function
# tol       Tolerance level for sup gradient
# mix       Initial estimate of the NPMLE
# method    Either CNM or CN1
# lsm       Line search method, either halving (step-halving) or optim (optimal)

cnm.normmix = function(x, verb=0, prec=1e-6,
    maxit=1000, check.points=100, tol=1e-6, mix=normmix(NA,va=1),
    method=c("cnm","cn1"), lsm=c("halving","optim") ) {

  method = match.arg(method)
  s = sqrt(mix$va[1]) ## The standard deviation
  if( length(x) < 2 )
    return( list(mix=normmix(x,1,s^2), ll=dnorm(x,x,sd=s,log=TRUE),
                 num.iterations=0, max.gradient=0, convergence=0) )
  lsm = match.arg(lsm)

  if( length(check.points) == 1 ) {
    check.points = max(check.points, range(x) / s * 20)
    check.points = seq(min(x)-.05*diff(range(x)),
      max(x)+.05*diff(range(x)), length=check.points)
  }
  if( any(is.na(mix$mu)) ) {
    k = max(ceiling( (max(x) - min(x)) / (5*sqrt(mix$va[1])) ), 10)
    ## why is only the first variance being taken? so minimum of 10
    ## starting points?
    mu = seq(min(x), max(x), len=k)
    ## Initial support point of mu?
    counts = hist(x, mu, plot=FALSE)$counts
    mix = normmix(mu = c(mu[1], mu[-1][counts != 0]), va=mix$va)
  }
  n = length(x)

  ll.mix = logLik(mix, x)
  convergence = 1
  for( i in 1:maxit ) {
    mix1 = mix
    ll.mix1 = ll.mix
    g = maxima.gradient(mix, x, check.points=check.points, tol=-Inf)
    gradient = max(g$gradient)
    if( gradient < tol ) { convergence = 0; break }
    if(method == "cnm") nsp = g$mu
    else nsp = g$mu[which.max(g$gradient)]
    mix = normmix(c(mix$mu,nsp), c(mix$pi,rep(0,length(nsp))),
      mix$va[1])
    r = cgn.pi(mix, x, ll.mix=ll.mix, lsm=lsm, gamma=1e-3)
    mix = r$mix
   ll.mix = r$ll
    mix = unique(mix, prec=c(.01,0,.1))
    switch(verb,
           {cat("mu =", mix$mu, "\n")},
           print(mix),
           cat("log-likelihood =", logLik(mix, x), "\n"),
           cat("Maximum gradient = ", gradient, "\n"))
  }
  mix = unique(sort.normmix(mix), prec=prec)
  g = maxima.gradient(mix, x, check.points=check.points, tol=-Inf)
  list(mix=mix, ll = logLik(mix, x), num.iterations=i,
       max.gradient=max(g$gradient), convergence=convergence)
}

# Updates mix$pi by using the constrained Newton method #

# mix       Current mixing distribution
# x         Observations
# maxit     Maximum number of iterations
# gamma     A small value used for NNLSE solution (refer to the paper)
# lsm       Line search method; either step-halving (halving) or optimal (optim)

cgn.pi = function(mix, x, maxit=1, ll.mix, gamma=1e-3, eps=1e-3,
        lsm=c("halving","optim")) {
  convergence = 1
  if( missing(ll.mix) ) ll.mix = logLik(mix, x)
  for( i in 1:maxit ) {
    mix1 = mix
    ll.mix1 = ll.mix
    d = dnormmix(mix, x)
    d[d<1e-100] = 1e-100
    a = outer.dnorm( x, mix$mu, sqrt(mix$va) ) / d
    b = rep(2, length(x))
    r = nnlse(a, b, eps=eps)
    mix$pi = r$x / sum(r$x)
    r = line.normmix(mix1, mix, x, ll.mix1=ll.mix1, lsm=lsm)
    mix = r$mix
    ll.mix = r$ll.mix
  }
  list(mix=mix, ll.mix=ll.mix, num.iterations=i, convergence=convergence)
}


# ------------------------------------------------------------------ #
# Functions for handling mixtures of univariate normal distributions #
# ------------------------------------------------------------------ #

# The mixture components may have different variances

# Creates a "normmix" object

normmix = function(mu=NULL, pi=1/length(mu), va=1) {
  if( is.null(mu) ) mix = list(mu=NULL, pi=NULL, va=NULL)  # NULL mixture
  else {               # If mu == NA, one component with unknown mean
    k = max(length(mu), length(pi), length(va), na.rm=TRUE)
    mu = rep(mu, len=k)
    pi = rep(pi, len=k)
    va = rep(va, len=k)
    mix = list(mu=mu, pi=pi/sum(pi), va=va)
  }
  class(mix) = "normmix"
  mix
}

is.normmix = function(mix) class(mix) == "normmix"

is.null.normmix = function(mix) is.null(mix$mu)

print.normmix = function(mix, ...) {
  if( is.null(mix) ) b = matrix(nrow=0, ncol=3)
  else b = cbind(mix$mu, mix$pi, mix$va)
  colnames(b) = c("mu","pi","va")
  print(b, ...)
}

sort.normmix = function(mix) {
  if( is.null(mix) ) return(mix)
  index = order(mix$mu)
  mix$mu = mix$mu[index]
  mix$pi = mix$pi[index]
  mix$va = mix$va[index]
  mix
}

is.unsorted.normmix = function(mix) is.unsorted(mix$mu)

# unique normmix

unique.normmix = function(mix, prec=c(0.01,1e-6,0.01)) {
  if( length(mix$mu) == 1 ) return(mix)
  if( is.unsorted.normmix(mix) ) mix = sort.normmix(mix)
  if ( max(prec) == 0 ) return(mix)
  prec = rep(prec, len=3)
  mu2 = mu = mix$mu
  pi2 = pi = mix$pi
  va2 = va = mix$va
  j  = va <= 1e-10 * prec[3] | pi <= prec[2]
  mu = mu[!j]
  pi = pi[!j]
  va = va[!j]
  index = 0
  repeat {
    if( length(mu) == 0 ) break
    j = abs(mu[1] - mu) <=  (sqrt(va[1]) * prec[1]) &
        abs(va[1] - va) <=  (sqrt(va[1]) * prec[3])
    index = index + 1
    mu2[index] = weighted.mean(mu[j], pi[j])
    pi2[index] = sum( pi[j] )
    va2[index] = weighted.mean(va[j], pi[j])
    mu = mu[!j]
    va = va[!j]
    pi = pi[!j]
  }
  normmix(mu=mu2[1:index], pi=pi2[1:index], va=va2[1:index])
}

# Check whether two mixtures are approximately equal (or identical)

aeq.normmix = function(a, b, precision=c(1e-6,1e-6,1e-6)) {
   precision = rep(precision, len=3)
  if( is.null(a) && is.null(b) ) return( TRUE )
  if( is.null(a) || is.null(b) ) return( FALSE )
  if( length(a$mu) != length(b$mu) ) return( FALSE )
  if( max(abs(a$mu - b$mu) / sqrt(pmin(a$va, b$va))) >= precision[1] )
    return( FALSE )
  if( max(abs(a$pi - b$pi)) >= precision[2]) return( FALSE )
  if( max(abs((a$va - b$va) / sqrt(pmin(a$va, b$va)))) >= precision[3] )
    return( FALSE )
  TRUE
}

# Computes the values of normal densities

# mu and sd are for the same distributions and should be of the same length

outer.dnorm = function(x, mu=0, sd=1, log=FALSE) {
  if( log )
    sweep(dnorm(sweep(outer(x, mu, "-"), 2, sd, "/"), log=TRUE), 2, log(sd),
          "-")
  else
    sweep(dnorm(sweep(outer(x, mu, "-"), 2, sd, "/")), 2, sd, "/")
}

# Density function of a normal mixture

dnormmix = function( mix, x, log=FALSE ) {
  if( log ) {
    logd = outer.dnorm( x, mix$mu, sqrt(mix$va), log=TRUE )
    ma = apply(logd, 1, max)
    ma + log( rowSums(sweep(exp(sweep(logd, 1, ma, "-")), 2, mix$pi, "*")) )
  }
  else rowSums( sweep(outer.dnorm(x, mix$mu, sqrt(mix$va)), 2, mix$pi, "*") )
}

# Returns a random sample of a mixture #

# n        Sample size
# mu       Means
# pi       Proportions for all components
# va       Variance
# fixed    = TRUE, proportions of observations drawn from each component
#          is fixed (as much as possible)

rnormmix = function(n=50, mu=c(0,4), pi=1, va=1, mix, fixed=TRUE) {
  if( ! missing(mix) )
    return( rnormmix( n=n, mu=mix$mu, pi=mix$pi, va=mix$va, fixed=fixed ) )

  if( n == 0 ) return( numeric(0) )
  k = length(mu)
  va = rep(va, len=k)
  pi = rep(pi, len=k)
  pi = pi / sum(pi)
  if( fixed ) {
    nj = floor(n * pi)
    x = rnormmix(n - sum(nj), mu=mu, pi=pi, fixed=FALSE)
  }
  else {
    cs = colSums(outer(runif(n), cumsum(pi), "<"))
    nj = c(cs[1], diff(cs))
    x = numeric(0)
  }
  for( j in 1:k )
    x = c(x, rnorm(nj[j], mu[j], sqrt(va[j])))
  if( n > 1 ) sample(x)
  else x
}

# Log-likelihood

logLik.normmix = function( mix, x, w=1 ) sum( w * dnormmix(mix, x, log=TRUE) )

# Line search

# mix1      Current normmix object
# mix2      Next normmix object; same number of components as in mix1
# x         Observations
# ll.mix1   log-likelihood of mix1
# lsm       Line search method, either step-halving (havling) or optimal (optim)

line.normmix = function(mix1, mix2, x, ll.mix1=NULL, tol=1e-10,
            lsm=c("halving","optim") ) {
  ll.alpha = function(alpha) {
    m = normmix( (1-alpha) * mix1$mu + alpha * mix2$mu,
      (1-alpha) * mix1$pi + alpha * mix2$pi,
      ((1-alpha) * sqrt(mix1$va) + alpha * sqrt(mix2$va))^2 )
    logLik(m, x)
  }
  if( is.null(ll.mix1) ) ll.mix1 = logLik(mix1, x)
  convergence = 0
  alpha = 1
  repeat {
    if( lsm == "optim") {
      m = normmix( (1-alpha) * mix1$mu + alpha * mix2$mu,
        (1-alpha) * mix1$pi + alpha * mix2$pi,
        ((1-alpha) * sqrt(mix1$va) + alpha * sqrt(mix2$va))^2 )
      d = c(m$pi - mix1$pi, m$mu-mix1$mu)
      gll = sum(gradient.normmix(m, x)$gradient * d)
      if( is.na(gll) ) alpha = alpha / 2
      else if(gll < 0)
        alpha = optimize(ll.alpha, lower=0, upper=alpha, max=TRUE,
          tol=1e-3)$maximum
    }
    ll.mix = ll.alpha(alpha)
    if( ll.mix > ll.mix1 ) break
    if(alpha < tol) {convergence = 1; print("Error: Line search failed"); break}
    alpha = alpha / 2
  }
  list(mix = normmix((1-alpha) * mix1$mu + alpha * mix2$mu,
         (1-alpha) * mix1$pi + alpha * mix2$pi,
         ((1-alpha) * sqrt(mix1$va) + alpha * sqrt(mix2$va))^2 ),
       ll.mix=ll.mix, alpha=alpha, convergence=convergence)
}

# Gradient of the log-likelihood (not the gradient function)

# mix          Normmix object
# x            Observed points
# individual   Returns the gradients of individual log-likelihoods

gradient.normmix = function(mix, x, individual=FALSE) {
  logd = outer.dnorm( x, mix$mu, sqrt(mix$va), log=TRUE )
  pi.j = exp(sweep(logd, 1, apply(logd, 1, max), "-"))    # avoids overflow
#   pi.j = outer.dnorm( x, mix$mu, sqrt(mix$va) )         # can cause overflow
  p.pi = pi.j / as.vector(pi.j %*% mix$pi)
  p = sweep(p.pi , 2, mix$pi, "*")
  p[p<1e-100] = 1e-100
  x.mu = outer(x, mix$mu, "-")
  x.mu.va = sweep(x.mu, 2, mix$va, "/")
  p.x.mu.va = p * x.mu.va
  gradient.i = cbind(p.pi, p.x.mu.va)
  if( individual ) list( gradient = gradient.i )
  else list( gradient = colSums(gradient.i) )
}

# The gradient function

normmix.gradient = function(mu, mix, x, sd=NULL) {
  if( is.null(sd) ) sd = sqrt(mix$va[1])
  if( ! is.vector(mix) ) mix = dnormmix(mix, x)
  mix[mix<1e-100] = 1e-100
  colSums(outer.dnorm(x, mu, sd) / mix) - length(x)
}

# First Direcvative of the gradient function

normmix.dgradient = function(mu, mix, x, sd=NULL) {
  if( is.null(sd) ) sd = sqrt(mix$va[1])
  if( ! is.vector(mix) ) mix = dnormmix(mix, x)
  mix[mix<1e-100] = 1e-100
  colSums(outer.dnorm(x, mu, sd) * outer(x, mu, "-") / sd^2 / mix)
}

# Second direcvative of the gradient function

normmix.d2gradient = function(mu, mix, x, sd=NULL) {
  if( is.null(sd) ) sd = sqrt(mix$va[1])
  if( ! is.vector(mix) ) mix = dnormmix(mix, x)
  mix[mix<1e-100] = 1e-100
  colSums(outer.dnorm(x, mu, sd) * (outer(x, mu, "-")^2-sd^2) / sd^4 / mix)
}

# Find all local maxima of the gradient function g(theta, mix, x)

# mix     A normmix object or its likelihood vector at x
# x       Sample
# sd      Standard deviation; if =null, SD of the first component will be used.
# check.points    Points at which the gradient values are checked
# tol     Tolerance level
# maxit   Maximum number of iterations

maxima.gradient = function(mix, x, sd=NULL, check.points=100, tol=0,
          maxit=100) {
  if( length(check.points) == 1 )
    check.points = seq(min(x)-.1, max(x)+.1, length=check.points)
  np = length(check.points)
  dg = normmix.dgradient(check.points, mix, x, sd=sd)
  jmax = (1:(np-1))[dg[1:(np-1)] > 0 & dg[2:np] < 0]
  if( length(jmax) < 1 ) return
  mu = (check.points[jmax] + check.points[jmax+1])/2
  left = check.points[jmax]
  right = check.points[jmax+1]
  for( i in 1:maxit ) {
    mu.old = mu
    d1 = normmix.dgradient(mu, mix, x, sd=sd)
    d2 = normmix.d2gradient(mu, mix, x, sd=sd)
    left[d1>0] = mu[d1>0]
    right[d1<0] = mu[d1<0]
    mu = mu - d1 / d2
    j = is.na(mu) | mu < left | mu > right
    mu[j] = (left[j] + right[j]) / 2
    if( max(abs(mu - mu.old)) < 1e-6 ) break
  }
  g = normmix.gradient(mu, mix, x, sd=sd)
  j = g >= tol
  list(mu=mu[j], gradient=g[j])
}

# -------------- #
# NNLS and NNLSE #
# -------------- #

# Least squares solution under nonnegativity constraint:
#
#        Minimize     ||Ax - b||
#        subject to   x >= 0

# This is an interface function to the Fortran subroutine NNLS that is
# downloaded from http://www.netlib.org/lawson-hanson

nnls = function(a,b) {
  m = as.integer(dim(a)[1])
  n = as.integer(dim(a)[2])
  storage.mode(a) = "double"
  storage.mode(b) = "double"
  x = double(n)
  rnorm = double(1)
  w = x
  zz = b
  index = integer(n)
  mode = integer(1)
  .Fortran("nnls",a,m,m,n,b,x=x,rnorm=rnorm,w,zz,index=index,
           mode=mode)[c("x","rnorm","index","mode")]
}

# Least squares solution under both nonnegativity and equality constraints
#
#           Minimize    || A x - b ||
#           subject to  C x = d
#           and         x >= 0

nnlse = function(a, b, c=rep(1,ncol(a)), d=1, eps=1e-3)
  nnls(rbind(c, a * eps), c(d, b * eps))




# ============================== #
# Functions for Poisson mixtures #
# ============================== #


# ----------------------------------------------------------------------- #
# The constrained Newton method (CNM) for computing the NPMLE of a mixing #
# distribution with Poisson components                                    #
# ----------------------------------------------------------------------- #

# x         Observations
# verb      Verbosity level (0-4)
# prec      Precision level to combine components
# maxit     Maximum number of iterations
# check.points   Number of grid points for computing the derivatives of
#           the gradient function
# tol       Tolerance level for sup gradient
# mix       Initial estimate of the NPMLE
# method    Either CNM or CN1
# lsm       Line search method, either halving (step-halving) or optim (optimal)

cnm.poismix = function(x, w=1, verb=0, prec=1e-6,
    maxit=1000, check.points=100, tol=1e-6,
    mix=poismix(NA,pi=1), lsm=c("halving","optim"), ... ) {

  lsm = match.arg(lsm)
  if( length(check.points) == 1 )
    check.points = seq(max(min(x),1e-100), max(x), length=check.points)
  if( any(is.na(mix$lambda)) )
    mix = poismix( lambda=seq(min(x), max(x), by=5) )
  prec = c(1e-2, 0)
  mix$lambda[mix$lambda==0] = 1e-100
  w = rep(w, len=length(x))
  u = unique(sort(x))
  w = as.vector(outer(u, x, "==") %*% w)
  x = u
  n = sum(w)
  ll.mix = -Inf
  convergence = 1
  gradient = Inf
  for( i in 1:maxit ) {
    mix1 = mix
    ll.mix1 = ll.mix
    g = maxima.gradient.poismix(x=x, w=w, mix=mix, check=check.points,
      tol=-Inf)
    if( max(g$gradient) < tol ) { convergence = 0; break }
    nsp = g$lambda
    mix = poismix(c(mix$lambda,nsp), c(mix$pi,rep(0,length(nsp))))
    r = cgn.poismix2(mix, x, w=w, lsm=lsm,
      ll.mix=ll.mix)
    mix = r$mix
    ll.mix = r$ll
    mix = unique(mix, prec=c(.01,0))
    switch(verb,
           {cat("lambda = \n"); print(mix$lambda)},
           print(mix),
           cat("log-likelihood =", logLik(mix, x, w), "\n"),
           cat("Maximum gradient = ", max(g$gradient), "\n"))
  }
  mix = unique(sort.poismix(mix), prec=prec)
  g = maxima.gradient.poismix(mix=mix, x=x, w=w, check.points=check.points,
    tol=-Inf)
  mix$lambda[mix$lambda < 1e-10] = 0
  mix = unique.poismix(sort.poismix(mix), prec=prec)
  list(mix=mix, ll = logLik(mix, x, w), num.iterations=i,
       max.gradient=max(g$gradient), convergence=convergence)
}

# Updates the mixing proportions

cgn.poismix2 = function(mix, x, w, maxit=1,
                        lsm=c("halving","optim"), ll.mix) {
  lsm = match.arg(lsm)
  convergence = 1
  k = length(mix$lambda)
  if( missing(ll.mix) ) ll.mix = logLik(mix, x, w)
  for( i in 1:maxit ) {
    mix1 = mix
    ll.mix1 = ll.mix
    dij = outer(x, mix$lambda, dpois)
    di = as.vector(dij %*% mix$pi)
    dl.pi = sqrt(w) * dij / di
    a = dl.pi
    b = 2 * sqrt(w)
    r = nnlse(a, b, eps=1e-3)
    mix$pi = r$x / sum(r$x)
    r = line.poismix(mix1, mix, x, w=w, ll.mix1=ll.mix, lsm=lsm)
    mix = r$mix
    ll.mix = r$ll.mix
  }
  list(mix=mix, ll=ll.mix, num.iterations=i)
}

# The gradient function

poismix.gradient = function(lambda, mix, x, w=1) {
  w = rep(w, len=length(x))
  colSums(w * outer(x, lambda, dpois) / dpoismix(mix, x)) -
    sum( w )
}

# Direcvative of the gradient function

poismix.dgradient = function(lambda, mix, x, w=1) {
  w = rep(w, len=length(x))
  colSums(w * outer(x, lambda, dpois) * (outer(x, lambda, "/")-1) /
          dpoismix(mix, x)  )
}

# Second direcvative of the gradient

poismix.d2gradient = function(lambda, mix, x, w=1) {
  colSums(w * outer(x, lambda, dpois) *
          ((outer(x, lambda, "/")-1)^2 - outer(x, lambda^2, "/")) /
          dpoismix(mix, x)  )
}

# Find all local maxima of the gradient function g(theta, mix, x)

# L             Likelihood vector
# x             Data
# mix           A "poismix" object
# check.points: If an integer, number of check points at which the first
#               derivative of the gradient will be first examined.
#               If a vector, the actual check points.
# tol:          The lower bound of the gradient value for a local maximum

maxima.gradient.poismix = function(L, x, w, mix, check.points=100, tol=0) {
  if( length(check.points) == 1 )
    check.points = seq(max(min(x),1e-100), max(x), length=check.points)
  np = length(check.points)
  if( missing(L) ) L = dpoismix(mix, x)
  dg = colSums(w * outer(x, check.points, dpois) *
    (outer(x, check.points, "/")-1) / L)
  jmax = (1:(np-1))[dg[1:(np-1)] > 0 & dg[2:np] < 0]
  g1 = colSums(w * outer(x, check.points[1], dpois) / L) - sum( w )
  if( length(jmax) < 1 ) {
    if(dg[1] < 0 && g1 > tol )
      return( list(lambda=check.points[1], gradient=g1) )
    else return
  }
  lambda = (check.points[jmax] + check.points[jmax+1])/2
  repeat {
    lambda.old = lambda
    lambda = lambda -
      colSums(w * outer(x, lambda, dpois) * (outer(x, lambda, "/")-1) / L) /
        colSums(w * outer(x, lambda, dpois) *
                ((outer(x, lambda, "/")-1)^2 - outer(x, lambda^2, "/")) / L)
    lambda[lambda<check.points[1]] = check.points[1]
    if( max(abs(lambda - lambda.old)) < 1e-6 ) break
  }
  g = colSums(w * outer(x, lambda, dpois) / L) - sum(w)
  j = g >= tol
  g1 = colSums(w * outer(x, check.points[1], dpois) / L) - sum( w )
  if(dg[1] < 0 && g1 > 0 )
    list(lambda=c(check.points[1],lambda[j]), gradient=c(g1,g[j]))
  else list(lambda=lambda[j], gradient=g[j])
}


# --------------------------------- #
# Mixtures of Poisson Distributions #
# --------------------------------- #

# Creates a "poismix" object

poismix = function(lambda=NULL, pi=1/length(lambda)) {
  if( is.null(lambda) ) a = list(lambda=NULL, pi=NULL)  # NULL mixture
  else {               # If lambda == NA, one component with unknown mean
    k = max(length(lambda), length(pi), na.rm=TRUE)
    lambda = rep(lambda, len=k)
    pi = rep(pi, len=k)
    a = list(lambda=lambda, pi=pi/sum(pi))
  }
  class(a) = "poismix"
  a
}

is.null.poismix = function(a) is.null(a$lambda)

print.poismix = function(a) {
  if( is.null(a$lambda) ) b = matrix(nrow=0, ncol=2)
  else b = cbind(a$lambda, a$pi)
  colnames(b) = c("lambda","pi")
  print(b)
}

sort.poismix = function(a) {
  if( is.null(a) ) return(a)
  i = order(a$lambda)
  a$lambda = a$lambda[i]
  a$pi = a$pi[i]
  a
}

is.unsorted.poismix = function(a) is.unsorted(a$lambda)

# Unique

unique.poismix = function(a, precision=c(0.01,1e-6)) {
  if( length(a$lambda) == 1 ) return(a)
  if( is.unsorted.poismix(a) ) a = sort.poismix(a)
  if ( max(precision) == 0 ) return(a)
  precision = rep(precision, len=3)
  count = i1 = i2 = 1
  for( i in 2:length(a$lambda) ) {
    if( a$lambda[i] - a$lambda[count] > precision[1] &&
       min(max(a$pi[count:(i-1)]), a$pi[i]) > precision[2] ) {
      a$lambda[count] = weighted.mean(a$lambda[i1:i2], a$pi[i1:i2])
      a$pi[count] = sum(a$pi[i1:i2])
      count = count + 1
      i1 = i2 = i
      a$lambda[count] = a$lambda[i]
    }
    else i2 = i
  }
  a$lambda[count] = weighted.mean(a$lambda[i1:i2], a$pi[i1:i2])
  a$pi[count] = sum(a$pi[i1:i2])
  poismix(lambda=a$lambda[1:count], pi=a$pi[1:count])
}

# Approximately equal

aeq.poismix = function(a, b, precision=c(1e-6,1e-6)) {
  if( is.null(a) && is.null(b) ) return( TRUE )
  if( is.null(a) || is.null(b) ) return( FALSE )
  if( length(a$lambda) != length(b$lambda) ) return( FALSE )
  if( max(abs(a$lambda - b$lambda) / sqrt(pmin(a$lambda, b$lambda)+1e-20))
     >= precision[1] )
    return( FALSE )
  if( max(abs(a$pi - b$pi)) >= precision[2]) return( FALSE )
  TRUE
}

# Returns a random sample of a Poisson mixture

# n        Sample size
# lambda   Means
# pi       Proportions for all components
# fixed    = TRUE, proportion of observations from each component is fixed

rpoismix = function(n=50, lambda=c(1,4), pi=1, mix, fixed=TRUE) {
  if( ! missing(mix) ) {lambda =  mix$lambda; pi=mix$pi}
  if( n == 0 ) return( numeric(0) )
  k = length(lambda)
  pi = rep(pi, len=k)
  pi = pi / sum(pi)
  if( fixed ) {
    nj = floor(n * pi)
    x = rpoismix(n - sum(nj), lambda=lambda, pi=pi, fixed=FALSE)
  }
  else {
    cs = colSums(outer(runif(n), cumsum(pi), "<"))
    nj = c(cs[1], diff(cs))
    x = numeric(0)
  }
  for( j in 1:k ) x = c(x, rpois(nj[j], lambda[j]))
  sample(x)
}

# ---------------------- #
# Distribution functions #
# ---------------------- #

# Density function of a Poisson mixture

dpoismix = function(a, x, log=FALSE) {
  log.dpois = function(x, lambda) dpois(x, lambda, log=TRUE)
  if (log) {
    logd = outer(x, a$lambda, log.dpois)
    ma = apply(logd, 1, max)
    ma + log(rowSums(sweep(exp(sweep(logd, 1, ma, "-")), 2, a$pi, "*")))
  }
  else rowSums(sweep(outer(x, a$lambda, dpois), 2, a$pi, "*"))
}

# log-likelihood

logLik.poismix = function( a, x, w=1 ) sum( w * dpoismix(a, x, log=TRUE) )

# Line search

# mix1      Current poismix object
# mix2      New poismix object
# x         Observations
# w         Frequency
# ll.mix1   log-likelihood of mix1
# lsm       Line search method, either step-halving (halving) or optimal (optim)

line.poismix = function(mix1, mix2, x, w, ll.mix1=NULL,
       lsm=c("halving","optim")) {
  lsm = match.arg(lsm)
  ll.alpha = function(alpha) {
    m = poismix( (1-alpha) * mix1$lambda + alpha * mix2$lambda,
      (1-alpha) * mix1$pi + alpha * mix2$pi )
    logLik(m, x, w)
  }
  if( is.null(ll.mix1) ) ll.mix1 = logLik(mix2, x, w)
  convergence = 0
  alpha = 1
  repeat {
    if( lsm == "optim" ) {
      m = poismix( (1-alpha) * mix1$lambda + alpha * mix2$lambda,
        (1-alpha) * mix1$pi + alpha * mix2$pi )
      d = m$pi - mix1$pi
      gll = sum(gradient.poismix(m, x, w)$gradient * d)
      if( gll < 0 ) alpha = optimize(ll.alpha, lower=0, upper=alpha, max=TRUE,
                      tol=1e-3)$maximum
    }
    ll.mix = ll.alpha(alpha)
    if( ll.mix > ll.mix1 ) break
    if( alpha < 1e-10 ) {convergence = 1; break}
    alpha = alpha / 2
  }
  list(mix = poismix((1-alpha) * mix1$lambda + alpha * mix2$lambda,
         (1-alpha) * mix1$pi + alpha * mix2$pi),
       ll.mix=ll.mix, alpha=alpha, convergence=convergence)
}

# Gradient of the log-likelihood

gradient.poismix = function(mix, x, w, individual=FALSE,
         hessian=FALSE) {
  dij = outer(x, mix$lambda, dpois)
  di = as.vector(dij %*% mix$pi)
  gradient.i = w * dij / di
  if( individual ) list( gradient = gradient.i )
  else list( gradient = colSums(gradient.i) )
}
