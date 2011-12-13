# =========================== #
# The common variance problem #
# =========================== #

# ========== #
# Class cvp2 #
# ========== #

# The common variance problem

rcvp = function(k=100, ni=3:5, mu=c(0,4), pr=c(.7,.3), sd=3, i1=1) {
  if(length(mu) == 1) mu = rep(mu, k)
  else mu = sample(mu, k, replace=TRUE, prob=pr)
  ni = rep(ni, len=k)
  nicum = c(0,cumsum(ni))
  n = nicum[k+1]
  data = matrix(nrow=n, ncol=2)
  dimnames(data) = list(NULL, c("group", "x"))
  for(i in 1:k) {
    j = (nicum[i]+1):(nicum[i+1])
    data[j,1] = i1 + i - 1          # sample(k, 1)
    data[j,2] = rnorm(ni[i], mu[i], sd)
  }
  data
}

# The common variance problem with data stored in the sufficient statistics
# (ni, mi, ri) for each group i, denoting the number of observations, the
# mean and the residual sum of squares, respectively.

cvp2 = function(x) {
  ni = tapply(x[,2], x[,1], length)
  group = as.numeric(names(ni))
  ni = as.vector(ni)
  mi = as.vector(tapply(x[,2], x[,1], mean))
  ri = as.vector(tapply(x[,2], x[,1], var) * (ni-1))
  ri[is.na(ri)] = 0
  names(ni) = names(mi) = names(ri) = NULL
  x2 = list(group=group, ni=ni, mi=mi, ri=ri)
  class(x2) = "cvp2"
  x2
}

# cvp3 = function(x) {
#   group = x[,1]
#   d = x[,2]
#   d = data.frame(d)
#   ni = aggregate(d, by=list(group=group), FUN=length)
#   mi = aggregate(d, by=list(group=group), FUN=mean)[,2]
#   ri = aggregate(d, by=list(group=group), FUN=var)[,2] * (ni[,2]-1)
#   ri[is.na(ri)] = 0
#   x2 = list(group=as.numeric(as.character(ni[,1])), ni=ni[,2],
#     mi=mi, ri=ri)
#   class(x2) = "cvp2"
#   x2
# }

rcvp2 = function(k=100, ni=3:5, mu=c(0,4), pr=c(.7,.3), sd=3)
  cvp2( rcvp(k=k, ni=ni, mu=mu, pr=pr, sd=sd) )

# Random duplications

# x      Data
# d      Mixing distribution
# sd     Standard deviation
# nd     Number of duplications
#

# rdup.cvp2 = function(x, d, sd, nd)
#   rcvp2(nd*length(x), x$ni, d$pt, d$pr, sd)

print.cvp2 = function(x, ...)
  print(cbind(group=x$group, ni=x$ni, mi=x$mi, ri=x$ri), ...)

length.cvp2 = function(x) length(x$ni)

range.cvp2 = function(x, beta, ...) range(x$mi, ...)

initial.cvp2 = function(x, beta=NULL, mix=NULL, kmax=NULL) {
  if(is.null(beta)) beta = sqrt(sum(x$ri) / (sum(x$ni) - length(x)))
  if(length(beta) != 1 || beta <= 0) stop("initial beta incorrect")
  if(is.null(kmax)) kmax = 10
  if(is.null(mix) || is.null(mix$pt))
    mix = dden(unique(quantile(x$mi, p=seq(0,1,len=kmax), type=1)))
  list(beta=beta, mix=mix)
}

valid.cvp2 = function(x, beta) beta > 0

# Compute the log density and its derivatives, for each combination of x[i]
# and pt[j].

logd.cvp2 = function(x, beta, pt, which=c(1,0,0,0)) {
  dl = vector("list", 4)
  names(dl) = c("ld","db1","dt1","dt2")
  if(which[1] == 1)
    dl$ld = - x$ni / 2 * log(2*pi*beta^2) -
      (x$ri + x$ni * outer(x$mi, pt, "-")^2) / beta^2 / 2
  if(which[2] == 1)
    dl$db1 = array(- x$ni / beta +
      (x$ri + x$ni * outer(x$mi, pt, "-")^2) / beta^3,
      dim=c(length(x$mi), length(pt), 1))
  if(which[3] == 1)
    dl$dt1 = x$ni / beta^2 * outer(x$mi, pt, "-")
  if(which[4] == 1)
    dl$dt2 = matrix(-x$ni / beta^2, nrow=length(x$ni), ncol=length(pt))
  dl
}



# ======================================






# parameter values must be valid in the functions below

lf1.cvp2 = function(x, beta, pt) {
  if( beta <= 0 ) rep(-Inf, nrow=length(x$mi))
  else - x$ni / 2 * log(2*pi*beta^2) -
    (x$ri + x$ni * (x$mi - pt)^2) / beta^2 / 2
}

dlb1.cvp2 = function(x, beta, pt, order=1) {
  dl = vector("list", length(order))
  names(dl) = paste("d", order, sep="")
  if(any(order == 0))
    dl$d0 = - x$ni / 2 * log(2*pi*beta^2) -
      (x$ri + x$ni * (x$mi - pt)^2) / beta^2 / 2
  if(any(order == 1))
    dl$d1 = - x$ni / beta + (x$ri + x$ni * (x$mi - pt)^2) / beta^3
#   if(any(order == 2)) {
#     dl$d2 = array(x$ni / beta^2 -
#       3 * (x$ri + x$ni * outer(x$mi, pt, "-")^2) /  beta^4,
#       dim=c(length(x$mi), length(pt), 1, 1))
#   }
  dl
}

dlt1.cvp2 = function(x, beta, pt, order=1) {
  dl = vector("list", length(order))
  names(dl) = paste("d", order, sep="")
  if(any(order == 0))
    dl$d0 = - x$ni / 2 * log(2*pi*beta^2) -
      (x$ri + x$ni * (x$mi - pt)^2) / beta^2 / 2
  if(any(order == 1)) dl$d1 = x$ni / beta^2 * (x$mi - pt)
  if(any(order == 2)) dl$d2 = -x$ni / beta^2
  dl
}

lf.cvp2 = function(x, beta, pt) {
  if( beta <= 0 ) matrix(-Inf, nrow=length(x$mi), ncol=length(pt))
  else - x$ni / 2 * log(2*pi*beta^2) -
    (x$ri + x$ni * outer(x$mi, pt, "-")^2) / beta^2 / 2
}

dlb.cvp2 = function(x, beta, pt, order=1) {
  dl = vector("list", length(order))
  names(dl) = paste("d", order, sep="")
  if(any(order == 0))
    dl$d0 = - x$ni / 2 * log(2*pi*beta^2) -
      (x$ri + x$ni * outer(x$mi, pt, "-")^2) / beta^2 / 2
  if(any(order == 1)) {
    dl$d1 = array(- x$ni / beta +
      (x$ri + x$ni * outer(x$mi, pt, "-")^2) / beta^3,
      dim=c(length(x$mi), length(pt), 1))
  }
#   if(any(order == 2)) {
#     dl$d2 = array(x$ni / beta^2 -
#       3 * (x$ri + x$ni * outer(x$mi, pt, "-")^2) /  beta^4,
#       dim=c(length(x$mi), length(pt), 1, 1))
#   }
  dl
}

dlt.cvp2 = function(x, beta, pt, order=1) {
  dl = vector("list", length(order))
  names(dl) = paste("d", order, sep="")
  if(any(order == 0))
    dl$d0 = - x$ni / 2 * log(2*pi*beta^2) -
      (x$ri + x$ni * outer(x$mi, pt, "-")^2) / beta^2 / 2
  if(any(order == 1))
    dl$d1 = x$ni / beta^2 * outer(x$mi, pt, "-")
  if(any(order == 2))
    dl$d2 = matrix(-x$ni / beta^2, nrow=length(x$ni), ncol=length(pt))
  dl
}

plot.cvp2 = function(x, beta, mix, hist.col=NULL, len=200, col="red",
     comp.col=col, breaks=30, ...) {
  data = rep(x$mi,x$ni)
  hist(data, freq=FALSE, main="", breaks=min(breaks, round(length(data)/4)),
       col=hist.col, ...)
  xlim = range(x) + c(-1,1) * beta
  y = seq(xlim[1], xlim[2], len=len)
  d = 0
  n = sum(x$ni)
  for(i in 1:length(x)) {
    d = d + dnorm(y, x$mi[i], sd=beta/sqrt(x$ni[i])) * x$ni[i] / n
  }
  lines(y, d, col="blue")
#   f = mix$pr * exp(lf(x, beta, mix$pt))

  m = length(mix$pt)
  d = 0
  for( i in 1:m ) {
    di = mix$pr[i] * dnorm(y, mix$pt[i], beta/sqrt(4))
    lines(y, di, col=comp.col, lty=2)
    d = d + di
  }
  lines(y, d, col=col, ...)
  points(mix$pt, rep(0,m), col=comp.col)
}

# d2lt.cvp2 = function(x, beta, pt)
#   matrix(-x$ni / beta^2, nrow=length(x$ni), ncol=length(pt))





