######################################################################
## Reproduceable example
######################################################################

#library(compiler)
library(Rcpp)
#library(rbenchmark)
library(abind)
library(fGarch)
source("~/Dropbox/Uni Work/Masters Thesis/Codes/Constraint Newton Method/cnm.R")
source("~/Dropbox/Uni Work/Masters Thesis/Codes/Constraint Newton Method/dden.R")

#source("~/work/spmle/R/old/cnm.R")
#source("~/work/spmle/R/old/dden.R")

## New function after improvement
logd.mgarch <- function(xt, beta, pt, which){
    xt <- as.numeric(xt)
    T <- length(xt)
    lpt <- length(pt)
    lb <- length(beta)
    dl <- vector("list", length = 4)
    names(dl) <- c("ld", "db1", "dt1", "dt2")
    betaSum <- as.numeric(filter(xt[-T]^2, beta[2], "recursive"))
    sigma.t <-
        sqrt(c(beta[4]^2,
               beta[1] * ((1 - cumprod(rep.int(beta[2], T - 1))))/
               (1 - beta[2]) +
               cumprod(rep.int(beta[2], T - 1)) * beta[4]^2 +
               beta[3] * betaSum))
    if(which[1] == 1){
        dl$ld <- matrix(-(0.5 * log(2 * pi) + log(sigma.t) +
                          rep(log(pt), each = T) +
                          ((xt - beta[5])^2/sigma.t^2)/rep(2 * pt^2, each = T)),
                        nr = T, nc = lpt)
    }
    if(which[2] == 1){
        dldsigma <- matrix((xt - beta[5])^2/sigma.t^3/
                           rep(pt^2, each = T) - 1/sigma.t,
                           nr = T, nc = lpt)
        sig.vec <- 2 * sigma.t
        convFilter <- 1:(T - 2) * cumprod(c(1, rep(beta[2], T - 3)))
        myconvolve <- beta[3] * convolve(convFilter, rev(xt^2),
                                      type = "open")[1:(T - 2)]
        cp.beta <- c(0, 0, myconvolve)
        cp.beta[1:2] <- 0
        dsigmadalpha0 <- cumsum(c(0, beta[2]^(0:(T-2)))) #1
        dsigmadalpha1 <-
            (beta[1] *
             cumsum(c(0, 0, 1:(T - 2) * beta[2]^(0:(T - 3))))) +
                 c(0, 1:(T - 1) * beta[2]^(0:(T - 2)) * beta[4]^2) + cp.beta #2
        dsigmadbeta1 <- c(0, betaSum) #3
        dsigmadsigma <- 2 * beta[4] * beta[2]^(0:(T - 1)) #4
        dldmu <- (xt - beta[5])/sigma.t^2 / rep(pt^2, each = T) #5
        dim(dldmu) <- c(T, lpt)
        dbvec <- array(c(rep.int(dsigmadalpha0, lpt),
                         rep.int(dsigmadalpha1, lpt),
                         rep.int(dsigmadbeta1, lpt),
                         rep.int(dsigmadsigma, lpt)),
                       dim = c(T, lpt, lb - 1))
        dldsigma <- array(dldsigma/sig.vec, dim = c(T, lpt, lb - 1))
        dl$db1 <- abind(dldsigma * dbvec, dldmu, along = 3)
    }
    if(which[3] == 1){
      dt1 = (xt - beta[5])^2/sigma.t^2 / rep(pt^3, each=T) - 1/rep(pt, each=T)
      dim(dt1) = c(T,lpt)
      dl$dt1 = dt1
    }
    if(which[4] == 1){
        dt2 <- (-3 * (xt - beta[5])^2)/sigma.t^2 /rep(pt^4, each = T) +
            1/(rep(pt^2, each = T))
        dim(dt2) = c(T, lpt)
        dl$dt2 = dt2
    }
    dl
}

logd.mgarch <- cmpfun(logd.mgarch)

valid.mgarch <- function(x, beta, mix){
    beta[1] > 0 &&
    beta[2] >= 0 &&
    beta[3] >= 0 &&
    beta[4] >= 0 &&
    mix$pt > 0
}

valid.snpmle <- function(x, beta, mix)
  valid(x, beta, mix) && all(mix$pr >= 0) && all(mix$pt >= 0)


initial.mgarch <- function(x, beta = NULL, mix = NULL, kmax = NULL){
    if(is.null(beta)){
        cgf <- coef(garchFit(data = as.numeric(x), trace = FALSE))
        beta <- c(cgf[2], cgf[4], cgf[3], sd(x), mean(x))
        names(beta) <- c("omega", "beta1", "alpha1", "sigma0", "mu")
        mix <- dden(1, 1)
        list(beta = beta, mix = mix)
    }
}


######################################################################
## modify some of the original function
######################################################################

maxgrad <- function(x, beta, dmix, ma, grid=100, tol=-Inf, maxit=100) {
  if(length(grid) == 1){
        ##rth <- range.mgarch(x, beta)
    grid <- seq(0.1, 10, length = grid)
  }
  np = length(grid)
  # print(grid)
  dg = grad(x, grid, beta, dmix, ma, order=1)$d1

  # d0 = grad(x, grid, beta, dmix, ma, order=0)$d0
  # print(rbind(grid,d0, dg))
  # stop()

  jmax = (1:(np-1))[dg[1:(np-1)] > 0 & dg[2:np] < 0]
  if( length(jmax) < 1 ) return
  pt = (grid[jmax] + grid[jmax+1]) * .5
  left = grid[jmax]
  right = grid[jmax+1]
  if(length(pt) != 0) {
    pt.old = left
    d1.old = grad(x, left, beta, dmix, ma, order=1)$d1
    d2 = rep(-1, length(pt))  # or d2 = rep(1, length(pt))
    for( i in 1:maxit ) {
      d1 = grad(x, pt, beta, dmix, ma, order=1)$d1
      d2t = (d1 - d1.old) / (pt - pt.old)
      jd = !is.na(d2t) & d2t < 0
      d2[jd] = d2t[jd]
      left[d1>0] = pt[d1>0]
      right[d1<0] = pt[d1<0]
      pt.old = pt
      d1.old = d1
      pt = pt - d1 / d2
      j = is.na(pt) | pt < left | pt > right
      pt[j] = (left[j] + right[j]) * .5
      # print(pt)
      if( max(abs(pt - pt.old)) <= 1e-14 * diff(range(grid))) break
    }
  }
  else i = 0
  # print(i)
  if(dg[np] >= 0) pt = c(grid[np], pt)
  if(dg[1] <= 0) pt = c(grid[1], pt)
  if(length(pt) == 0) stop("no new support point found") # should not happen
  g = grad(x, pt, beta, dmix, ma, order=0)$d0
  names(pt) = names(g) = NULL
  j = g >= tol
  list(pt=pt[j], grad=g[j], num.iterations=i)
}


plotgrad <- function(x, beta, mix, ma, len=500,
      xlab=expression(theta), ylab,
      cex=1, pch=1, order=0, lower, upper, ...) {
  if(missing(ylab)) {
    ylab = switch(order+1,
    expression(d(theta * "; " * G, beta)),
    expression(d[1](theta * "; " * G, beta)),
    expression(d[2](theta * "; " * G, beta))  )
  }
  if( missing(lower) || missing(upper) ) {
    rth = c(0, 20)
    if( missing(lower) ) lower = rth[1] # - .05 * diff(rth)
    if( missing(upper) ) upper = rth[2] # + .05 * diff(rth)
  }
  pt = seq(lower, upper, len=len)
  g = switch(order+1,
    grad(x, pt, beta, mix, ma, order=order)$d0,
    grad(x, pt, beta, mix, ma, order=order)$d1,
    grad(x, pt, beta, mix, ma, order=order)$d2)
  plot(pt, g, type="l", col="blue", xlab=xlab, ylab=ylab,
       cex = cex, cex.axis = cex, cex.lab = cex, ... )
  if(is.dden(mix)) {
    j = mix$pr != 0
    points(mix$pt[j], rep(0,length(mix$pt[j])), pch=pch, col="red")
    abline(v=mix$pt[j], lty=3, col="red")
  }
  lines(c(lower, upper), c(0,0), col="black")
}


cnmms <- function(x=rcvp2(), init=NULL, maxit=1000,
                  model=c("spmle","npmle"),
                  tol=1e-10, grid=100, kmax=Inf,
                  plot=c("null", "gradient", "prob","dden"),
                  plotorder=0, verb=0, llt=NULL) {
    plot = match.arg(plot)
    model = match.arg(model)
    k = length(x)
    if(kmax == Inf) init = initial.snpmle(x, init)
    else init = initial.snpmle(x, init, kmax=kmax)
    beta = init$beta
    nb = length(beta)
    mix = init$mix
    ll1 = -Inf
    convergence = 1
    for(i in 1:maxit) {
        #cat(paste("\n", "Iteration:", i, "\n", sep = " "))
        l = logd(x, beta, mix$pt, which=c(1,0,0,0))$ld
        ma = apply(l, 1, max)
        dmix = drop(exp(l - ma) %*% mix$pr) + 1e-100
        switch(plot,
               "gradient" = plotgrad(x, beta, mix, ma,
                pch=19, order=plotorder),
               "prob" = plot(x, mix, beta),
               "dden" = plot(mix) )
        #if(plot == "gradient") points(x$mi, rep(0,length(x$mi)),
         #  pch="|", cex=.5)
        if(length(mix$pt) < kmax) {
            gridpoints = seq(0.1, 20, length = grid)
            g = maxgrad(x, beta, dmix, ma, grid=gridpoints, tol=-Inf)
            # g = maxgrad2(x, beta, dmix, ma, mix$pt, tol=-Inf)
            if(plot=="gradient") points(g$pt, g$grad, pch=20, col="blue")
            gradient = max(g$grad)
            kpt = min(kmax - length(mix$pt), length(g$pt))
            jpt = order(g$grad, decreasing=TRUE)
            mix = dden(c(mix$pt,g$pt[jpt][1:kpt]), c(mix$pr,rep(0,kpt)))
        }
        lpt = logd(x, beta, mix$pt, which=c(1,0,0,0))$ld
        dpt = pmin(exp(lpt - ma), 1e100)
        a = cbind(dpt/dmix - drop(rep(2,k)))
        r = nnls(rbind(a, rep(1,length(mix$pt))), c(rep(0,nrow(a)),1))
        sol = r$x / sum(r$x)
        r = lsch(mix, beta, dden(mix$pt,sol), beta, x, which=c(1,0,0))
        mix = collapse.snpmle(r$mix, beta, x)
        r = switch(model,
        spmle = bfgs(mix, beta, x, which=c(1,1,1)),
        npmle = bfgs(mix, beta, x, which=c(1,1,0)))

        ## Scale the error distribution and the beta
        sc <- sqrt(sum(r$mix$pr * r$mix$pt^2))
        #print("Before Scale:")
        #print.snpmle(verb, x, r$mix, r$beta, gradient)
        #print(paste("Convergence Code: ", r$conv, sep = ""))
        #print(r$beta)
        #print(r$mix)
        #print(r$num.iter)
        if(abs(sc - 1) > tol){
            new.sc <- ifelse(valid(x, c(r$beta[1] * sc,
                                        r$beta[2],
                                        r$beta[3] * sc,
                                        r$beta[4] * sc), r$mix),
                             sc, (1/(r$beta[2] + r$beta[3]) - 1e-7))
            #print(new.sc)
            #print(r$mix)
            r$mix$pt <- r$mix$pt/new.sc
            #print(r$mix)
            r$beta[c(1, 3)] <- r$beta[c(1, 3)] * new.sc^2
            r$beta[4] <- r$beta[4] * new.sc
            r$conv <- 4
            #next
        }
        #print("After Scale:")
        #print.snpmle(verb, x, r$mix, r$beta, gradient)
        #print(sc)
        #print(sqrt(sum(r$mix$pr * r$mix$pt^2)))
        #print(r$mix)
#        if(r$conv == 3) {convergence = r$conv; break}
        beta = r$beta
        mix = r$mix
        #print(r$ll)
        #print(ll1)
        if(is.null(llt))
        {if(r$ll >= ll1 && r$ll <= ll1 + tol) {convergence = 0; break}}
        else {
            if(r$ll >= llt) {convergence = 0; break}
            else if(r$ll >= ll1 && r$ll <= ll1 + 1e-16) {convergence = 1; break}
        }
        ll1 = r$ll
        print.snpmle(verb, x, mix, beta, gradient)
    }
    list(mix=mix, beta=beta, num.iterations=i,
         ll=r$ll, grad=r$grad,
                                        # max.gradient=gradient,
         convergence=convergence)
}
