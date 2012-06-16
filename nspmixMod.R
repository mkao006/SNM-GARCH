########################################################################
## This script contains function of the nspmix packge which are modified
########################################################################

## NOTE (Michael): I need to up date the xt with every iteration of the
##                 bfgs rather than updating it at the end of the bfgs.
##
## NOTE (Michael): The calculation of mean for the distribution is not
##                 included in every position and thus the
##                 log-likelihood calculated is incorrect for many
##                 places in the algorithm.

## This is the cnmms function with the scale and the function
## demeaned, a class is also added to the result for other function to
## work.

## For all initialisation
cnmms <- function (x, init = NULL, maxit = 1000,
                   model = c("spmle", "npmle"), tol = 1e-8,
                   grid = 100, kmax = Inf,
                   plot = c("null", "gradient", "prob"), verb = 0){
  initMethod = attr(x, "initialisation")
  
  T <- length(x)
  .StartTime = Sys.time()
  plot = match.arg(plot)
  model = match.arg(model)
  if (kmax == Inf)
    init = initial.snpmle(x, init)
  else init = initial.snpmle(x, init, kmax = kmax)
  beta = init$beta
  if (is.null(beta) || is.na(beta))
    model = "npmle"
  nb = length(beta)
  mix = init$mix
  gradient = "Not computed"

  ## Set parameter for different parametisation
  if(initMethod == "Smooth"){
    sigma1 = beta[5]
    k = attr(x, "smoothWindow")
  } else {
    sigma1 = attr(x, "sigma1")
    k = 1
  }
  
  switch(plot, gradient = plotgrad(x, beta, mix, pch = 1),
         prob = plot(x, mix, beta))
  ll1 = -Inf
  
  convergence = 1
  w = weights(x)
  wr = sqrt(w)
  for (i in 1:maxit) {
    s = cond.sd(x, beta)
    mu = c((sqrt(2 * outer(s^2, mix$pt^2)/pi) *
      (beta[4]^2 - 1/beta[4]^2)/(beta[4] + 1/beta[4])) %*% mix$pr)
    l = logd(x - mu, beta, mix$pt, which = c(1, 0, 0))$ld
    ma = apply(l, 1, max)
    dmix = drop(exp(l - ma) %*% mix$pr) + 1e-100
    
    ## Finding and adding a new support point
    if (length(mix$pt) < kmax) {
      gp = gridpoints(x, beta, grid)
      g = maxgrad(x, beta, dmix, ma, grid = gp, tol = -Inf)
      gradient = max(g$grad)
      kpt = min(kmax - length(mix$pt), length(g$pt))
      jpt = order(g$grad, decreasing = TRUE)
      mix = disc(c(mix$pt, g$pt[jpt][1:kpt]), c(mix$pr,
        rep(0, kpt)))
    }
    
    lpt = logd(x - mu, beta, mix$pt, which = c(1, 0, 0))$ld
    dpt = pmin(exp(lpt - ma), 1e+100)
    a = wr * (dpt/dmix - 2)
    grad.support = colSums(w * (dpt/dmix - 1))
    
    ## find the mixing proportion
    r = nnls(rbind(a, rep(1, length(mix$pt))), c(rep(0, nrow(a)), 1))
    sol = r$x/sum(r$x)
    r = lsch(mix, beta, disc(mix$pt, sol), beta, x, which = c(1, 0, 0))
    mix = collapse.snpmle(r$mix, beta, x)
    if (max(grad.support) < 1e+05) {
      r = switch(model, spmle = bfgs(mix, beta, x, which = c(1, 1, 1)),
        npmle = bfgs(mix, beta, x, which = c(1, 1, 0)))
      beta = r$beta
      mix = collapse.snpmle(r$mix, beta, x)
    }
    
    
    if(initMethod == "BackFilter"){
      ## Update the initial values.
      updateSigma1 <- as.numeric(sqrt(abs((s[2]^2 - beta[1] -
                                           beta[3] * x[1]^2)/beta[2])))
      if(abs(updateSigma1 - attr(x, "sigma1")) > tol){
        attr(x, "sigma1") <- updateSigma1
        next
      }
    }
    
    switch(plot, gradient = plotgrad(x, beta, mix, pch = 1),
           prob = plot(x, mix, beta))
    
    if (r$ll >= ll1 - tol * abs(r$ll) &&
        r$ll <= ll1 + tol * abs(r$ll)) {
      convergence = 0
      break
    }
    ll1 = r$ll
    print.snpmle(verb, x, mix, beta, gradient)
  }
  
  m = length(mix$pt)
  ## It seems that k is undefined and should be replaced by m.
  if (m < length(r$mix$pt)) {
    d = dll.snpmle(x, mix, beta, which = c(0, 1, 1, 1))
    grad = c(d$dp, d$dt, d$db)
    names(grad) = c(paste("pr", 1:m, sep = "."),
           paste("pt", 1:m, sep = "."),
           paste("beta", 1:length(beta), sep = "."))
  }
  else grad = r$grad
  grad[1:m] = grad[1:m] - sum(rep(w, len = length(x)))
  
  ## Start of scaling
  sc = sqrt(mssmnorm(varmix = mix, xi = beta[4], moments = 2))
  if(abs(sc - 1) > tol){
    if(verb != 0){
      print("Likelihood Before Scaling:")
      print.snpmle(verb, x, mix, beta, gradient)
      print(paste("Scaled by :", sc, sep = ""))
      mix$pt <- mix$pt/sc
      beta[c(1, 3)] <- beta[c(1, 3)] * sc^2
      if(initMethod == "Smooth"){
        beta[5] <- beta[5] * sc
      } else {
        attr(x, "sigma1") <- attr(x, "sigma1") * sc
      }
      print("Likelihood After Scale:")
      print.snpmle(verb, x, mix, beta, gradient)
    }
  }
  ## End of scaling
  
  .EndTime = Sys.time()
  Time = .EndTime - .StartTime
  
  ## Compute the conditional standard deviation
  if(initMethod == "Smooth"){
    sigma1 = beta[5]
  } else {
    sigma1 = attr(x, "sigma1")
  }
  
  betaSum <- as.numeric(filter(x[c(k:(T - 1))]^2, beta[2], "recursive"))
  sigma.t <-
    sqrt(c(rep(sigma1^2, k),
           beta[1] * (1 - cumprod(rep.int(beta[2], T - k)))/
           (1 - beta[2]) +
           cumprod(rep.int(beta[2], T - k)) * sigma1^2 +
           beta[3] * betaSum))
  
  result <- list(mix = mix, beta = beta, num.iterations = i,
                 ll = attr(mix,"ll")[1], grad = grad,
                 convergence = convergence,
                 sigma.t = sigma.t)
  attr(result, "class") <- "mgarch"
  attr(result, "time") <- Time
  result
}

## For smoothing initialisation
## cnmms <- function (x, init = NULL, maxit = 1000,
##                    model = c("spmle", "npmle"), tol = 1e-5,
##                    grid = 100, kmax = Inf,
##                    plot = c("null", "gradient", "prob"), verb = 0){
##     .StartTime = Sys.time()
##     plot = match.arg(plot)
##     model = match.arg(model)
##     if (kmax == Inf)
##         init = initial.snpmle(x, init)
##     else init = initial.snpmle(x, init, kmax = kmax)
##     beta = init$beta
##     if (is.null(beta) || is.na(beta))
##         model = "npmle"
##     nb = length(beta)
##     mix = init$mix
##     gradient = "Not computed"
##     switch(plot, gradient = plotgrad(x, beta, mix, pch = 1),
##            prob = plot(x, mix, beta))
##     ll1 = -Inf
##
##     convergence = 1
##     w = weights(x)
##     wr = sqrt(w)
##     for (i in 1:maxit) {
##         s = cond.sd(x, beta)
##         mu = c((sqrt(2 * outer(s^2, mix$pt^2)/pi) *
##           (beta[5]^2 - 1/beta[5]^2)/(beta[5] + 1/beta[5])) %*% mix$pr)
##         l = logd(x - mu, beta, mix$pt, which = c(1, 0, 0))$ld
##         ma = apply(l, 1, max)
##         dmix = drop(exp(l - ma) %*% mix$pr) + 1e-100
##
##         ## Finding and adding a new support point
##         if (length(mix$pt) < kmax) {
##             gp = gridpoints(x, beta, grid)
##             g = maxgrad(x, beta, dmix, ma, grid = gp, tol = -Inf)
##             gradient = max(g$grad)
##             kpt = min(kmax - length(mix$pt), length(g$pt))
##             jpt = order(g$grad, decreasing = TRUE)
##             mix = disc(c(mix$pt, g$pt[jpt][1:kpt]), c(mix$pr,
##                        rep(0, kpt)))
##         }
##
##         lpt = logd(x - mu, beta, mix$pt, which = c(1, 0, 0))$ld
##         dpt = pmin(exp(lpt - ma), 1e+100)
##         a = wr * (dpt/dmix - 2)
##         grad.support = colSums(w * (dpt/dmix - 1))
##
##         ## find the mixing proportion
##         r = nnls(rbind(a, rep(1, length(mix$pt))), c(rep(0, nrow(a)), 1))
##         sol = r$x/sum(r$x)
##         r = lsch(mix, beta, disc(mix$pt, sol), beta, x, which = c(1, 0, 0))
##         mix = collapse.snpmle(r$mix, beta, x)
##         if (max(grad.support) < 1e+05) {
##             r = switch(model, spmle = bfgs(mix, beta, x, which = c(1, 1, 1)),
##                        npmle = bfgs(mix, beta, x, which = c(1, 1, 0)))
##             beta = r$beta
##             mix = collapse.snpmle(r$mix, beta, x)
##         }
##
##         print(beta)
##         switch(plot, gradient = plotgrad(x, beta, mix, pch = 1),
##                prob = plot(x, mix, beta))
##
##         if (r$ll >= ll1 - tol * abs(r$ll) &&
##             r$ll <= ll1 + tol * abs(r$ll)) {
##             convergence = 0
##             break
##         }
##         ll1 = r$ll
##         print.snpmle(verb, x, mix, beta, gradient)
##     }
##
##     m = length(mix$pt)
##     ## It seems that k is undefined and should be replaced by m.
##     if (m < length(r$mix$pt)) {
##       d = dll.snpmle(x, mix, beta, which = c(0, 1, 1, 1))
##       grad = c(d$dp, d$dt, d$db)
##       names(grad) = c(paste("pr", 1:m, sep = "."), paste("pt",
##             1:m, sep = "."), paste("beta", 1:length(beta), sep = "."))
##     }
##     else grad = r$grad
##     grad[1:m] = grad[1:m] - sum(rep(w, len = length(x)))
##
##     ## Start of scaling
##     print(mix)
##     print(beta)
##     sc = sqrt(mssmnorm(varmix = mix, xi = beta[5], moments = 2))
##     if(abs(sc - 1) > tol){
##       if(verb != 0){
##         print("Likelihood Before Scaling:")
##         print.snpmle(verb, x, mix, beta, gradient)
##         print(paste("Scaled by :", sc, sep = ""))
##         mix$pt <- mix$pt/sc
##         beta[c(1, 3)] <- beta[c(1, 3)] * sc^2
##         beta[4] <- beta[4] * sc
##         print("Likelihood After Scale:")
##         print.snpmle(verb, x, mix, beta, gradient)
##       }
##     }
##     ## End of scaling
##
##     .EndTime = Sys.time()
##     Time = .EndTime - .StartTime
##
##     k = attr(x, "smoothWindow")
##     T <- length(x)
##     betaSum <- as.numeric(filter(x[c(k:(T - 1))]^2, beta[2], "recursive"))
##     sigma.t <-
##         sqrt(c(rep(beta[4]^2, k),
##                beta[1] * (1 - cumprod(rep.int(beta[2], T - k)))/
##                (1 - beta[2]) +
##                cumprod(rep.int(beta[2], T - k)) * beta[4]^2 +
##                beta[3] * betaSum))
##
##      result <- list(mix = mix, beta = beta, num.iterations = i,
##                    ll = attr(mix,"ll")[1], grad = grad,
##                    convergence = convergence,
##                    sigma.t = sigma.t)
##     attr(result, "class") <- "mgarch"
##     attr(result, "time") <- Time
##     result
## }


## For back estimation and unconditional variance
## cnmms <- function (x, init = NULL, maxit = 1000,
##                    model = c("spmle", "npmle"), tol = 1e-5,
##                    grid = 100, kmax = Inf,
##                    plot = c("null", "gradient", "prob"), verb = 0){
##     .StartTime = Sys.time()
##     plot = match.arg(plot)
##     model = match.arg(model)
##     if (kmax == Inf)
##         init = initial.snpmle(x, init)
##     else init = initial.snpmle(x, init, kmax = kmax)
##     beta = init$beta
##     if (is.null(beta) || is.na(beta))
##         model = "npmle"
##     nb = length(beta)
##     mix = init$mix
##     gradient = "Not computed"
##     switch(plot, gradient = plotgrad(x, beta, mix, pch = 1),
##            prob = plot(x, mix, beta))
##     ll1 = -Inf
##
##     convergence = 1
##     w = weights(x)
##     wr = sqrt(w)
##     for (i in 1:maxit) {
##         s = cond.sd(x, beta)
##         mu = c((sqrt(2 * outer(s^2, mix$pt^2)/pi) *
##           (beta[4]^2 - 1/beta[4]^2)/(beta[4] + 1/beta[4])) %*% mix$pr)
##         l = logd(x - mu, beta, mix$pt, which = c(1, 0, 0))$ld
##         ma = apply(l, 1, max)
##         dmix = drop(exp(l - ma) %*% mix$pr) + 1e-100
##
##         ## Finding and adding a new support point
##         if (length(mix$pt) < kmax) {
##             gp = gridpoints(x, beta, grid)
##             g = maxgrad(x, beta, dmix, ma, grid = gp, tol = -Inf)
##             gradient = max(g$grad)
##             kpt = min(kmax - length(mix$pt), length(g$pt))
##             jpt = order(g$grad, decreasing = TRUE)
##             mix = disc(c(mix$pt, g$pt[jpt][1:kpt]), c(mix$pr,
##                        rep(0, kpt)))
##         }
##
##         lpt = logd(x - mu, beta, mix$pt, which = c(1, 0, 0))$ld
##         dpt = pmin(exp(lpt - ma), 1e+100)
##         a = wr * (dpt/dmix - 2)
##         grad.support = colSums(w * (dpt/dmix - 1))
##
##         ## find the mixing proportion
##         r = nnls(rbind(a, rep(1, length(mix$pt))), c(rep(0, nrow(a)), 1))
##         sol = r$x/sum(r$x)
##         r = lsch(mix, beta, disc(mix$pt, sol), beta, x, which = c(1, 0, 0))
##         mix = collapse.snpmle(r$mix, beta, x)
##         if (max(grad.support) < 1e+05) {
##             r = switch(model, spmle = bfgs(mix, beta, x, which = c(1, 1, 1)),
##                        npmle = bfgs(mix, beta, x, which = c(1, 1, 0)))
##             beta = r$beta
##             mix = collapse.snpmle(r$mix, beta, x)
##         }
##
##         switch(plot, gradient = plotgrad(x, beta, mix, pch = 1),
##                prob = plot(x, mix, beta))
##
##         ## Update the initial values.
##         sigma1 <- as.numeric(sqrt(abs((s[2]^2 - beta[1] -
##                                        beta[3] * x[1]^2)/beta[2])))
##         if(abs(sigma1 - attr(x, "sigma1")) > tol){
##           attr(x, "sigma1") <- sigma1
##           next
##         }
##
##         if (## max(r$grad) < tol &&
##             r$ll >= ll1 - tol * abs(r$ll) &&
##             r$ll <= ll1 + tol * abs(r$ll)) {
##             convergence = 0
##             break
##         }
##         ll1 = r$ll
##         print.snpmle(verb, x, mix, beta, gradient)
##     }
##
##     m = length(mix$pt)
##     ## It seems that k is undefined and should be replaced by m.
##     if (m < length(r$mix$pt)) {
##       d = dll.snpmle(x, mix, beta, which = c(0, 1, 1, 1))
##       grad = c(d$dp, d$dt, d$db)
##       names(grad) = c(paste("pr", 1:m, sep = "."), paste("pt",
##             1:m, sep = "."), paste("beta", 1:length(beta), sep = "."))
##     }
##
##     else grad = r$grad
##     grad[1:m] = grad[1:m] - sum(rep(w, len = length(x)))
##
##     ## Start of scaling
##     sc = sqrt(mssmnorm(varmix = mix, xi = beta[4], moments = 2))
##     if(abs(sc - 1) > tol){
##       if(verb != 0){
##         print("Likelihood Before Scaling:")
##         print.snpmle(verb, x, mix, beta, gradient)
##         print(paste("Scaled by :", sc, sep = ""))
##         mix$pt <- mix$pt/sc
##         beta[c(1, 3)] <- beta[c(1, 3)] * sc^2
##         attr(x, "sigma1") <- attr(x, "sigma1") * sc
##         print("Likelihood After Scale:")
##         print.snpmle(verb, x, mix, beta, gradient)
##       }
##     }
##     ## End of scaling
##
##     .EndTime = Sys.time()
##     Time = .EndTime - .StartTime
##
##     T <- length(x)
##     betaSum <- as.numeric(filter(x[c(1:(T - 1))]^2, beta[2], "recursive"))
##     sigma.t <-
##         sqrt(c(sigma1^2,
##                beta[1] * (1 - cumprod(rep.int(beta[2], T - 1)))/
##                (1 - beta[2]) +
##                cumprod(rep.int(beta[2], T - 1)) * sigma1^2 +
##                beta[3] * betaSum))
##
##      result <- list(mix = mix, beta = beta, num.iterations = i,
##                    ll = attr(mix,"ll")[1], grad = grad,
##                    convergence = convergence,
##                    sigma.t = sigma.t)
##     attr(result, "class") <- "mgarch"
##     attr(result, "time") <- Time
##     result
## }

lsch <- function (mix1, beta1, mix2, beta2, x, maxit = 100,
                  which = c(1, 1, 1), brkt = FALSE){
    k = length(mix1$pt)
    convergence = 1
    dl1 = dll.snpmle(x, mix1, beta1, which = c(1, which))
    lla = ll1 = dl1$ll
    names.grad = c(if (which[1]) paste("pr", 1:k, sep = ".") else NULL,
        if (which[2]) paste("pt", 1:k, sep = ".") else NULL,
        if (which[3]) paste("beta", 1:length(beta1), sep = ".") else NULL)
    grad1 = c(if (which[1]) dl1$dp else NULL, if (which[2]) dl1$dt else NULL,
        if (which[3]) dl1$db else NULL)
    names(grad1) = names.grad
    d1 = c(if (which[1]) mix2$pr - mix1$pr else NULL, if (which[2]) mix2$pt -
        mix1$pt else NULL, if (which[3]) beta2 - beta1 else NULL)
    d1.norm = sqrt(sum(d1 * d1))
    s = d1/d1.norm
    g1d1 = sum(grad1 * d1)
    dla = g1s = g1d1/d1.norm
    if (d1.norm == 0 || g1s <= 0) {
        return(list(mix = mix1, beta = beta1, grad = grad1, ll = ll1,
            convergence = 3))
    }
    ## print(mix2)
    ## print(beta2)
    a = 0
    b = 1
    if (which[1] && any(mix2$pr == 0))
        brkt = FALSE
    for (i in 1:maxit) {
        for (j in 1:3000) {
            m = disc((1 - b) * mix1$pt + b * mix2$pt, (1 - b) *
                mix1$pr + b * mix2$pr)
            beta = if (is.null(beta1))
                NULL
            else (1 - b) * beta1 + b * beta2
            if (valid.snpmle(x, beta, m))
                break
            brkt = FALSE
            b = 0.5 * a + 0.5 * b
        }
        if (j == 3000)
            warning("Can not produce valid interior point in lsch()")
        ## print(beta)
        dl = dll.snpmle(x, m, beta, which = c(1, which))
        ll = dl$ll
        grad = c(if (which[1]) dl$dp else NULL, if (which[2]) dl$dt else NULL,
            if (which[3]) dl$db else NULL)
        gs = sum(grad * s)
        zzz = c(gs, g1s, ll, ll1, g1d1, b)
        ## print(zzz)
        if(any( is.na(zzz) )) {print(m); print(beta); print(dl); print(zzz);
                               stop("=== here ===")}

        if (brkt && gs > g1s * 0.5 && ll >= ll1 + g1d1 * b * 0.33) {
          a = b
          b = 2 * b
          lla = ll
          dla = gs
        }
        else break
    }
    if (i == maxit)
        brkt = FALSE
    alpha = b
    llb = ll
    dlb = gs
    for (i in 1:maxit) {
        g1d = g1d1 * alpha
        if (ll >= ll1 - 1e-15 * abs(ll1) && g1d <= 1e-15 * abs(ll)) {
            convergence = 2
            break
        }
        if (brkt) {
            if (ll >= ll1 + g1d * 0.33 && abs(gs) <= g1s * 0.5) {
                convergence = 0
                break
            }
            if (ll >= ll1 + g1d * 0.33 && gs > g1s * 0.5) {
                a = alpha
                lla = ll
                dla = gs
            }
            else {
                b = alpha
                llb = ll
                dlb = gs
            }
        }
        else {
            if (ll >= ll1 + g1d * 0.33) {
                convergence = 0
                break
            }
            else {
                b = alpha
                llb = ll
                dlb = gs
            }
        }
        alpha = (a + b) * 0.5
        m = disc((1 - alpha) * mix1$pt + alpha * mix2$pt, (1 -
            alpha) * mix1$pr + alpha * mix2$pr)
        beta = if (is.null(beta1))
            NULL
        else (1 - alpha) * beta1 + alpha * beta2
        dl = dll.snpmle(x, m, beta, which = c(1, which))
        ll = dl$ll
        grad = c(if (which[1]) dl$dp else NULL, if (which[2]) dl$dt else NULL,
            if (which[3]) dl$db else NULL)
        gs = sum(grad * s)
    }
    names(grad) = names.grad
    beta = if (is.null(beta1))
        NULL
    else (1 - alpha) * beta1 + alpha * beta2
    # print(alpha)
    list(mix = disc((1 - alpha) * mix1$pt + alpha * mix2$pt,
        (1 - alpha) * mix1$pr + alpha * mix2$pr), beta = beta,
        grad = grad, ll = ll, convergence = convergence, num.iterations = i)
}





