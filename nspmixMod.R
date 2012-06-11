########################################################################
## This script contains function of the nspmix packge which are modified
########################################################################

library(nspmix)

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
cnmms <- function (x, init = NULL, maxit = 1000,
                   model = c("spmle", "npmle"), tol = 1e-5,
                   grid = 100, kmax = Inf,
                   plot = c("null", "gradient", "prob"), verb = 0){
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
        ## print(mu)
        ## print(mssmnorm(varmix = mix, xi = beta[4], moments = 1))
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

        switch(plot, gradient = plotgrad(x, beta, mix, pch = 1),
               prob = plot(x, mix, beta))

        ## Update the initial values.
        attr(x, "sigma1") <- as.numeric(sqrt(abs((s[2]^2 - beta[1] -
                                    beta[3] * x[1]^2)/beta[2])))
        
        ll2 <- logLik.snpmle(x, beta, mix, attr = FALSE)
        print(r$ll)
        print(ll2)
        if(abs(r$ll - ll2) > tol) next
        

        ## Recompute the likelihood given the updated beta[4], the
        ## likelihood will no necessary increase since the solution is
        ## restricted rather than freely estimated. Need to think of a
        ## way to avoid this!
        ##
        ## r$ll = logLik.snpmle(x, beta, mix)[1]
        
        if (## max(r$grad) < tol &&
            r$ll >= ll1 - tol * abs(r$ll) &&
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
      names(grad) = c(paste("pr", 1:m, sep = "."), paste("pt",
            1:m, sep = "."), paste("beta", 1:length(beta), sep = "."))
    }
    ## if (m < length(r$mix$pt)) {
    ##     d = dll.snpmle(x, mix, beta, which = c(0, 1, 1, 1))
    ##     grad = c(d$dp, d$dt, d$db)
    ##     names(grad) = c(paste("pr", 1:k, sep = "."), paste("pt",
    ##         1:k, sep = "."), paste("beta", 1:length(beta), sep = "."))
    ## }
    else grad = r$grad
    grad[1:m] = grad[1:m] - sum(rep(w, len = length(x)))
    
    ## ## Start of scaling
    sc = sqrt(mssmnorm(varmix = mix, xi = beta[4], moments = 2))
    if(abs(sc - 1) > tol){
      if(verb != 0){
        print("Likelihood Before Scaling:")
        print.snpmle(verb, x, mix, beta, gradient)
        print(paste("Scaled by :", sc, sep = ""))
        mix$pt <- mix$pt/sc
        beta[c(1, 3)] <- beta[c(1, 3)] * sc^2
        attr(x, "sigma1") <- attr(x, "sigma1") * sc
        print("Likelihood After Scale:")
        print.snpmle(verb, x, mix, beta, gradient)
      }
    }
    ## ## End of scaling
    
    .EndTime = Sys.time()
    Time = .EndTime - .StartTime
    sigma.t = cond.sd(x, beta)
    result <- list(mix = mix, beta = beta, num.iterations = i,
                   ll = attr(mix,"ll")[1], grad = grad,
                   convergence = convergence,
                   sigma.t = sigma.t)
    attr(result, "class") <- "mgarch"
    attr(result, "time") <- Time
    result
}


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





bfgs <- function (mix, beta, x, tol = 1e-15, maxit = 1000, which = c(1, 
    1, 1), D = NULL) {
    w = weights(x)
    k1 = if (which[1]) 
        length(mix$pr)
    else 0
    k2 = if (which[2]) 
        length(mix$pt)
    else 0
    k3 = if (which[3]) 
        length(beta)
    else 0
    if (sum(which) == 0) 
        stop("No parameter specified for updating in bfgs()")
    dl = dll.snpmle(x, mix, beta, which = c(1, which), ind = TRUE)
    ll = sum(w * dl$ll)
    grad = c(if (which[1]) colSums(w * dl$dp) else NULL, if (which[2]) colSums(w * 
        dl$dt) else NULL, if (which[3]) colSums(w * dl$db) else NULL)
    r = list(mix = mix, beta = beta, ll = ll, grad = grad, convergence = 1)
    par = c(if (which[1]) r$mix$pr else NULL, if (which[2]) r$mix$pt else NULL, 
        if (which[3]) r$beta else NULL)
    if (is.null(D)) 
        D = diag(-1, nrow = k1 + k2 + k3)
    else if (nrow(D) != k1 + k2 + k3) 
        stop("D provided has incompatible dimensions")
    bs = suppspace(x, beta)
    for (i in 1:maxit) {
        old.r = r
        old.par = par
        d1 = -drop(D %*% r$grad)
        if (which[1]) {
            D1 = D[, 1:k1, drop = FALSE]
            D11 = D1[1:k1, , drop = FALSE]
            lambda = -sum(r$grad * rowSums(D1))/sum(D11)
            d1 = d1 - lambda * rowSums(D1)
        }
        par2 = par + d1
        if (which[2]) {
            theta1 = par[k1 + 1:k2]
            theta2 = par2[k1 + 1:k2]
            dtheta = d1[k1 + 1:k2]
            if (any(theta2 < bs[1], theta2 > bs[2])) {
                out = pmax(c(bs[1] - theta2, theta2 - bs[2]), 
                  0)
                j = out > 0 & (theta1 == bs[1] | theta1 == bs[2])
                if (any(j)) {
                  jj = (which(j) - 1)%%k2 + 1
                  D[k1 + jj, ] = D[, k1 + jj] = 0
                  d1 = -drop(D %*% r$grad)
                  par2 = par + d1
                }
                else {
                  ratio = out/abs(dtheta)
                  ratio[dtheta == 0] = 0
                  jmax = which.max(ratio)
                  alpha = 1 - ratio[jmax]
                  d1 = alpha * d1
                  par2 = par + d1
                  if (jmax <= k2) 
                    par2[k1 + jmax] = bs[1]
                  else par2[k1 + jmax - k2] = bs[2]
                }
            }
        }
        if (which[1]) {
            if (any(par2[1:k1] < 0)) {
                step = d1[1:k1]
                ratio = pmax(-par2[1:k1], 0)/abs(step)
                jmax = which.max(ratio)
                alpha = 1 - ratio[jmax]
                d1 = alpha * d1
                par2 = par + d1
                par2[jmax] = 0
            }
        }
        mix2 = disc(if (which[2]) 
            par2[k1 + 1:k2]
        else r$mix$pt, if (which[1]) 
            par2[1:k1]
        else r$mix$pr)
        beta2 = if (which[3]) 
            par2[k1 + k2 + 1:k3]
        else r$beta
        r = lsch(r$mix, r$beta, mix2, beta2, x, which = which, 
            brkt = TRUE)
        if (any(r$mix$pr == 0)) {
            j = r$mix$pr == 0
            r$mix = disc(r$mix$pt[!j], r$mix$pr[!j])
            j2 = which(j)
            if (which[2]) 
                j2 = c(j2, k1 + j2)
            D = D[-j2, -j2]
            return(bfgs(r$mix, r$beta, x, tol, maxit, which, 
                D = D))
        }
        if (r$conv != 0) 
            break
        if (r$ll >= old.r$ll - tol * abs(r$ll) && r$ll <= old.r$ll + 
            tol * abs(r$ll)) {
            convergence = 0
            break
        }
        g = r$grad - old.r$grad
        par = c(if (which[1]) r$mix$pr else NULL, if (which[2]) r$mix$pt else NULL, 
            if (which[3]) r$beta else NULL)
        d = par - old.par
        dg = sum(d * g)
        if (dg < 0) {
            nd = length(d)
            dod = d * rep(d, each = nd)
            dim(dod) = c(nd, nd)
            dog = d * rep(g, each = nd)
            dim(dog) = c(nd, nd)
            dogD = dog %*% D
            D = D + (1 + drop(t(g) %*% D %*% g)/dg) * dod/dg - 
                (dogD + t(dogD))/dg
        }
    }
    r$num.iterations = i
    r
}


## initial.snpmle <- function (x, init = NULL, kmax = NULL, grid = 100) {
##     if (is.null(init) || is.null(init$mix)) {
##         init = initial(x, init$beta, init$mix, kmax = kmax)
##         init$ll = logLik.snpmle(x, init$beta, init$mix)
##     }
##     else {
##         init0 = initial(x, init$beta, kmax = kmax)
##         init0$ll = logLik.snpmle(x, init0$beta, init0$mix)
##         ll = logLik.snpmle(x, init0$beta, init$mix)
##         if (ll > init0$ll) 
##             init = list(beta = init0$beta, mix = init$mix, ll = ll)
##         else init = init0
##     }
##     if (!valid.snpmle(x, init$beta, init$mix)) 
##         stop("Invalid initial values!")
##     init
## }


## ## Finds the support point set which satisfies the NPMLE and generalised
## ## equilibrium theorem
## maxgrad <- function (x, beta, dmix, ma, grid = 100, tol = -Inf, maxit = 100) {
##     if (length(grid) == 1) 
##         grid = gridpoints(x, beta, grid)
##     np = length(grid)
##     dg = grad(x, grid, beta, dmix, ma, order = 1)$d1
##     jmax = which(dg[-np] > 0 & dg[-1] < 0)
##     if (length(jmax) < 1) 
##         return
##     pt = (grid[jmax] + grid[jmax + 1]) * 0.5
##     left = grid[jmax]
##     right = grid[jmax + 1]
##     if (length(pt) != 0) {
##         pt.old = left
##         d1.old = grad(x, left, beta, dmix, ma, order = 1)$d1
##         d2 = rep(-1, length(pt))
##         for (i in 1:maxit) {
##             d1 = grad(x, pt, beta, dmix, ma, order = 1)$d1
##             d2t = (d1 - d1.old)/(pt - pt.old)
##             jd = !is.na(d2t) & d2t < 0
##             d2[jd] = d2t[jd]
##             left[d1 > 0] = pt[d1 > 0]
##             right[d1 < 0] = pt[d1 < 0]
##             pt.old = pt
##             d1.old = d1
##             pt = pt - d1/d2
##             j = is.na(pt) | pt < left | pt > right
##             pt[j] = (left[j] + right[j]) * 0.5
##             if (max(abs(pt - pt.old)) <= 1e-14 * diff(range(grid))) 
##                 break
##         }
##     }
##     else i = 0
##     if (dg[np] >= 0) 
##         pt = c(grid[np], pt)
##     if (dg[1] <= 0) 
##         pt = c(grid[1], pt)
##     if (length(pt) == 0) 
##         stop("no new support point found")
##     g = grad(x, pt, beta, dmix, ma, order = 0)$d0
##     names(pt) = names(g) = NULL
##     j = g >= tol
##     list(pt = pt[j], grad = g[j], num.iterations = i)
## }


## ## Computes the gradient o each mixture, returns the same length as the
## ## mixture.
## grad <- function (x, pt, beta, dmix, ma, order = 0) {
##     w = weights(x)
##     if (is.disc(dmix)) {
##       l = logd(x, beta, dmix$pt, which = c(1, 0, 0))$ld
##       ma = apply(l, 1, max)
##       dmix = drop(exp(l - ma) %*% dmix$pr) + 1e-100
##     }
##     g = vector("list", length(order))
##     names(g) = paste("d", order, sep = "")
##     which = c(1, 0, 0)
##     if (any(order >= 1)) 
##         which[3:max(order + 2)] = 1
##     s = cond.sd(x, beta)
##     dl = logd(x, beta, pt, which = which)
##     ws = w * pmin(exp(dl$ld - ma), 1e+100)/dmix
##     if (0 %in% order) 
##         g$d0 = colSums(ws - w)
##     if (1 %in% order) 
##         g$d1 = colSums(ws * dl$dt)
##     g
## }






## logLik.snpmle <- function (x, beta, mix, attr = TRUE) {
##   ld = logd(x, beta, mix$pt, which = c(1, 0, 0))$ld + rep(log(mix$pr), 
##                                each = length(x))
##   ma = apply(ld, 1, max)
##   pid = exp(ld - ma)
##   pis = rowSums(pid)
##   ll = sum(weights(x) * (log(pis) + ma))
##   if (attr) {
##     attr(ll, "p") = pid/pis
##     attr(ll, "ma") = ma
##   }
##   ll
## }



## dll.snpmle <- function (x, mix, beta, which = c(1, 0, 0, 0), ind = FALSE) {
##     w = weights(x)
##     r = list()
##     dl = logd(x, beta, mix$pt, which = c(1, which[4:3]))
##     lpt = dl$ld
##     ma = apply(lpt, 1, max)
##     if (which[1] == 1) {
##         r$ll = log(drop(exp(lpt - ma) %*% mix$pr)) + ma
##         if (!ind) 
##             r$ll = sum(w * r$ll)
##     }
##     if (sum(which[2:4]) == 0) 
##         return(r)
##     dmix = drop(exp(lpt - ma) %*% mix$pr) + 1e-100
##     dpt = pmin(exp(lpt - ma), 1e+100)
##     dp = dpt/dmix
##     if (which[2] == 1) {
##         if (ind) 
##             r$dp = dp
##         else r$dp = colSums(w * dp)
##     }
##     if (sum(which[3:4]) == 0) 
##         return(r)
##     p = dp * rep(mix$pr, each = nrow(dp))
##     if (which[3] == 1) {
##         r$dt = p * dl$dt
##         if (!ind) 
##             r$dt = colSums(w * r$dt)
##     }
##     if (which[4] == 0) 
##         return(r)
##     dl1 = dl$db
##     r$db = apply(sweep(dl1, c(1, 2), p, "*"), c(1, 3), sum)
##     if (!ind) 
##         r$db = colSums(w * r$db)
##     r
## }

