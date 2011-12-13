MySVfilter <- function(num, y, phi0, phi1, sQ, alpha, sR0, mu1, sR1){
    ## Initialise the parameters
    y <- as.matrix(y)
    like <- double(length(y))
    Q <- sQ^2
    R0 <- sR0^2
    R1 <- sR1^2
    xf <- 0     	     #  <- h_0^0
    Pf <- sQ^2/(1-phi1)     #  <- P_0^0
    Pf[Pf<0] <- 0           # make sure Pf not negative
    xp <- matrix(0,num,1)   #  <- h_t^t-1
    Pp <- matrix(0,num,1)   #  <- P_t^t-1
    pi1 <- .5    #initial mix probs
    pi0 <- .5
    fpi1 <- .5
    fpi0 <- .5
    like <- 0                  # -log(likelihood)
    for (i in 1:num){
        xp[i] <- phi1 * xf + phi0
        Pp[i] <- phi1 * Pf * phi1 + Q
        sig1 <- Pp[i] + R1     #innov var
        sig0 <- Pp[i] + R0
        k1 <- Pp[i]/sig1
        k0 <- Pp[i]/sig0
        e1 <- y[i] - xp[i] - mu1 - alpha
        e0 <- y[i] - xp[i] - alpha

        den1 <- (1/sqrt(sig1)) * exp(-0.5 * e1^2/sig1)
        den0 <- (1/sqrt(sig0)) * exp(-0.5 * e0^2/sig0)
        denom <- pi1 * den1 + pi0 * den0
        fpi1 <- pi1 * den1/denom
        fpi0 <- pi0 * den0/denom

        xf <- xp[i] + fpi1 * k1 * e1 + fpi0 * k0 * e0
        Pf <- fpi1 * (1 - k1) * Pp[i] + fpi0 * (1 - k0) * Pp[i]
        like[i] <- - 0.5 * log(pi1 * den1 + pi0 * den0)
    }
    list(xp = xp, Pp = Pp, like = like)
}

######################################################################
## Temporary Kalman filter function for mixture model
######################################################################

Kfilter0 <- function(num, y, A, mu0, Sigma0, Phi, cQ, cR){
      # NOTE: must give cholesky decomp: cQ = chol(Q), cR = chol(R)
    Q = t(cQ)%*% cQ
    R = t(cR)%*% cR
      # y is num by q  (time = row series = col)
      # A is a q by p matrix
      # R is q by q
      # mu0 is p by 1
      # Sigma0, Phi, Q are p by p
    Phi = as.matrix(Phi)
    pdim = nrow(Phi)
    y = as.matrix(y)
    qdim = ncol(y)
    xp = array(NA, dim = c(pdim, 1, num))         # xp = x_t^{t-1}
    Pp = array(NA, dim = c(pdim, pdim, num))      # Pp = P_t^{t-1}
    xf = array(NA, dim = c(pdim, 1, num))         # xf = x_t^t
    Pf = array(NA, dim = c(pdim, pdim, num))      # Pf = x_t^t
    innov = array(NA, dim = c(qdim, 1, num))      # innovations
    sig = array(NA, dim = c(qdim, qdim, num))     # innov var-cov matrix
       # initialize (because R can't count from zero)
    x00 = as.matrix(mu0, nrow = pdim, ncol = 1)
    P00 = as.matrix(Sigma0, nrow = pdim, ncol = pdim)
    xp[, , 1] = Phi%*%x00
    Pp[, , 1] = Phi%*%P00%*%t(Phi)+Q
    sigtemp = A%*%Pp[, , 1]%*%t(A)+R
    sig[, , 1] = (t(sigtemp)+sigtemp)/2
       # innov var - make sure it's symmetric
    siginv = solve(sig[, , 1])
    K = Pp[, , 1]%*%t(A)%*%siginv
    innov[, , 1] = y[1, ]-A%*%xp[, , 1]
    xf[, , 1] = xp[, , 1]+K%*%innov[, , 1]
    Pf[, , 1] = Pp[, , 1]-K%*%A%*%Pp[, , 1]
    sigmat = as.matrix(sig[, , 1], nrow = qdim, ncol = qdim)
    like  =  log(det(sigmat)) +
        t(innov[, , 1])%*%siginv%*%innov[, , 1]   # -log(likelihood)

    ########## start filter iterations ###################
    for (i in 2:num){
        if (num < 2) break
        xp[, , i] = Phi%*%xf[, , i-1]
        Pp[, , i] = Phi%*%Pf[, , i-1]%*%t(Phi)+Q
        sigtemp = A%*%Pp[, , i]%*%t(A)+R
        sig[, , i] = (t(sigtemp)+sigtemp)/2
           # innov var - make sure it's symmetric
        siginv = solve(sig[, , i])
        K = Pp[, , i]%*%t(A)%*%siginv
        innov[, , i] = y[i, ]-A%*%xp[, , i]
        xf[, , i] = xp[, , i]+K%*%innov[, , i]
        Pf[, , i] = Pp[, , i]-K%*%A%*%Pp[, , i]
        sigmat = as.matrix(sig[, , i],  nrow = qdim, ncol = qdim)
        like =  like + log(det(sigmat)) +
            t(innov[, , i])%*%siginv%*%innov[, , i]
    }
    like = 0.5*like
    list(xp = xp, Pp = Pp, xf = xf, Pf = Pf, like = like,
         innov = innov, sig = sig, Kn = K)
}
