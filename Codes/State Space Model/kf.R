######################################################################
## Attempt to re-write the kf function so that it incorporates the
## non-normal density estimation
######################################################################

kf <- function (yt, Zt, Tt, Rt, Ht, Qt, a1, P1, P1inf = 0, optcal = c(TRUE,
    TRUE, TRUE, TRUE), tol = 1e-07)
{
    if (!is.array(yt)) {
        if (!is.matrix(yt))
            yt <- array(yt, dim = c(1, length(yt)))
        else yt <- array(yt, dim = dim(yt))
    }
    p <- dim(yt)[1]
    n <- dim(yt)[2]
    m <- length(a1)
    if (is.vector(Qt))
        r <- 1
    else r <- dim(as.array(Qt))[2]
    tv <- array(0, dim = 5)
    tv[1] <- !(is.na(dim(as.array(Tt))[3]) || dim(as.array(Tt))[3] ==
        1)
    tv[2] <- !(is.na(dim(as.array(Rt))[3]) || dim(as.array(Rt))[3] ==
        1)
    tv[3] <- !(is.na(dim(as.array(Qt))[3]) || dim(as.array(Qt))[3] ==
        1)
    tv[4] <- !(is.na(dim(as.array(Ht))[3]) || dim(as.array(Ht))[3] ==
        1)
    tv[5] <- !(is.na(dim(as.array(Zt))[3]) || dim(as.array(Zt))[3] ==
        1)
    ymiss <- is.na(yt)
    ydimt <- array(p, dim = n)
    tv[5] <- max(tv[4], tv[5])
    if (sum(ymiss) > 0) {
        tv[4] <- tv[5] <- 1
    }
    H <- array(Ht, c(p, p, (n - 1) * tv[4] + 1))
    Z <- array(Zt, dim = c(p, m, (n - 1) * tv[5] + 1))
    y <- yt
    if (sum(ymiss) > 0) {
        for (i in 1:n) {
            ydimt[i] <- sum(!ymiss[1:p, i])
            if (ydimt[i] != p && ydimt[i] != 0) {
                y[1:ydimt[i], i] <- yt[!ymiss[, i], i]
                H[1:ydimt[i], 1:ydimt[i], i] <- H[!ymiss[, i],
                  !ymiss[, i], i]
                Z[1:ydimt[i], , i] <- Z[!ymiss[, i], , i]
            }
        }
    }
    at <- array(0, dim = c(m, n + 1))
    Pt <- Pinf <- Pstar <- array(0, dim = c(m, m, n + 1))
    Kt <- Ktuni <- Kinf <- Kstar <- Kinfuni <- Kstaruni <- array(0,
        dim = c(m, p, n))
    vt <- vtuni <- Ftuni <- Finfuni <- Fstaruni <- array(0, dim = c(p,
        n))
    Ft <- Finf <- Fstar <- array(0, dim = c(p, p, n))
    Lt <- Linf <- Lstar <- array(0, dim = c(m, m, n))
    Pinf[, , 1] <- P1inf
    lik <- j <- d <- 0
    info <- rep(0, 4)
    storage.mode(d) <- storage.mode(j) <- storage.mode(p) <- storage.mode(m) <- storage.mode(r) <- storage.mode(n) <- storage.mode(tv) <- storage.mode(info) <- storage.mode(optcal) <- storage.mode(ydimt) <- storage.mode(j) <- "integer"
    kfout <- NULL
    kfout <- .Fortran("kf", PACKAGE = "KFAS", NAOK = TRUE, yt = y,
        ydimt = ydimt, tv = tv, Zt = Z, Tt = array(Tt, c(m, m,
            (n - 1) * tv[1] + 1)), Rt = array(Rt, c(m, r, (n -
            1) * tv[2] + 1)), Ht = H, Qt = array(Qt, c(r, r,
            (n - 1) * tv[3] + 1)), a1 = array(a1, c(m)), P1 = array(P1,
            c(m, m)), at = at, Pt = Pt, vtuni = vtuni, Ftuni = Ftuni,
        Ktuni = Ktuni, Pinf = Pinf, Pstar = Pstar, Finfuni = Finfuni,
        Fstaruni = Fstaruni, Kinfuni = Kinfuni, Kstaruni = Kstaruni,
        d = d, j = j, p = p, m = m, r = r, n = n, lik = lik,
        optcal = optcal, info = info, vt = vt, Ft = Ft, Kt = Kt,
        Lt = Lt, Finf = Finf, Fstar = Fstar, Kinf = Kinf, Kstar = Kstar,
        Linf = Linf, Lstar = Lstar, tol = tol)
    kfout$tv[4] <- !(is.na(dim(as.array(Ht))[3]) || dim(as.array(Ht))[3] ==
        1)
    kfout$tv[5] <- !(is.na(dim(as.array(Zt))[3]) || dim(as.array(Zt))[3] ==
        1)
    kfout$Pinf <- array(kfout$Pinf[, , 1:(kfout$d + 1)], c(m,
        m, (kfout$d + 1) * (kfout$d > 0)))
    kfout$Pstar <- array(kfout$Pstar[, , 1:(kfout$d + 1)], c(m,
        m, (kfout$d + 1) * (kfout$d > 0)))
    kfout$Finfuni <- array(kfout$Finfuni[, 1:kfout$d], c(p, kfout$d))
    kfout$Fstaruni <- array(kfout$Fstaruni[, 1:kfout$d], c(p,
        kfout$d))
    kfout$Kinfuni <- array(kfout$Kinfuni[, , 1:kfout$d], c(m,
        p, kfout$d))
    kfout$Kstaruni <- array(kfout$Kstaruni[, , 1:kfout$d], c(m,
        p, kfout$d))
    kfout$yt <- yt
    kfout$Tt <- Tt
    kfout$Rt <- Rt
    kfout$Qt <- Qt
    kfout$Zt <- Zt
    kfout$Ht <- Ht
    Zt <- array(Zt, c(p, m, (n - 1) * kfout$tv[5] + 1))
    Ht <- array(Ht, c(p, p, (n - 1) * kfout$tv[4] + 1))
    if (kfout$d > 0) {
        kfout$Kt[, , 1:kfout$d] <- NA
        kfout$Ft[, , 1:kfout$d] <- NA
        kfout$Pt[, , 1:kfout$d] <- NA
        for (i in 1:kfout$d) {
            if (ydimt[i] != p) {
                kfout$Finfuni[(ydimt[i] + 1):p, i] <- kfout$Fstaruni[(ydimt[i] +
                  1):p, i] <- kfout$vtuni[(ydimt[i] + 1):p, i] <- NA
                if (sum(optcal) > 0) {
                  if (optcal[2] == 1) {
                    kfout$Fstar[, , i] <- matrix(Zt[, , (i -
                      1) * kfout$tv[5] + 1], p, m) %*% kfout$Pstar[,
                      , i] %*% t(matrix(Zt[, , (i - 1) * kfout$tv[5] +
                      1], p, m)) + Ht[, , (i - 1) * kfout$tv[4] +
                      1]
                    kfout$Finf[, , i] <- matrix(Zt[, , (i - 1) *
                      kfout$tv[5] + 1], p, m) %*% kfout$Pinf[,
                      , i] %*% t(matrix(Zt[, , (i - 1) * kfout$tv[5] +
                      1], p, m))
                  }
                  if (optcal[1]) {
                    kfout$vt[!ymiss[, i], i] <- kfout$vt[1:ydimt[i],
                      i]
                    kfout$vt[ymiss[, i], i] <- NA
                  }
                  if (optcal[3]) {
                    kfout$Kinf[, !is.na(yt[, i]), i] <- kfout$Kinf[,
                      1:kfout$ydimt[i], i]
                    kfout$Kinf[, is.na(yt[, i]), i] <- NA
                    kfout$Kstar[, !is.na(yt[, i]), i] <- kfout$Kstar[,
                      1:kfout$ydimt[i], i]
                    kfout$Kstar[, is.na(yt[, i]), i] <- NA
                  }
                }
            }
        }
    }
    if (kfout$d < n) {
        for (i in (kfout$d + 1):n) {
            if (ydimt[i] != p) {
                kfout$Ftuni[(ydimt[i] + 1):p, i] <- kfout$vtuni[(ydimt[i] +
                  1):p, i] <- NA
                if (sum(optcal) > 0) {
                  if (optcal[2] == 1) {
                    kfout$Ft[, , i] <- matrix(Zt[, , (i - 1) *
                      kfout$tv[5] + 1], p, m) %*% kfout$Pt[,
                      , i] %*% t(matrix(Zt[, , (i - 1) * kfout$tv[5] +
                      1], p, m)) + Ht[, , (i - 1) * kfout$tv[4] +
                      1]
                  }
                  if (optcal[1]) {
                    kfout$vt[!ymiss[, i], i] <- kfout$vt[1:ydimt[i],
                      i]
                    kfout$vt[ymiss[, i], i] <- NA
                  }
                  if (optcal[3]) {
                    kfout$Kt[, !is.na(yt[, i]), i] <- kfout$Kt[,
                      1:kfout$ydimt[i], i]
                    kfout$Kt[, is.na(yt[, i]), i] <- NA
                  }
                }
            }
        }
    }
    if (optcal[1] == 0)
        kfout$vt <- NULL
    if (optcal[2] == 0) {
        kfout$Ft <- kfout$Finf <- kfout$Fstar <- NULL
    }
    else {
        kfout$Finf <- array(kfout$Finf[, , 1:kfout$d], c(p, p,
            kfout$d))
        kfout$Fstar <- array(kfout$Fstar[, , 1:kfout$d], c(p,
            p, kfout$d))
    }
    if (optcal[3] == 0) {
        kfout$Kt <- kfout$Kinf <- kfout$Kstar <- NULL
    }
    else {
        kfout$Kinf <- array(kfout$Kinf[, , 1:kfout$d], c(m, p,
            kfout$d))
        kfout$Kstar <- array(kfout$Kstar[, , 1:kfout$d], c(m,
            p, kfout$d))
    }
    if (optcal[4] == 0) {
        kfout$Lt <- kfout$Linf <- kfout$Lstar <- NULL
    }
    else {
        kfout$Linf <- array(kfout$Linf[, , 1:kfout$d], c(m, m,
            kfout$d))
        kfout$Lstar <- array(kfout$Lstar[, , 1:kfout$d], c(m,
            m, kfout$d))
    }
    if (sum(kfout$info) != 0) {
        if (kfout$info[1]) {
            kfout$lik <- -Inf
            warning("Could not diagonalize Ht")
        }
        if (kfout$info[2]) {
            warning("Could not compute multivariate Kstar because multivariate Fstar is singular")
        }
        if (kfout$info[3]) {
            warning("Could not compute multivariate Kinf because multivariate Finf is singular")
        }
        if (kfout$info[4]) {
            warning("Could not compute multivariate Kt because multivariate Ft is singular")
        }
    }
    return(kfout)
}
