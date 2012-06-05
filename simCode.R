########################################################################
## This script contains functions for simulation study
########################################################################

library(fGarch)

## NOTE (Michael): Need to check on these formulas
AIC <- function(ll, k){
  -2 * (ll - k)
}

AICC <- function(ll, n, k){
  AIC(ll,  k) + (2 * k * (k + 1))/(n - k + 1)
}

BIC <- function(ll, n,  k){
  -2 * (ll - k * log(n))
}


## NOTE (Michael): The violation is assumed to be independent, and thus
##                 the likelihood ratio test can be bootstrapped.
##
## NOTE (Michael): Do the AIC outside the loop after obtaining all the
##                 likelihood, this will make the computation much more
##                 efficient. Will have to save a vector of number of
##                 parameters then.
##
## NOTE (Michael): Write a print function for the mgarch class, need to
##                 change class. (second priority)
##
## NOTE (Michael): Also return which model failed
##
## NOTE (Michael): The VaRtest fails if there are no violation.
##
## NOTE (Michael): Need to correct the VaRtest, current calculates the
##                 values, but will need to change the code as the true
##                 likelihood ratio statistics are actually in row 1, 3,
##                 5, and 2, 4, 6 for each respective alpha level.
##
modelSim <- function(Data, windowSize = 1000, alpha = c(0.05, 0.01),
                     n.ahead = 1, trace = FALSE, plot = FALSE){
  if(n.ahead != 1)
    warning("prediction for skewed scale mixture normal is incorrect")
  
  .StartTime = Sys.time()
  myinit <- NULL
  Data <- c(data.matrix(Data))
  T = length(Data)
  n.alpha <- length(alpha)
  if(windowSize > T - 1)
    stop("time series not sufficiently long enough for this window size")
  n.slide = T - windowSize
  ll <- matrix(0, nc = 7, nr = n.slide)
  colnames(ll) <- c("Normal", "Standard-t", "Ged", "S.Normal",
                    "S.Standard-t", "S.Ged", "S.Mixture Normal")
  n.param = AIC = AICC = BIC = pll = ll
  n.param[, 1] = 3
  n.param[, 2] = 4
  n.param[, 3] = 4
  n.param[, 4] = 4
  n.param[, 5] = 5
  n.param[, 6] = 5
  n.param[, 7] = NA
  
  norm.vio = matrix(0, nc = n.alpha, nr = n.slide)
  std.vio = matrix(0, nc = n.alpha, nr = n.slide)
  ged.vio = matrix(0, nc = n.alpha, nr = n.slide)
  snorm.vio = matrix(0, nc = n.alpha, nr = n.slide)
  sstd.vio = matrix(0, nc = n.alpha, nr = n.slide)
  sged.vio = matrix(0, nc = n.alpha, nr = n.slide)
  ssmnorm.vio = matrix(0, nc = n.alpha, nr = n.slide)

  norm.coef = matrix(NA, nc = 3, nr = n.slide)
  std.coef = matrix(NA, nc = 4, nr = n.slide)
  ged.coef = matrix(NA, nc = 4, nr = n.slide)
  snorm.coef = matrix(NA, nc = 4, nr = n.slide)
  sstd.coef = matrix(NA, nc = 5, nr = n.slide)
  sged.coef = matrix(NA, nc = 5, nr = n.slide)
  ssmnorm.coef = matrix(NA, nc = 5, nr = n.slide)
  ssmnorm.mix = vector("list", n.slide)


  ## Determine whether plot should be plotted.
  cnmPlot <- ifelse(plot, "gradient", "null")
  
  ## Start subset the time series and slide
  for(i in 1:n.slide){
    cat(paste("Window frame ", i, " of ", n.slide, "\n", sep = ""))
    subData <- Data[1:windowSize + i - 1]
    

    ## Fit the model
    norm.fit <- try(garchFit(data = subData, cond.dist = "norm",
      include.mean = FALSE, trace = trace))
    std.fit <- try(garchFit(data = subData, cond.dist = "std",
      include.mean = FALSE, trace = trace))
    ged.fit <- try(garchFit(data = subData, cond.dist = "ged",
      include.mean = FALSE, trace = trace))
    snorm.fit <- try(garchFit(data = subData, cond.dist = "snorm",
      include.mean = FALSE, trace = trace))
    sstd.fit <- try(garchFit(data = subData, cond.dist = "sstd",
      include.mean = FALSE, trace = trace))
    sged.fit <- try(garchFit(data = subData, cond.dist = "sged",
      include.mean = FALSE, trace = trace))
    ssmnorm.fit <- try(cnmms(as.mgarch(subData), grid = 3000,
      plot = cnmPlot, init = myinit))

    ## Update initial value
    if(!inherits(ssmnorm.fit, "try-error") &&
       valid.mgarch(x = Data[1:windowSize + i], beta = coef(ssmnorm.fit),
                    mix = ssmnorm.fit$mix)){
      myinit = list(beta = coef(ssmnorm.fit), mix = ssmnorm.fit$mix)
    } else {
      myinit <- NULL
    }
    
    ## Save the likelihood and the information criteria
    if(!inherits(norm.fit, "try-error"))
      ll[i, 1] = -norm.fit@fit$llh
    if(!inherits(std.fit, "try-error"))
      ll[i, 2] = -std.fit@fit$llh
    if(!inherits(ged.fit, "try-error"))
      ll[i, 3] = -ged.fit@fit$llh
    if(!inherits(snorm.fit, "try-error"))
      ll[i, 4] = -snorm.fit@fit$llh
    if(!inherits(sstd.fit, "try-error"))
      ll[i, 5] = -sstd.fit@fit$llh
    if(!inherits(sged.fit, "try-error"))
      ll[i, 6] = -sged.fit@fit$llh
    if(!inherits(ssmnorm.fit, "try-error"))
      ll[i, 7] = ssmnorm.fit$ll

    ## Save the coefficients
    if(!inherits(norm.fit, "try-error"))
      norm.coef[i, ] = coef(norm.fit)
    if(!inherits(std.fit, "try-error"))
      std.coef[i, ] = coef(std.fit)
    if(!inherits(ged.fit, "try-error"))
      ged.coef[i, ] = coef(ged.fit)
    if(!inherits(snorm.fit, "try-error"))
      snorm.coef[i, ] = coef(snorm.fit)
    if(!inherits(sstd.fit, "try-error"))
      sstd.coef[i, ] = coef(sstd.fit)
    if(!inherits(sged.fit, "try-error"))
      sged.coef[i, ] = coef(sged.fit)
    if(!inherits(ssmnorm.fit, "try-error"))
      ssmnorm.coef[i, ] = coef(ssmnorm.fit)
    if(!inherits(ssmnorm.fit, "try-error")){
      ssmnorm.mix[[i]] = ssmnorm.fit$mix
      n.param[i, 7] = 2 * NROW(ssmnorm.mix[[i]]) + 1
    }

    ## Calculate the VaR bound
    if(!inherits(norm.fit, "try-error")){
      norm.psd <- predict(norm.fit, n.ahead = n.ahead)$standardDeviation
      norm.VaR <-  norm.psd * qnorm(alpha)
    }
    if(!inherits(std.fit, "try-error")){
      std.psd <- predict(std.fit, n.ahead = n.ahead)$standardDeviation
      std.VaR <- std.psd * qstd(alpha, nu = coef(std.fit)["shape"])
    }
    if(!inherits(ged.fit, "try-error")){
      ged.psd <- predict(ged.fit, n.ahead = n.ahead)$standardDeviation
      ged.VaR <- ged.psd * qged(alpha, nu = coef(ged.fit)["shape"])
    }
    if(!inherits(snorm.fit, "try-error")){
      snorm.psd <- predict(snorm.fit, n.ahead = n.ahead)$standardDeviation
      snorm.VaR <-  snorm.psd * qsnorm(alpha, xi = coef(snorm.fit)["skew"])
    }
    if(!inherits(sstd.fit, "try-error")){
      sstd.psd <- predict(sstd.fit, n.ahead = n.ahead)$standardDeviation
      sstd.VaR <- sstd.psd * qsstd(alpha, xi = coef(sstd.fit)["skew"],
                                   nu = coef(sstd.fit)["shape"])
    }
    if(!inherits(sged.fit, "try-error")){
      sged.psd <- predict(sged.fit, n.ahead = n.ahead)$standardDeviation
      sged.VaR <-  sged.psd * qsged(alpha, xi = coef(sged.fit)["skew"],
                                   nu = coef(sged.fit)["shape"])
    }
    if(!inherits(ssmnorm.fit, "try-error")){
      mybeta <- coef(ssmnorm.fit)
      ssmnorm.psd <- sqrt(mybeta[1] + mybeta[3] * subData[windowSize]^2 +
                          mybeta[2]  * ssmnorm.fit$sigma.t[windowSize]^2)
      q <- try(qssmnorm(alpha, varmix = ssmnorm.fit$mix,
                        xi = coef(ssmnorm.fit)["xi"]))
      print(ssmnorm.fit$mix)
      print(coef(ssmnorm.fit)["xi"])
      ssmnorm.VaR <- ssmnorm.psd * q
    }
      
    ## Calculate the likelihood of the predicted value
    if(!inherits(norm.fit, "try-error"))
      pll[i, 1] = dnorm(Data[windowSize + n.ahead + i]/norm.psd)
    if(!inherits(std.fit, "try-error"))
      pll[i, 2] = dstd(Data[windowSize + n.ahead + i]/std.psd,
           nu = coef(std.fit)["shape"])
    if(!inherits(ged.fit, "try-error"))
      pll[i, 3] = dged(Data[windowSize + n.ahead + i]/ged.psd,
           nu = coef(std.fit)["shape"])
    if(!inherits(snorm.fit, "try-error"))
      pll[i, 4] = dsnorm(Data[windowSize + n.ahead + i]/snorm.psd,
           xi = coef(snorm.fit)["skew"])
    if(!inherits(sstd.fit, "try-error"))
      pll[i, 5] = dsstd(Data[windowSize + n.ahead + i]/sstd.psd,
           nu = coef(sstd.fit)["shape"], coef(sstd.fit)["skew"])
    if(!inherits(sged.fit, "try-error"))
      pll[i, 6] = dsged(Data[windowSize + n.ahead + i]/sged.psd,
           nu = coef(sged.fit)["shape"], coef(sged.fit)["skew"])
    if(!inherits(ssmnorm.fit, "try-error"))
      pll[i, 7] = dssmnorm(Data[windowSize + n.ahead + i]/ssmnorm.psd,
           varmix = ssmnorm.fit$mix, xi = coef(ssmnorm.fit)["xi"])

    ## Test if violates
    if(!inherits(norm.fit, "try-error"))
      norm.vio[i, ] <- Data[windowSize + n.ahead + i] <= norm.VaR
    if(!inherits(std.fit, "try-error"))
      std.vio[i, ] <- Data[windowSize + n.ahead + i] <= std.VaR
    if(!inherits(ged.fit, "try-error"))
      ged.vio[i, ] <- Data[windowSize + n.ahead + i] <= ged.VaR
    if(!inherits(snorm.fit, "try-error"))
      snorm.vio[i, ] <- Data[windowSize + n.ahead + i] <= snorm.VaR
    if(!inherits(sstd.fit, "try-error"))
      sstd.vio[i, ] <- Data[windowSize + n.ahead + i] <= sstd.VaR
    if(!inherits(sged.fit, "try-error"))
      sged.vio[i, ] <- Data[windowSize + n.ahead + i] <= sged.VaR
    if(!inherits(ssmnorm.fit, "try-error"))
      ssmnorm.vio[i, ] <- Data[windowSize + n.ahead + i] <= ssmnorm.VaR

    plotgrad(as.mgarch(subData), ssmnorm.coef[i, ], ssmnorm.mix[[i]])
  }


  AIC = -2 * (ll - n.param)
  AICC = AIC * (n.param * (n.param + 1)/(T - n.param + 1))
  BIC = -2 * ll + n.param * log(T) 
  
  
  ## Temporary function to recover model failure by replacing NA with
  ## binomial random variables.
  na.sub <- function(x, alpha){
    x[is.na(x)] <- rbinom(sum(is.na(x)), 1,
                          rep(alpha, each = sum(is.na(x))/length(alpha)))
    x
  }

  ## Think of a way of accounting for NA and model failures, also fix
  ## why model fails
  if(sum(is.na(norm.vio)) < length(norm.vio) * 0.25){
    norm.vio <- na.sub(norm.vio, alpha = alpha)
    norm.test <- try(apply(norm.vio, 2, VaRtest, alpha = alpha))
  }
  if(sum(is.na(std.vio)) < length(std.vio) * 0.25){
    std.vio <- na.sub(std.vio, alpha = alpha)
    std.test <- try(apply(std.vio, 2, VaRtest, alpha = alpha))
  }
  if(sum(is.na(ged.vio)) < length(ged.vio) * 0.25){
    ged.vio <- na.sub(ged.vio, alpha = alpha)
    ged.test <- try(apply(ged.vio, 2, VaRtest, alpha = alpha))
  }
  if(sum(is.na(snorm.vio)) < length(snorm.vio) * 0.25){
    snorm.vio <- na.sub(snorm.vio, alpha = alpha)
    snorm.test <- try(apply(snorm.vio, 2, VaRtest, alpha = alpha))
  }
  if(sum(is.na(sstd.vio)) < length(sstd.vio) * 0.25){
    sstd.vio <- na.sub(sstd.vio, alpha = alpha)
    sstd.test <- try(apply(sstd.vio, 2, VaRtest, alpha = alpha))
  }
  if(sum(is.na(sged.vio)) < length(sged.vio) * 0.25){
    sged.vio <- na.sub(sged.vio, alpha = alpha)
    sged.test <- try(apply(sged.vio, 2, VaRtest, alpha = alpha))
  }
  if(sum(is.na(ssmnorm.vio)) < length(ssmnorm.vio) * 0.25){
    ssmnorm.vio <- na.sub(ssmnorm.vio, alpha = alpha)
    ssmnorm.test <- try(apply(ssmnorm.vio, 2, VaRtest, alpha = alpha))
  }
  
  .EndTime = Sys.time()
  .timeDiff = .EndTime - .StartTime
  
  norm.fit <- list(norm.vio = norm.vio,
                   norm.test = norm.test,
                   norm.coef = norm.coef)
  std.fit <- list(std.vio = std.vio,
                  std.test = std.test,
                  std.coef = std.coef)
  ged.fit <- list(ged.vio = ged.vio,
                  ged.test = ged.test,
                  ged.coef = ged.coef)
  snorm.fit <- list(snorm.vio = snorm.vio,
                    snorm.test = snorm.test,
                    snorm.coef = snorm.coef)
  sstd.fit <- list(sstd.vio = sstd.vio,
                   sstd.test = sstd.test,
                   sstd.coef = sstd.coef)
  sged.fit <- list(sged.vio = sged.vio,
                   sged.test = sged.test,
                   sged.coef = sged.coef)
  ssmnorm.fit <- list(ssmnorm.vio = ssmnorm.vio,
                      ssmnorm.test = ssmnorm.test,
                      ssmnorm.coef = ssmnorm.coef)

  list(ll = ll, pll = pll, AIC = AIC, AICC = AICC, BIC = BIC,
       norm.fit = norm.fit, std.fit = std.fit, ged.fit = ged.fit,
       snorm.fit = snorm.fit, sstd.fit = sstd.fit, sged.fit = sged.fit,
       ssmnorm.fit = ssmnorm.fit, ssmnorm.mix = ssmnorm.mix,
       StartTime = .StartTime, EndTime = .EndTime, TotalTime = .timeDiff)
}





VaRtest <- function(ht, alpha){
  T = length(ht)
  sht = sum(ht)
  LRuc = 2 * (dbinom(sht, T, sht/T, log = TRUE) -
    dbinom(sht, T, alpha, log = TRUE))
  viovec = cbind(ht[-T], ht[-1])
  pivec = integer(NROW(viovec))
  for(i in 1:NROW(viovec)){
    if(viovec[i, 1] == 0 & viovec[i, 2] == 0){
      pivec[i] = 1
    } else if(viovec[i, 1] == 0 & viovec[i, 2] == 1){
      pivec[i] = 2
    } else if(viovec[i, 1] == 1 & viovec[i, 2] == 0){
      pivec[i] = 3
    } else {
      pivec[i] = 4
    }
  }
  tab = table(pivec)
  n00 = tab[1]
  n01 = tab[2]
  n10 = tab[3]
  n11 = ifelse(dim(tab) == 4, tab[4], 0)
  pi01 = n01/(n00 + n01)
  pi11 = n11/(n10 + n11)
  LRind = 2 * (log((1 - pi01)^n00 * pi01^n01 * (1 - pi11)^n10 *
    pi11^n11) - log((1 - alpha)^(n00 + n10) * alpha^(n01 + n11)))
  LRucStats = dchisq(LRuc, df = 1)
  LRindStats = dchisq(LRind, df = 1)
  LRcc = LRuc + LRind
  LRccStats = dchisq(LRcc, df = 2)
  c(LRucStats = LRucStats, LRindStats = LRindStats, LRccStats = LRccStats)
}



## NOTE (Michael): This is the old simulation code, which has no sliding.
## modelSim <- function(inSampData, outSampData, n = length(inSampData),
##                      n.sim = 10, cl = c(0.05, 0.01)){
##   .StartTime = Sys.time()
##   inSampData <- c(data.matrix(inSampData))
##   outSampData <- c(data.matrix(outSampData))
##   T1 = length(inSampData)
##   T2 = length(outSampData)
##   ks <- matrix(NA, nc = 7, nr = n.sim)
##   VaR <- matrix(NA, nc = 4, nr = length(cl) * 14)
##   rownames(VaR) = paste(rep(rep(c("in ", "out "), each = 7), length(cl)),
##             c("Normal", "Standard-t", "Ged", "S.Normal",
##             "S.Standard-t", "S.Ged", "S.Mixture Normal"), "(",
##             rep(cl, each = 14), ")", sep = "")
##   colnames(VaR) = c("Violation", "LRuc", "LRind", "LRcc")
##   ics <- matrix(NA, nc = 4, nr = 7)
##   colnames(ics) = c("Log-likelihood", "AIC", "AICC", "BIC")
##   rownames(ics) = c("Normal", "Standard-t", "Ged", "S.Normal", "S.Standard-t",
##          "S.Ged", "S.Mixture Normal")
##   ## Fit the model
##   norm.fit = try(garchFit(data = inSampData, cond.dist = "norm",
##     include.mean = FALSE))
##   std.fit = try(garchFit(data = inSampData, cond.dist = "std",
##     include.mean = FALSE))
##   ged.fit = try(garchFit(data = inSampData, cond.dist = "ged",
##     include.mean = FALSE))
##   snorm.fit = try(garchFit(data = inSampData, cond.dist = "snorm",
##     include.mean = FALSE))
##   sstd.fit = try(garchFit(data = inSampData, cond.dist = "sstd",
##     include.mean = FALSE))
##   sged.fit = try(garchFit(data = inSampData, cond.dist = "sged",
##     include.mean = FALSE))
##   ssmnorm.fit = try(cnmms(as.mgarch(inSampData), grid = 3000, plot = "gradient"))
##   ## Calculate the information criteria, need to check BIC and maybe
##   ## remove AICC since it is equivalent to AIC in large samples
##   if(!inherits(norm.fit, "try-error")){
##     ll = -norm.fit@fit$ll
##     ics[1, ] = c(ll, AIC(ll, k = 3), AICC(ll, n = T1, k = 3),
##          BIC(ll, n = T1, k = 3))
##   }
##   if(!inherits(std.fit, "try-error")){
##     ll = -std.fit@fit$ll
##     ics[2, ] = c(ll, AIC(ll, k = 4), AICC(ll, n = T1, k = 4),
##          BIC(ll, n = T1, k = 4))
##   }
##   if(!inherits(ged.fit, "try-error")){
##     ll = -ged.fit@fit$ll
##     ics[3, ] = c(ll, AIC(ll, k = 4), AICC(ll, n = T1, k = 4),
##          BIC(ll, n = T1, k = 4))
##   }
##   if(!inherits(snorm.fit, "try-error")){
##     ll = -snorm.fit@fit$ll
##     ics[4, ] = c(ll, AIC(ll, k = 4), AICC(ll, n = T1, k = 4),
##          BIC(ll, n = T1, k = 4))
##   }
##   if(!inherits(sstd.fit, "try-error")){
##     ll = -sstd.fit@fit$ll
##     ics[5, ] = c(ll, AIC(ll, k = 5), AICC(ll, n = T1, k = 5),
##          BIC(ll, n = T1, k = 5))
##   }
##   if(!inherits(sged.fit, "try-error")){
##     ll = -sged.fit@fit$ll
##     ics[6, ] = c(ll, AIC(ll, k = 5), AICC(ll, n = T1, k = 5),
##          BIC(ll, n = T1, k = 5))
##   }
##   if(!inherits(ssmnorm.fit, "try-error")){
##     ll = ssmnorm.fit$ll
##     if(length(ssmnorm.fit$mix$pt) == 1){
##       n.param = length(ssmnorm.fit$beta)
##     } else {
##       n.param = length(ssmnorm.fit$beta) + 2 * (length(ssmnorm.fit$mix$pt) - 1)
##     }
##     ics[7, ] = c(ll, AIC(ll, k = n.param), AICC(ll, n = T1, k = n.param),
##          BIC(ll, n = T1, k = n.param))
##   }
##   ## Don't run the simulation if the scale normal mixture fails
##   if(inherits(ssmnorm.fit, "try-error"))
##     stop("mixture normal fit failed")
##   ## Print the coefficient and the mixture, might want to save this
##   if(!inherits(norm.fit, "try-error"))
##     print(coef(norm.fit))
##   if(!inherits(std.fit, "try-error"))
##     print(coef(std.fit))
##   if(!inherits(ged.fit, "try-error"))
##     print(coef(ged.fit))
##   if(!inherits(snorm.fit, "try-error"))
##     print(coef(snorm.fit))
##   if(!inherits(sstd.fit, "try-error"))
##     print(coef(sstd.fit))
##   if(!inherits(sged.fit, "try-error"))
##     print(coef(sged.fit))
##   if(!inherits(ssmnorm.fit, "try-error"))
##     print(ssmnorm.fit$beta)
##   if(!inherits(ssmnorm.fit, "try-error"))
##     print(ssmnorm.fit$mix)
##   ## Generate the spec of each model for simulation
##   if(!inherits(norm.fit, "try-error"))
##     norm.spec = garchSpec(model = as.list(norm.fit@fit$par),
##       cond.dist = "norm")
##   if(!inherits(std.fit, "try-error"))
##     std.spec = garchSpec(model = as.list(std.fit@fit$par),
##       cond.dist = "std")
##   if(!inherits(ged.fit, "try-error"))
##     ged.spec = garchSpec(model = as.list(ged.fit@fit$par),
##       cond.dist = "ged")
##   if(!inherits(snorm.fit, "try-error"))
##     snorm.spec = garchSpec(model = as.list(snorm.fit@fit$par),
##       cond.dist = "snorm")
##   if(!inherits(sstd.fit, "try-error"))
##     sstd.spec = garchSpec(model = as.list(sstd.fit@fit$par),
##       cond.dist = "sstd")
##   if(!inherits(sged.fit, "try-error"))
##     sged.spec = garchSpec(model = as.list(sged.fit@fit$par),
##       cond.dist = "sged")
##   if(!inherits(ssmnorm.fit, "try-error"))
##     ssmnorm.spec = garchSpec(model = as.list(ssmnorm.fit$beta),
##       cond.dist = "ssmnorm")
##   dist = c("norm", "std", "ged", "snorm", "sstd", "sged", "ssmnorm")
##   color = c("red", "orange", "yellow", "green", "blue", "violet", "purple")
##   pass = c(!inherits(norm.fit, "try-error"),
##     !inherits(std.fit, "try-error"), !inherits(ged.fit, "try-error"),
##     !inherits(snorm.fit, "try-error"), !inherits(sstd.fit, "try-error"),
##     !inherits(sged.fit, "try-error"), !inherits(ssmnorm.fit, "try-error"))
##   ## Plot the conditional distributions
##   if(!inherits(norm.fit, "try-error"))
##     curve(dnorm(x), from = -5, to = 5, ylim = c(0, 0.5), col = "red",
##           main = "Conditional distribution")
##   if(!inherits(std.fit, "try-error"))
##     curve(dstd(x, nu = coef(std.fit)["shape"]), col = "orange", add = TRUE)
##   if(!inherits(ged.fit, "try-error"))
##     curve(dged(x, nu = coef(ged.fit)["shape"]), col = "yellow", add = TRUE)
##   if(!inherits(snorm.fit, "try-error"))
##     curve(dsnorm(x, xi = coef(snorm.fit)["skew"]), col = "green", add = TRUE)
##   if(!inherits(sstd.fit, "try-error"))
##     curve(dsstd(x, xi = coef(sstd.fit)["skew"],
##                 nu = coef(sstd.fit)["shape"]), col = "blue", add = TRUE)
##   if(!inherits(sged.fit, "try-error"))
##     curve(dsged(x, xi = coef(sged.fit)["skew"],
##                 nu = coef(sged.fit)["shape"]), col = "violet", add = TRUE)
##   if(!inherits(ssmnorm.fit, "try-error"))
##     curve(dmsnorm(x, varmix = ssmnorm.fit$mix, xi = ssmnorm.fit$beta["xi"]),
##           col = "purple", add = TRUE)
##   legend("topleft", legend = c(dist[pass]), lty = 1, lwd = 1,
##          col = color[pass], bty = "n")
##   for(j in 1:length(cl)){
##     ## In sample VaR
##     plot(inSampData, type = "l")
##     alpha = c(cl[j], 1 - cl[j])
##     if(!inherits(norm.fit, "try-error")){
##       k = qnorm(alpha)
##       VaR[1 + 14 * (j - 1), ] <- VaRtest(data = inSampData,
##                                          VaR = norm.fit@sigma.t * k[1],
##                                          alpha = alpha[1])
##       lines(norm.fit@sigma.t * k[1], col = "red")
##       lines(norm.fit@sigma.t * k[2], col = "red")
##     }
##     if(!inherits(std.fit, "try-error")){
##       k = qstd(alpha, nu = coef(std.fit)["shape"])
##       VaR[2 + 14 * (j - 1), ] <- VaRtest(data = inSampData,
##                                          VaR = std.fit@sigma.t * k[1],
##                                          alpha = alpha[1])
##       lines(std.fit@sigma.t * k[1], col = "orange")
##       lines(std.fit@sigma.t * k[2], col = "orange")
##     }
##     if(!inherits(ged.fit, "try-error")){
##       k = qged(alpha, nu = coef(ged.fit)["shape"])
##       VaR[3 + 14 * (j - 1), ] <- VaRtest(data = inSampData,
##                                          VaR = ged.fit@sigma.t * k[1],
##                                          alpha = alpha[1])
##       lines(ged.fit@sigma.t * k[1], col = "yellow")
##       lines(ged.fit@sigma.t * k[2], col = "yellow")
##     }
##     if(!inherits(snorm.fit, "try-error")){
##       k = qsnorm(alpha, xi = coef(snorm.fit)["skew"])
##       VaR[4 + 14 * (j - 1), ] <- VaRtest(data = inSampData,
##                                          VaR = snorm.fit@sigma.t * k[1],
##                                          alpha = alpha[1])
##       lines(snorm.fit@sigma.t * k[1], col = "green")
##       lines(snorm.fit@sigma.t * k[2], col = "green")
##     }
##     if(!inherits(sstd.fit, "try-error")){
##       k = qsstd(alpha, xi = coef(sstd.fit)["skew"],
##         nu = coef(sstd.fit)["shape"])
##       VaR[5 + 14 * (j - 1), ] <- VaRtest(data = inSampData,
##                                          VaR = sstd.fit@sigma.t * k[1],
##                                          alpha = alpha[1])
##       lines(sstd.fit@sigma.t * k[1], col = "blue")
##       lines(sstd.fit@sigma.t * k[2], col = "blue")
##     }
##     if(!inherits(sged.fit, "try-error")){
##       k = qsged(alpha, xi = coef(sged.fit)["skew"],
##         nu = coef(sged.fit)["shape"])
##       VaR[6 + 14 * (j - 1), ] <- VaRtest(data = inSampData,
##                                          VaR = sged.fit@sigma.t * k[1],
##                                          alpha = alpha[1])
##       lines(sged.fit@sigma.t * k[1], col = "violet")
##       lines(sged.fit@sigma.t * k[2], col = "violet")
##     }
##     if(!inherits(ssmnorm.fit, "try-error")){
##       k = qmsnorm(alpha, varmix = ssmnorm.fit$mix,
##         xi = ssmnorm.fit$beta["xi"])
##       VaR[7 + 14 * (j - 1), ] <- VaRtest(data = inSampData,
##                           VaR = cond.sd(inSampData, ssmnorm.fit$beta) * k[1],
##                           alpha = alpha[1])
##       lines(cond.sd(inSampData, ssmnorm.fit$beta) * k[1], col = "purple")
##       lines(cond.sd(inSampData, ssmnorm.fit$beta) * k[2], col = "purple")
##     }
##     ## Out of sample VaR
##     plot(outSampData, type = "l")
##     if(!inherits(norm.fit, "try-error")){
##       k = qnorm(alpha)
##       cs = cond.sd(c(inSampData, outSampData),
##         c(coef(norm.fit)[c(1, 3, 2)], norm.fit@sigma.t[1]))[(T1 + 1):(T1 + T2)]
##       VaR[8 + 14 * (j - 1), ] <- VaRtest(data = outSampData, VaR = cs * k[1],
##                           alpha = alpha[1])
##       lines(cs * k[1], col = "red")
##       lines(cs * k[2], col = "red")
##     }
##     if(!inherits(std.fit, "try-error")){
##       k = qstd(alpha, nu = coef(std.fit)["shape"])
##       cs = cond.sd(c(inSampData, outSampData),
##         c(coef(std.fit)[c(1, 3, 2)], std.fit@sigma.t[1]))[(T1 + 1):(T1 + T2)]
##       VaR[9 + 14 * (j - 1), ] <- VaRtest(data = outSampData, VaR = cs * k[1],
##                                          alpha = alpha[1])
##       lines(cs * k[1], col = "orange")
##       lines(cs * k[2], col = "orange")
##     }
##     if(!inherits(ged.fit, "try-error")){
##       k = qged(alpha, nu = coef(ged.fit)["shape"])
##       cs = cond.sd(c(inSampData, outSampData),
##         c(coef(ged.fit)[c(1, 3, 2)], ged.fit@sigma.t[1]))[(T1 + 1):(T1 + T2)]
##       VaR[10 + 14 * (j - 1), ] <- VaRtest(data = outSampData, VaR = cs * k[1],
##                                           alpha = alpha[1])
##       lines(cs * k[1], col = "yellow")
##       lines(cs * k[2], col = "yellow")
##     }
##     if(!inherits(snorm.fit, "try-error")){
##       k = qsnorm(alpha, xi = coef(snorm.fit)["skew"])
##       cs = cond.sd(c(inSampData, outSampData),
##         c(coef(snorm.fit)[c(1, 3, 2)], snorm.fit@sigma.t[1]))[(T1 + 1):(T1 + T2)]
##       VaR[11 + 14 * (j - 1), ] <- VaRtest(data = outSampData, VaR = cs * k[1],
##                                           alpha = alpha[1])
##       lines(cs * k[1], col = "green")
##       lines(cs * k[2], col = "green")
##     }
##     if(!inherits(sstd.fit, "try-error")){
##       k = qsstd(alpha, xi = coef(sstd.fit)["skew"],
##         nu = coef(sstd.fit)["shape"])
##       cs = cond.sd(c(inSampData, outSampData),
##         c(coef(sstd.fit)[c(1, 3, 2)], sstd.fit@sigma.t[1]))[(T1 + 1):(T1 + T2)]
##       VaR[12 + 14 * (j - 1), ] <- VaRtest(data = outSampData, VaR = cs * k[1],
##                                           alpha = alpha[1])
##       lines(cs * k[1], col = "blue")
##       lines(cs * k[2], col = "blue")
##     }
##     if(!inherits(sged.fit, "try-error")){
##       k = qsged(alpha, xi = coef(sged.fit)["skew"],
##         nu = coef(sged.fit)["shape"])
##       cs = cond.sd(c(inSampData, outSampData),
##         c(coef(sged.fit)[c(1, 3, 2)], sged.fit@sigma.t[1]))[(T1 + 1):(T1 + T2)]
##       VaR[13 + 14 * (j - 1), ] <- VaRtest(data = outSampData, VaR = cs * k[1],
##                                           alpha = alpha[1])
##       lines(cs * k[1], col = "violet")
##       lines(cs * k[2], col = "violet")
##     }
##     if(!inherits(ssmnorm.fit, "try-error")){
##       k = qmsnorm(alpha, varmix = ssmnorm.fit$mix,
##         xi = ssmnorm.fit$beta["xi"])
##       cs = cond.sd(c(inSampData, outSampData),
##         ssmnorm.fit$beta)[(T1 + 1):(T1 + T2)]
##       VaR[14 + 14 * (j - 1), ] <- VaRtest(data = outSampData, VaR = cs * k[1],
##                            alpha = alpha[1])
##       lines(cs * k[1], col = "purple")
##       lines(cs * k[2], col = "purple")
##     }
##   }
##   for(i in 1:n.sim){
##     print(paste("Simulation: ", i, sep = ""))
##     ## Simulate the time series from each model
##     if(!inherits(norm.fit, "try-error"))
##       norm.sim = garchSim(n = n, spec = norm.spec)
##     if(!inherits(std.fit, "try-error"))
##       std.sim = garchSim(n = n, spec = std.spec)
##     if(!inherits(ged.fit, "try-error"))
##       ged.sim = garchSim(n = n, spec = ged.spec)
##     if(!inherits(snorm.fit, "try-error"))
##       snorm.sim = garchSim(n = n, spec = snorm.spec)
##     if(!inherits(sstd.fit, "try-error"))
##       sstd.sim = garchSim(n = n, spec = sstd.spec)
##     if(!inherits(sged.fit, "try-error"))
##       sged.sim = garchSim(n = n, spec = sged.spec)
##     if(!inherits(ssmnorm.fit, "try-error"))
##       ssmnorm.sim = garchSim(n = n, spec = ssmnorm.spec)
##     hist(inSampData, freq = FALSE, breaks = n/10,
##          main = paste("Unconditional distribution of simulation ", i, sep = ""),
##          ylim = c(0, 2 * max(hist(inSampData, plot = FALSE,
##            breaks = n/10)$density)))
##     if(!inherits(norm.fit, "try-error"))
##       lines(density(norm.sim, kernel = "epanechnikov"), col = "red")
##     if(!inherits(std.fit, "try-error"))
##       lines(density(std.sim, kernel = "epanechnikov"), col = "orange")
##     if(!inherits(ged.fit, "try-error"))
##       lines(density(ged.sim, kernel = "epanechnikov"), col = "yellow")
##     if(!inherits(snorm.fit, "try-error"))
##       lines(density(snorm.sim, kernel = "epanechnikov"), col = "green")
##     if(!inherits(sstd.fit, "try-error"))
##       lines(density(sstd.sim, kernel = "epanechnikov"), col = "blue")
##     if(!inherits(sged.fit, "try-error"))
##       lines(density(sged.sim, kernel = "epanechnikov"), col = "violet")
##     if(!inherits(ssmnorm.fit, "try-error"))
##       lines(density(ssmnorm.sim, kernel = "epanechnikov"), col = "purple")
##     legend("topleft", legend = c(dist[pass]), lty = 1, lwd = 1,
##            col = color[pass], bty = "n")
##     ## Calculate the KS-statistic between the two sample (Twosided),
##     ## maybe we can implement just one side.
##     if(!inherits(norm.fit, "try-error"))
##       ks[, 1] = ks2Test(norm.sim, inSampData)@test$statistic[1]
##     if(!inherits(std.fit, "try-error"))
##       ks[, 2] = ks2Test(std.sim, inSampData)@test$statistic[1]
##     if(!inherits(ged.fit, "try-error"))
##       ks[, 3] = ks2Test(ged.sim, inSampData)@test$statistic[1]
##     if(!inherits(snorm.fit, "try-error"))
##       ks[, 4] = ks2Test(snorm.sim, inSampData)@test$statistic[1]
##     if(!inherits(sstd.fit, "try-error"))
##       ks[, 5] = ks2Test(sstd.sim, inSampData)@test$statistic[1]
##     if(!inherits(sged.fit, "try-error"))
##       ks[, 6] = ks2Test(sged.sim, inSampData)@test$statistic[1]
##     if(!inherits(ssmnorm.fit, "try-error"))
##       ks[, 7] = ks2Test(ssmnorm.sim, inSampData)@test$statistic[1]
##   }
##   meanKs = colMeans(ks, na.rm = TRUE)
##   names(meanKs) = c("Normal", "Standard-t", "Ged", "S.Norm", "S.Standard-t",
##             "S.Ged", "S.Mixture Normal")
##   .EndTime = Sys.time()
##   timediff = .EndTime - .StartTime
##   list(meanKs = meanKs, VaR = VaR, ics = ics, StartTime = .StartTime,
##        EndTime = .EndTime, TotalTime = timediff)
## }

## NOTE (Michael): This works with the old modelSim function
## VaRtest <- function(data, VaR, alpha){
##   data = c(data.matrix(data))
##   T = length(data)
##   ht = data <= VaR
##   sht = sum(ht)
##   LRuc = 2 * (dbinom(sht, T, sht/T, log = TRUE) -
##     dbinom(sht, T, alpha, log = TRUE))
##   viovec = cbind(ht[-T], ht[-1])
##   pivec = integer(NROW(viovec))
##   for(i in 1:NROW(viovec)){
##     if(viovec[i, 1] == 0 & viovec[i, 2] == 0){
##       pivec[i] = 1
##     } else if(viovec[i, 1] == 0 & viovec[i, 2] == 1){
##       pivec[i] = 2
##     } else if(viovec[i, 1] == 1 & viovec[i, 2] == 0){
##       pivec[i] = 3
##     } else {
##       pivec[i] = 4
##     }
##   }
##   tab = table(pivec)
##   n00 = tab[1]
##   n01 = tab[2]
##   n10 = tab[3]
##   n11 = ifelse(dim(tab) == 4, tab[4], 0)
##   pi01 = n01/(n00 + n01)
##   pi11 = n11/(n10 + n11)
##   LRind = 2 * (log((1 - pi01)^n00 * pi01^n01 * (1 - pi11)^n10 *
##     pi11^n11) - log((1 - alpha)^(n00 + n10) * alpha^(n01 + n11)))
##   LRucStats = dchisq(LRuc, df = 1)
##   LRindStats = dchisq(LRind, df = 1)
##   LRcc = LRuc + LRind
##   LRccStats = dchisq(LRcc, df = 2)
##   c(Violation = sht, LRucStats = LRucStats, LRindStats = LRindStats,
##     LRccStats = LRccStats)
## }


## ## NOTE (Michael): Just work for GARCH(1, 1) for now
## edSim <- function(n.samp = 100, spec){
##   norm.rlt <- matrix(NA, nc = 6, nr = 1)
##   colnames(norm.rlt) <- c("omega", "alpha1", "beta1", "hd", "kld", "mise")
##   std.rlt <- matrix(NA, nc = 6, nr = 1)
##   colnames(std.rlt) <- c("omega", "alpha1", "beta1", "hd", "kld", "mise")
##   ged.rlt <- matrix(NA, nc = 6, nr = 1)
##   colnames(ged.rlt) <- c("omega", "alpha1", "beta1", "hd", "kld", "mise")
##   mnorm.rlt <- matrix(NA, nc = 6, nr = 1)
##   colnames(mnorm.rlt) <- c("omega", "alpha1", "beta1", "hd", "kld", "mise")
##   sim <- garchSim(spec = spec, n = n.samp)
##   model = spec@model
##   if(spec@distribution == "norm")
##     trueDist <- function(x) dnorm(x)
##   if(spec@distribution == "std")
##     trueDist <- function(x) dstd(x, nu = model$shape)
##   if(spec@distribution == "ged")
##     trueDist <- function(x) dged(x, nu = model$shape)
##   if(spec@distribution == "cauchy")
##     trueDist <- function(x) dcauchy(x, location = 0, scale = 1)
##   if(spec@distribution == "mnorm")
##     trueDist <- function(x) dmnorm(x, mean = 0, varmix = model$shape)
##   if(spec@distribution == "snorm")
##     trueDist <- function(x) dsnorm(x, xi = model$skew)
##   if(spec@distribution == "sstd")
##     trueDist <- function(x) dsstd(x, nu = model$shape, xi = model$skew)
##   if(spec@distribution == "sged")
##     trueDist <- function(x) dsged(x, nu = model$shape, xi = model$skew)
##   if(spec@distribution == "scauchy")
##     trueDist <- function(x) dscauchy(x, xi = model$skew)
##   if(spec@distribution == "msnorm")
##     trueDist <- function(x) dmsnorm(x, varmix = model$shape, xi = model$skew)
##   curve(trueDist(x), -5, 5)
##   ## Normal distribution
##   ## -------------------
##   norm.fit <- try(garchFit(data = sim, cond.dist = "snorm", trace = FALSE,
##                            include.mean = FALSE))
##   if(!inherits(norm.fit, "try-error")){
##     fitDist <- function(x) dsnorm(x, mean = 0, sd = 1,
##                                   xi = coef(norm.fit)["skew"])
##     curve(fitDist(x), col = "red", add = TRUE)
##     ## Hellinger's distance
##     shdf <- function(x) (sqrt(trueDist(x)) - sqrt(fitDist(x)))^2
##     hd <- try(integrate(shdf, lower = -Inf, upper = Inf)$value)
##     if(!inherits(hd, "try-error")){
##       norm.rlt[, "hd"] <- 0.5 * hd
##     } else {
##       norm.rlt[, "hd"] <- NaN
##     }
##     ## Kullback-Leibner divergence
##     kldf <- function(x) trueDist(x) * (log(trueDist(x)) - log(fitDist(x)))
##     kld <- try(integrate(kldf, lower = -20, upper = 20)$value)
##     if(!inherits(kld, "try-error")){
##       norm.rlt[, "kld"] <- kld
##     } else {
##       norm.rlt[, "kld"] <- NaN
##     }
##     ## Mean integrated squared error
##     misef <- function(x) (trueDist(x) - fitDist(x))^2
##     mise <- try(integrate(misef, lower = -Inf, upper = Inf)$value)
##     if(!inherits(mise, "try-error")){
##       norm.rlt[, "mise"] <- mise
##     } else {
##       norm.rlt[, "mise"] <- NaN
##     }
##     norm.rlt[, "omega"] <- abs(model$omega - coef(norm.fit)["omega"])
##     norm.rlt[, "alpha1"] <- abs(model$alpha - coef(norm.fit)["alpha1"])
##     norm.rlt[, "beta1"] <- abs(model$beta - coef(norm.fit)["beta1"])
##   }
##   ## Student t distribution
##   ## ----------------------
##   std.fit <- try(garchFit(data = sim, cond.dist = "sstd", trace = FALSE,
##                            include.mean = FALSE))
##   if(!inherits(std.fit, "try-error")){
##     fitDist <- function(x) dsstd(x, mean = 0, sd = 1,
##                                  nu = coef(std.fit)["shape"],
##                                  xi = coef(std.fit)["skew"])
##     curve(fitDist(x), col = "green", add = TRUE)
##     ## Hellinger's distance
##     shdf <- function(x) (sqrt(trueDist(x)) - sqrt(fitDist(x)))^2
##     hd <- try(integrate(shdf, lower = -Inf, upper = Inf)$value)
##     if(!inherits(hd, "try-error")){
##       std.rlt[, "hd"] <- 0.5 * hd
##     } else {
##       std.rlt[, "hd"] <- NaN
##     }
##     ## Kullback-Leibner divergence
##     kldf <- function(x) trueDist(x) * (log(trueDist(x)) - log(fitDist(x)))
##     kld <- try(integrate(kldf, lower = -20, upper = 20)$value)
##     if(!inherits(kld, "try-error")){
##       std.rlt[, "kld"] <- kld
##     } else {
##       std.rlt[, "kld"] <- NaN
##     }
##     ## Mean integrated squared error
##     misef <- function(x) (trueDist(x) - fitDist(x))^2
##     mise <- try(integrate(misef, lower = -Inf, upper = Inf)$value)
##     if(!inherits(mise, "try-error")){
##       std.rlt[, "mise"] <- mise
##     } else {
##       std.rlt[, "mise"] <- NaN
##     }
##     std.rlt[, "omega"] <- abs(model$omega - coef(std.fit)["omega"])
##     std.rlt[, "alpha1"] <- abs(model$alpha - coef(std.fit)["alpha1"])
##     std.rlt[, "beta1"] <- abs(model$beta - coef(std.fit)["beta1"])
##   }
##   ## Generalised Error distribution
##   ## ------------------------------
##   ged.fit <- try(garchFit(data = sim, cond.dist = "sged", trace = FALSE,
##                            include.mean = FALSE))
##   if(!inherits(ged.fit, "try-error")){
##     fitDist <- function(x) dsged(x, mean = 0, sd = 1,
##                                nu = coef(ged.fit)["shape"],
##                                xi = coef(ged.fit)["skew"])
##   curve(fitDist(x), col = "orange", add = TRUE)
##     ## Hellinger's distance
##     shdf <- function(x) (sqrt(trueDist(x)) - sqrt(fitDist(x)))^2
##     hd <- try(integrate(shdf, lower = -Inf, upper = Inf)$value)
##     if(!inherits(hd, "try-error")){
##       ged.rlt[, "hd"] <- 0.5 * hd
##     } else {
##       ged.rlt[, "hd"] <- NaN
##     }
##     ## Kullback-Leibner divergence
##     kldf <- function(x) trueDist(x) * (log(trueDist(x)) - log(fitDist(x)))
##     kld <- try(integrate(kldf, lower = -20, upper = 20)$value)
##     if(!inherits(kld, "try-error")){
##       ged.rlt[, "kld"] <- kld
##     } else {
##       ged.rlt[, "kld"] <- NaN
##     }
##     ## Mean integrated squared error
##     misef <- function(x) (trueDist(x) - fitDist(x))^2
##     mise <- try(integrate(misef, lower = -Inf, upper = Inf)$value)
##     if(!inherits(mise, "try-error")){
##       ged.rlt[, "mise"] <- mise
##     } else {
##       ged.rlt[, "mise"] <- NaN
##     }
##     ged.rlt[, "omega"] <- abs(model$omega - coef(ged.fit)["omega"])
##     ged.rlt[, "alpha1"] <- abs(model$alpha - coef(ged.fit)["alpha1"])
##     ged.rlt[, "beta1"] <- abs(model$beta - coef(ged.fit)["beta1"])
##   }
##   ## Scale Normal Mixture distribution
##   ## ---------------------------------
##   mnorm.fit <- try(cnmms(as.mgarch(sim), plot = "gradient", grid = 1000,
##                          verb = 0, tol = 1e-5))
##   if(!inherits(mnorm.fit, "try-error")){
##     print(summary.mgarch(mnorm.fit))
##     fitDist <- function(x) dmsnorm(x, xi = mnorm.fit$beta["xi"],
##                                    varmix = disc(pt = mnorm.fit$mix$pt,
##                                                  pr = mnorm.fit$mix$pr))
##     curve(fitDist(x), add = TRUE, col = "blue")
##     ## Hellinger's distance
##     shdf <- function(x) (sqrt(trueDist(x)) - sqrt(fitDist(x)))^2
##     hd <- try(integrate(shdf, lower = -Inf, upper = Inf)$value)
##     if(!inherits(hd, "try-error")){
##       mnorm.rlt[, "hd"] <- 0.5 * hd
##     } else {
##       mnorm.rlt[, "hd"] <- NaN
##     }
##     ## Kullback-Leibner divergence
##     kldf <- function(x) trueDist(x) * (log(trueDist(x)) - log(fitDist(x)))
##     kld <- try(integrate(kldf, lower = -20, upper = 20)$value)
##     if(!inherits(kld, "try-error")){
##       mnorm.rlt[, "kld"] <- kld
##     } else {
##       mnorm.rlt[, "kld"] <- NaN
##     }
##     ## Mean integrated squared error
##     misef <- function(x) (trueDist(x) - fitDist(x))^2
##     mise <- try(integrate(misef, lower = -Inf, upper = Inf)$value)
##     if(!inherits(mise, "try-error")){
##       mnorm.rlt[, "mise"] <- mise
##     } else {
##       mnorm.rlt[, "mise"] <- NaN
##     }
##     mnorm.rlt[, "omega"] <- abs(model$omega - coef(mnorm.fit)["omega"])
##     mnorm.rlt[, "alpha1"] <- abs(model$alpha - coef(mnorm.fit)["alpha1"])
##     mnorm.rlt[, "beta1"] <- abs(model$beta - coef(mnorm.fit)["beta1"])
##   }
##   legend("topleft", legend = c("real", "norm", "t", "ged", "mnorm"),
##          col = c("black", "red", "green", "orange", "blue"), bty = "n",
##          lty = 1)
##   rbind(norm.rlt, std.rlt, ged.rlt, mnorm.rlt)
## }
