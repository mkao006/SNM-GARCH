# source("functions.R")

######################################################################
## DEM2GBP
######################################################################

## Example for the model

# package(fGarch)

tdg <- as.numeric(data.matrix(dem2gbp))[1:100]
tmdg <- tdg
class(tmdg) <- "mgarch"
# (dg.mg <- cnmms(tmdg, plot = "gradient", grid = 100, verb = 4))


# dg.t0 <- garchFit(data = tdg, cond.dist = "norm")
# dg.t1 <- garchFit(data = tdg, cond.dist = "ged")
# dg.t2 <- garchFit(data = tdg, cond.dist = "std")



######################################################################
## check
######################################################################

check <- logd.mgarch(1:5/10, pt = c(1, 2), beta = rep(0.5, 4),
                     which = rep(1, 4))


## Check the derivatives
# incremt <- 1e-10

## dl/dalpha0
# (logd.mgarch(nyse[1:5], pt = c(1, 2), beta = c(0.5 + incremt, 0.5, 0.5, 0.5),
#             which = c(1, 0, 0, 0))$ld -
#     logd.mgarch(nyse[1:5], pt = c(1, 2), beta = rep(0.5, 4),
#                 which = c(1, 0, 0, 0))$ld)/incremt

# logd.mgarch(nyse[1:5], pt = c(1, 2), beta = rep(0.5, 4),
#            which = c(0, 1, 0, 0))$db1[,,1]

test.diff = function(x=tdg[1:5], beta=c(.3,.4,.5,.7), incremt=1e-10, pt=1:2) {
  index0 = rep(0,4)
  for(i in 1:4) {
    index = index0
    index[i] = 1
    d1 = (logd.mgarch(x, pt=pt, beta=beta+index*incremt, which=c(1,0,0,0))$ld -
          logd.mgarch(x, pt=pt, beta=beta, which=c(1,0,0,0))$ld ) / incremt
    d2 = logd.mgarch(x, pt=pt, beta=beta, which=c(0,1,0,0))$db1[,,i]
    print(max(abs(d1-d2)))
  }
}

# (logd.mgarch(tdg[1:5], pt = c(1, 2), beta = c(0.5 + incremt, 0.5, 0.5, 0.5),
#             which = c(1, 0, 0, 0))$ld -
#     logd.mgarch(tdg[1:5], pt = c(1, 2), beta = rep(0.5, 4),
#                 which = c(1, 0, 0, 0))$ld)/incremt
# 
# logd.mgarch(tdg[1:5], pt = c(1, 2), beta = rep(0.5, 4),
#             which = c(0, 1, 0, 0))$db1[,,1]
# 
# 
# 
# ## dl/dalpha1
# (logd.mgarch(tdg[1:5], pt = c(1, 2), beta = c(0.5, 0.5 + incremt, 0.5, 0.5),
#             which = c(1, 0, 0, 0))$ld -
#     logd.mgarch(tdg[1:5], pt = c(1, 2), beta = rep(0.5, 4),
#                 which = c(1, 0, 0, 0))$ld)/incremt
# 
# 
# logd.mgarch(tdg[1:5], pt = c(1, 2), beta = rep(0.5, 4),
#             which = c(0, 1, 0, 0))$db1[,,2]
# 
# 
# ## dl/dbeta1
# (logd.mgarch(tdg[1:5], pt = c(1, 2), beta = c(0.5, 0.5, 0.5 + incremt, 0.5),
#             which = c(1, 0, 0, 0))$ld -
#     logd.mgarch(tdg[1:5], pt = c(1, 2), beta = rep(0.5, 4),
#                 which = c(1, 0, 0, 0))$ld)/incremt
# 
# 
# logd.mgarch(tdg[1:5], pt = c(1, 2), beta = rep(0.5, 4),
#             which = c(0, 1, 0, 0))$db1[,,3]
# 
# 
# ## dl/dsigma0
# (logd.mgarch(tdg[1:5], pt = c(1, 2), beta = c(0.5, 0.5, 0.5, 0.5 + incremt),
#             which = c(1, 0, 0, 0))$ld -
#     logd.mgarch(tdg[1:5], pt = c(1, 2), beta = rep(0.5, 4),
#                 which = c(1, 0, 0, 0))$ld)/incremt
# 
# 
# logd.mgarch(tdg[1:5], pt = c(1, 2), beta = rep(0.5, 4),
#             which = c(0, 1, 0, 0))$db1[,,4]
# 
