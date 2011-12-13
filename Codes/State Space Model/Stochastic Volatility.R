######################################################################
## Example 6.18 Analysis of the New Tork Stock Exchange Returns
######################################################################

## Load the data and functions
load("tsa3.rda")

## Plot the original and the transformed time series
par(mfrow = c(2, 1))
plot(nyse)
plot(log(nyse^2))

## Stochastic volatility model of NYSE return
y <- log(nyse^2)
num <- length(y)

## Initiate the parameters
phi0 = 0
phi1 = 0.95
sQ = 0.2
alpha = mean(y)
sR0 = 1
mu1 = -2
sR1 = 2
init.par = c(phi0, phi1, sQ, alpha, sR0, mu1, sR1)

## The likelihood estimation based on Kalman filter
Linn <- function(para){
    phi0 = para[1]
    phi1 = para[2]
    sQ = para[3]
    alpha = para[4]
    sR0 = para[5]
    mu1 = para[6]
    sR1 = para[7]
    sv <- SVfilter(num, y, phi0, phi1, sQ, alpha, sR0, mu1, sR1)
    sv$like
}

## Maximise the likelihood (MLE)
(est <- optim(init.par, Linn, NULL, method = "BFGS",
              hessian = TRUE, control = list(trace = 1, REPORT = 1)))
SE <- sqrt(diag(solve(est$hessian)))
u <- cbind(estimates = est$par, SE)
rownames(u) <- c("phi0", "phi1", "sQ", "alpha", "sigv0", "mu1", "sigv1")
u


# Graphics (need filters at the estimated parameters)
phi0=est$par[1]
phi1=est$par[2]
sQ=est$par[3]
alpha=est$par[4]
sR0=est$par[5]
mu1=est$par[6]
sR1=est$par[7]

## The filtered series, I suspect that the mean alpha was not added
## back
sv = SVfilter(num,y,phi0,phi1,sQ,alpha,sR0,mu1,sR1)
mysv = MySVfilter(num,y,phi0,phi1,sQ,alpha,sR0,mu1,sR1)

## The fit
#Time = 801:1000
Time = 1:2000
plot(Time, y[Time], type="l", main="log(Squared NYSE Returns)",
     ylim = c(-20, 2), lwd = 2)
lines(Time, sv$xp[Time] + alpha,type="l", main="Predicted log-Volatility",
      ylab="", xlab="", col = "red", lty = 2, lwd = 2)
lines(Time, sv$xp[Time] + alpha + 2*sqrt(sv$Pp[Time]),
      lty="dashed", col = "orange")
lines(Time, sv$xp[Time] + alpha - 2*sqrt(sv$Pp[Time]),
      lty="dashed", col = "orange")



res <- y - sv$xp - alpha
hist(res, freq = FALSE, breaks = 300, ylim = c(0, 0.35))
#densities plot (f is chi-sq, fm is fitted mixture)
x = seq(-15,6,by=.01)
f = exp(-.5*(exp(x)-x))/(sqrt(2*pi))
f0 = exp(-.5*(x^2)/sR0^2)/(sR0*sqrt(2*pi))
f1 = exp(-.5*(x-mu1)^2/sR1^2)/(sR1*sqrt(2*pi))
fm = (f0+f1)/2
lines(x, f, lwd = 2)
lines(x, fm, lty = 2, lwd = 2, col = "red")
lines(x, f0, col = "red", lty = 2)
lines(x, f1, col = "red", lty = 2)


## The total fit
plot(y, type="l", main="log(Squared NYSE Returns)",
     ylim = range(y), lwd = 1)
lines(sv$xp + alpha,type="l", main="Predicted log-Volatility",
      ylab="", xlab="", col = "red", lty = 2)
lines(sv$xp + alpha + 2*sqrt(sv$Pp[Time]), lty="dashed",
      lwd = 2, col = "orange")
lines(sv$xp + alpha - 2*sqrt(sv$Pp[Time]), lty="dashed",
      lwd = 2, col = "orange")

## It appears that the noise to signal ratio is high


## Lets try it with the mixture density

MyKfilter0 <- function(num,y,A,mu0,Sigma0,Phi,cQ,cR){
  #
  # NOTE: must give cholesky decomp: cQ=chol(Q), cR=chol(R)
 Q=t(cQ)%*%cQ
 R=t(cR)%*%cR
   # y is num by q  (time=row series=col)
   # A is a q by p matrix
   # R is q by q
   # mu0 is p by 1
   # Sigma0, Phi, Q are p by p
 Phi=as.matrix(Phi)
 pdim=nrow(Phi)
 y=as.matrix(y)
 qdim=ncol(y)
 xp=array(NA, dim=c(pdim,1,num))         # xp=x_t^{t-1}
 Pp=array(NA, dim=c(pdim,pdim,num))      # Pp=P_t^{t-1}
 xf=array(NA, dim=c(pdim,1,num))         # xf=x_t^t
 Pf=array(NA, dim=c(pdim,pdim,num))      # Pf=x_t^t
 innov=array(NA, dim=c(qdim,1,num))      # innovations
 sig=array(NA, dim=c(qdim,qdim,num))     # innov var-cov matrix
# initialize (because R can't count from zero)
   x00=as.matrix(mu0, nrow=pdim, ncol=1)
   P00=as.matrix(Sigma0, nrow=pdim, ncol=pdim)
   xp[,,1]=Phi%*%x00
   Pp[,,1]=Phi%*%P00%*%t(Phi)+Q
     sigtemp=A%*%Pp[,,1]%*%t(A)+R
   sig[,,1]=(t(sigtemp)+sigtemp)/2     # innov var - make sure it's symmetric
   siginv=solve(sig[,,1])
   K=Pp[,,1]%*%t(A)%*%siginv
   innov[,,1]=y[1,]-A%*%xp[,,1]
   xf[,,1]=xp[,,1]+K%*%innov[,,1]
   Pf[,,1]=Pp[,,1]-K%*%A%*%Pp[,,1]
   sigmat=as.matrix(sig[,,1], nrow=qdim, ncol=qdim)
   like = log(det(sigmat)) + t(innov[,,1])%*%siginv%*%innov[,,1]   # -log(likelihood)
########## start filter iterations ###################
 for (i in 2:num){
   if (num < 2) break
  xp[,,i]=Phi%*%xf[,,i-1]
  Pp[,,i]=Phi%*%Pf[,,i-1]%*%t(Phi)+Q
     sigtemp=A%*%Pp[,,i]%*%t(A)+R
   sig[,,i]=(t(sigtemp)+sigtemp)/2     # innov var - make sure it's symmetric
   siginv=solve(sig[,,i])
  K=Pp[,,i]%*%t(A)%*%siginv
  innov[,,i]=y[i,]-A%*%xp[,,i]
  xf[,,i]=xp[,,i]+K%*%innov[,,i]
  Pf[,,i]=Pp[,,i]-K%*%A%*%Pp[,,i]
    sigmat=as.matrix(sig[,,i], nrow=qdim, ncol=qdim)
  like= like + dnormmix(mix, x, log = TRUE)
  }
    list(xp=xp,Pp=Pp,xf=xf,Pf=Pf,like=like,innov=innov,sig=sig,Kn=K)
}


## The likelihood estimation based on Kalman filter
Linn <- function(para, mix){
    phi0 = para[1]
    phi1 = para[2]
    sQ = para[3]
    alpha = para[4]
    sR0 = para[5]
    mu1 = para[6]
    sR1 = para[7]
    sv <- Kfilter0(num, y, phi0, phi1, sQ, alpha, sR0, mu1, sR1, mix)
    sv$like
}

iter = 10
ss.est <- matrix(0, nc = 7, nr = iter)
for(i in 1:iter){
    cnm.normmix(


## Maximise the likelihood (MLE)
(est <- optim(init.par, Linn, NULL, method = "BFGS",
              hessian = TRUE, control = list(trace = 1, REPORT = 1)))
SE <- sqrt(diag(solve(est$hessian)))
u <- cbind(estimates = est$par, SE)
rownames(u) <- c("phi0", "phi1", "sQ", "alpha", "sigv0", "mu1", "sigv1")
u
