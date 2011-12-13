# ======================================== #
# Multivariate Scale Mixture Normal #
# ======================================== #
library(mvtnorm)

## constructor
msmn = function(x, sigma=NULL){
    if(!is.numeric(x)) stop("x must be numeric")
    if(is.vector(x)) x = matrix(x)
    if(class(x)!="matrix") stop("x must be vector or matirx")
    if(ncol(x)==1) sigma=matrix(1)
    else sigma=cov(x)
    data = list(x=x, sigma=sigma, sigma.inv=solve(sigma))
    class(data) = "msmn"
    data
}


## update sigma
setSigma = function(x, sigma){
    if(ncol(sigma)!=nrow(sigma)) stop("sigma must be squared matrix")
    if(ncol(sigma)!=dim(x)) stop("dimension of sigma must be same as dim(x)")
    if(class(try(solve(sigma)))=="try-error") stop("sigma is singular")
    x$sigma = sigma
    x$sigma.inv = solve(sigma)
    x
}

'[.msmn' = function(obj, i) msmn(x=obj$x[i,], sigma=obj$sigma)
mean.msmn = function(obj) apply(obj$x, 2, mean)
dim.msmn = function(obj) ncol(obj$x)
length.msmn = function(obj) nrow(obj$x)
cov.msmn = function(obj) cov(obj$x)


## information of 1-D msmn object,
## include sample mean, sample variance, fitted mean, fitted variance
info.msmn = function(x, beta, mix){
    if(dim.msmn(x)==1){
        tab = matrix(ncol=2, nrow=2)
        dimnames(tab) = list(c("Sample", "Estimate"), c("mean", "var"))
        tab[1,] = c(mean.msmn(x), cov.msmn(x))
        tab[2,] = c(beta, sum(mix$pt*mix$pr))
    } else{
        stop("dimension not 1")
    }
    tab
}

## generate random number of msmn object
rmsmn = function(n=100, mu=3,
mix=dden(pr=c(.7,.3), pt=c(3, 8)), sigma=matrix(1)){
    data = numeric(0);
    for(i in 1:length(mix$pr)){
        ni = ceiling((n*mix$pr[i]))
        data = rbind(data, rmvnorm(n=ni, mean=mu, sigma=mix$pt[i]*sigma))
    }
    data = data[sample(1:n),]
    msmn(data, sigma)
}


## cumulative distribution function of msmn
pmsmn = function(x, beta, mix){
    n = length(x)
    y = numeric(n)
    for(i in 1:n)
        y[i] = sum(mix$pr*pnorm(x[i], beta, sqrt(mix$pt)))
    y
}


## quantile function of msmn object
qmsmn = function(p, beta, mix){
    n = length(p)
    y = numeric(n)
    SD = sqrt(sum(mix$pr*mix$pt))
    for(i in 1:n){
        r = beta + c(-100*SD, 100*SD)
        while(TRUE){
            m = pmsmn(mean(r), beta, mix)
            if(abs(m-p[i])<1e-15) break
            if(p[i]<m) r = c(r[1], mean(r))
            else r = c(mean(r), r[2])
        }
        y[i] = mean(r)
    }
    y
}

## density of msmn object
dmsmn = function(x, mu, mix, sigma=x$sigma){
    d = numeric(length(x))
    for(i in 1:length(mix$pr)){
        temp = mix$pr[i]*dmvnorm(x$x, mean=mu, sigma=mix$pt[i]*sigma)
        if(all(temp==0)){
            print(paste("(pr, pt) =(", mix$pr[i], ",", mix$pt[i], ")"))
            warning("precision error?")
        }
        d = d + temp
    }
    d
}


## range of mixing variances
range.msmn = function(x, ...){
    if(dim.msmn(x)==1){
        lower = cov.msmn(x)*1e-8
        upper = cov.msmn(x)*1e8
    }
    else{
        lower = 1e-8
        upper = 1e8
    }
    c(lower, upper)
}


## initial value
initial.msmn = function(x, beta=NULL, mix=NULL, kmax=NULL) {
    if(is.null(beta)) beta = mean(x)
    if(is.null(kmax)) kmax = 10
    if(is.null(mix) || is.null(mix$pt)){
        pt.rep = seq(range(x, beta)[1], range(x, beta)[2], length.out=10)
        mix = dden(unique(quantile(pt.rep, p=seq(0,1,len=kmax), type=1)))
    }
    list(beta=beta, mix=mix)
}


## validation
valid.msmn = function(x, beta, mix){
    all(mix$pt > 0) &&
    all(mix$pt >= range(x)[1]) &&
    all(mix$pt <= range(x)[2])
}


## log likelihood expression
logd.msmn = function(x, beta, pt, which=c(1,0,0,0)){
    p = length(beta)
    dl = vector("list", 4)
    d = mahalanobis(x$x, beta, x$sigma)
    names(dl) = c("ld","db1","dt1","dt2")
    if(which[1] == 1){
        dl$ld = outer(rep(-1/2, nrow(x$x)),
        log((2*pi*pt)^p*det(x$sigma)), "*") -
            outer(d, 1/(2*pt), "*")
    }
    if(which[2] == 1) {
        dl$db1 = array(dim=c(nrow(x$x), length(pt), length(beta)))
        x.temp = t(x$sigma.inv%*%t(sweep(x$x, 2, beta)))
        for(i in 1:length(pt))
            dl$db1[,i,] = x.temp/pt[i]
    }
    if(which[3] == 1)
        dl$dt1 = outer(rep(p, nrow(x$x)), -1/(2*pt), "*") +
            outer(d, 1/(2*pt^2), "*")
    if(which[4] == 1)
        dl$dt2 = outer(rep(p, nrow(x$x)), 1/(2*pt^2), "*") -
            outer(d, 1/(pt^3), "*")
    dl
}


## find new sigma during the EM algorithm
newSigma = function(x, beta, mix, sigma){
    p = matrix(nrow=length(x), ncol=length(mix$pt))
    z = dmsmn(x, beta, mix, sigma)
    for(k in 1:ncol(p))
        p[,k] = mix$pr[k]*dmvnorm(x$x, mean=beta, sigma=mix$pt[k]*sigma)
    p = p/z
    sigma.new = 0
    p.new = sweep(p, 2, mix$pt, "/")
    if(is.vector(p.new)) p.new = matrix(p.new)
    else p.new = t(p.new)
    x.new = t(sweep(x$x, 2, beta))
    for(i in 1:nrow(p)){
        sigma.new = sigma.new + x.new[,i]%*%t(x.new[,i])*sum(p.new[,i])
    }
    sigma.new/length(x)
}


## update sigma using the EM algorithm
updateSigma = function(x, beta, mix, sigma, tol=1e-6){
    while(TRUE){
        sigma.new = newSigma(x, beta, mix, sigma)
        z = abs(sum(log(dmsmn(x, beta, mix, sigma.new))) -
        sum(log(dmsmn(x, beta, mix, sigma))))
        sigma = sigma.new
        if(z<tol) break
    }
    sigma
}


## find semiparametric maximum likelihood,
## using the combination of the EM algorithm and the CNM-based algorithm
msmnFit = function(x, phi=initial(x), cnmFUN=cnmms, plot="g", tol=1e-10){
    if(dim.msmn(x)==1){
        phi = cnmFUN(x, plot=plot)
        phi$mix = truncMix(x, phi$beta, phi$mix)
    } else{
        phi$sigma = x$sigma
        while(TRUE){
            print(phi)
            sigma.new = updateSigma(x, phi$beta, phi$mix, phi$sigma, tol=tol)
            x = setSigma(x, sigma.new)
            print("Sigma Updated :")
            print(x$sigma)
            fit = cnmFUN(x, init=list(beta=phi$beta, mix=phi$mix), plot=plot)
                                        #fit = cnmFUN(x, plot=plot)
            phi.new = list(beta=fit$beta, mix=fit$mix, sigma=sigma.new)
            print("(beta, mix) Upated")
            z = abs(sum(log(dmsmn(x, phi$beta, phi$mix, phi$sigma))) -
            sum(log(dmsmn(x, phi.new$beta, phi.new$mix, phi.new$sigma))))
            print(cat("difference of ll is", z, "\n"))
            phi = phi.new
            phi$mix = truncMix.msmn(x, phi$beta, phi$mix)
            if(max(phi$mix$pt)==max(range(x)) || min(phi$mix$pt)==min(range(x)))
                if(max(phi$mix$pt)/min(phi$mix$pt)<max(range(x))/min(range(x))){
                    phi$sigma = phi$sigma*min(phi$mix$pt)
                    phi$mix$pt = phi$mix$pt/min(phi$mix$pt)
                }
            print(phi)
            if(z<tol) break
        }
        phi$sigma = phi$sigma*min(phi$mix$pt)
        phi$mix$pt = phi$mix$pt/min(phi$mix$pt)
    }
    phi$ll.max = sum(log(dmsmn(x, phi$beta, phi$mix, phi$sigma)))
    phi
}


## d=mahalanobis distance
## plot mahalanobis distance for bivariate normal distribution
plot.ellipse = function(d, mu, sigma, col="red", lwd=3, add=FALSE, ...){
    n=1e3
    xy = matrix(ncol=length(mu), nrow=n)
    theta = seq(0, 2*pi, length.out=n)
    eg = eigen(sigma)
    for(i in 1:n)
        xy[i,] = d * eg$vectors %*% (sqrt(eg$values)*c(cos(theta[i]),
          sin(theta[i]))) + mu
    if(add) lines(xy[,1], xy[,2], col=col, lwd=lwd, ...)
    else plot(xy, type="l", col=col, lwd=lwd, ...)
}

## plot for msmn object, only works for 1-D and 2-D
## generate density plot in 1-D case,
## generate mahalanobis distance plot in 2-D case
msmn.plot = function(mu, mix, d=1, sigma=NULL, add=FALSE, each=TRUE,
lwd=2, leg=TRUE, ...){
    if(d==1){
        if(each && lwd==1) stop("lwd cannot be 1 when each=TRUE")
        if(length(mu)!=1) stop("mu is not 1-D")
        x.range = mu + c(-1, 1)*sqrt(max(mix$pt))*20
        x = seq(x.range[1], x.range[2], length.out=1e5)
        k = length(mix$pt)
        y = matrix(nrow=length(x), ncol=k)
        for(i in 1:k)
            y[,i] = mix$pr[i]*dnorm(x, mu, sqrt(mix$pt[i]))
        y = cbind(y, apply(y, 1, sum))
        if(!add){
            plot.new(); plot.window(xlim=x.range, ylim=c(0, max(y)))
            box(); axis(1); axis(2)
        }
        lines(x, y[,k+1], lwd=lwd, ...)
        if(each){
            color = rainbow(k)
            for(i in 1:k)
                lines(x, y[,i], lty=2, lwd=2, col=color[i])
            if(leg){
                legend("topright", legend=paste("(",signif(mix$pt,4),", ",
                                   signif(mix$pr, 4),")",sep=""),
                       title=expression(theta),
                       lwd=lwd-1, col=color, lty=2)
            }
        }
    }else if(d==2){
        if(length(mu)!=2) stop("mu is not 2-D")
        k = length(mix$pt)
        color = rainbow(k)
        for(i in k:1){
            plot.ellipse(d=1, mu, mix$pt[i]*sigma, col=color[i], lty=2, add=add)
            add = TRUE
        }
        if(leg){
            legend("topright", legend=paste("(",signif(mix$pt,4),", ",
                               signif(mix$pr, 4),")",sep=""),
                   title=expression(theta),
                   lwd=2, col=color, lty=2)
        }
    }
}
