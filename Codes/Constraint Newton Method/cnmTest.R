library(lsei)
source("cnm.R")


## using the example in Yong's website
mix.efron <- normmix(mu = c(-10.9, -7, -4.9, -1.9, -1.1, 0, 2.4, 6.1),
                     pi = c(1.5, 1.3, 5.6, 12.3, 13.6, 60.9, 2.7, 2.2))
set.seed(1)

x <- rnormmix(n = 1000, mix = mix.efron)
fit <- cnm.normmix(x, tol = 1e-5)
abline(v = fit$mix$mu, col ="red")


thai <- data.frame(x=c(0:21, 23, 24),
                   freq=c(120, 64, 69, 72, 54, 35, 36, 25, 25, 19, 18,
                   18, 13, 4, 3, 6, 6, 5, 1, 3, 1, 2, 1, 2))
cnm.poismix(x=thai$x, w=thai$freq)



x <- rnorm(100)
y <- 2 + x + c(rnorm(100, 0, rep(c(2, 10), 50)))

cnm.normmix(y - x, tol = 1e-5)

## Normmix only works for different mean.
y <- c(rnorm(100, rep(c(8, 10), 50)), 1)
cnm.normmix(y)

myll <- function(x, a, b) dnorm(b - x * a)


myll <- function(x, a, b) sum(dnorm(b - x * a))
optimize(myll, c(0, 1), a = x, b = y)


## This is the usual maximum likelihood estimation
x <- rnorm(100)
y <- 2 + 0.8 * x + rnorm(100)

myll <- function(beta){
    -sum(dnorm(y - beta[1] - beta[2] * x))
}

test.mle <- optim(c(1, 1), myll)
test.mle
plot(x, y)
abline(test.mle$par[1], test.mle$par[2], col = "red")
test.lm <- lm(y ~ x)
abline(coef(test.lm)[1], coef(test.lm)[2], col = "blue")
test.glm <- lm(y ~ x)
abline(coef(test.glm)[1], coef(test.glm)[2], col = "green")


x1 <- rnorm(100)
x2 <- rnorm(100)
y <- 1 + 2 * x1 + 3 * x2 + rnorm(100)

## Extended to 2 variable
myll2 <- function(beta){
    -sum(dnorm(y - beta[1] - beta[2] * x1 - beta[3] * x2))
}
optim(c(0, 0, 0), myll2)

t3d.df <- data.frame(y, x1, x2)
library(boot)
myll3 <- function(data){
    ll <- -sum(dnorm(data[, 1] - beta[1] - beta[2] *
               data[, 2] - beta[3] * data[, 3]))
    optim(c(0, 0, 0), ll)$par
}

myll3 <- function(data){
    ll <- function(beta)
        -sum(dnorm(data[, 1] - beta[1] - beta[2] *
                     data[, 2] - beta[3] * data[, 3]))
    optim(c(0, 0, 0), ll)$par
}


boot(t3d.df, myll3, R = 1000)

######################################################################
## Discrete density estimation from mlogit.R
######################################################################
## Generate the dataset first
rpr <- c(0.27, 0.13, 0.068, 0.532)
rpt <- c(-3.245, -2.981, -0.705, 0.886)
sim.df <- rmlogit(k = 20, gi = 1, ni = 20, beta = 1,
                  pr = rpr, pt = rpt)

## Two different likelihood function, not sure why there are such
## instances
lf1.mlogit(sim.df, beta = 1, pt = rpt)
lf.mlogit(sim.df, beta = 1, pt = rpt)

## I suspect this is the likelihood function and the first order of
## the likelihood function
dlb1.mlogit(sim.df, 1, pt = rpt, order = c(0, 1))

######################################################################
## Common variance problem
######################################################################

## Generate the same data structure as presented in the paper
nsSim.df <- rcvp(ni = 5)

## The cvp2 function saves the data in the form of sufficient
## statistic
nss <- cvp2(nsSim.df)

logd.cvp2(nss, beta = 3, pt = 4)


######################################################################
## using the new cnm.R
######################################################################

## This part work through the new CNM function step by step to see
## what is going on in detail.

sim.df <- rcvp2()

for(i in 1:10){
    dev.new()
## The initial function sets where the CNM should start.
plotorder = 0
    verb = 0
k = i
maxit = 1000
grid = 100
tol = 1e-6
k = length(sim.df)
init <- initial.snpmle(sim.df, NULL)
beta = init$beta
nb = length(beta)
mix = init$mix
ll = init$ll
lll = -Inf
convergence = 1
for(i in 1:k){
    l <- logd(sim.df, inits$beta, inits$mix$pt,
                   which = c(1, 0, 0, 0))$ld
    ma = apply(l, 1, max)
    dmix = drop(exp(l - ma) %*% mix$pr) + 1e-100
    plotgrad(sim.df, beta, mix, ma, pch=19, order=plotorder)
    points(sim.df$mi, rep(0,length(sim.df$mi)), pch="|", cex=.5)
    rth = range(sim.df, beta)
    gridpoints = seq(rth[1], rth[2], length=grid)

    ## Finding the maxx of the gradient function
    g = maxgrad(sim.df, beta, dmix, ma, grid = gridpoints, tol = -Inf)
    gradient = max(g$grad)
    lll = ll
    mix1 = mix
    mix = dden(c(mix$pt, g$pt), c(mix$pr, rep(0, length(g$pt))))
    lpt = logd(sim.df, beta, mix$pt, which = c(1, 0, 0, 0))$ld
    dpt = pmin(exp(lpt - ma), 1e100)
    a = dpt/dmix - 2

    r = nnls(rbind(a, rep(1,length(mix$pt))), c(rep(0,nrow(a)),1))
    sol = r$x / sum(r$x)
    r = lsch(mix, beta, dden(mix$pt,sol), beta, sim.df, which=c(1,0,0))
    mix = collapse.snpmle(r$mix, beta, sim.df)
    ll = attr(mix, "ll")
    print.snpmle(verb, x, mix, beta, gradient)


      list(mix=mix, beta=beta, num.iterations=i, ll=ll[1],
       gradient=gradient, convergence=convergence)
}
}


######################################################################
## skewed normal distribution
######################################################################

snorm <- function(x, alpha) 2 * dnorm(x) * pnorm(alpha * x)

curve(dnorm(x), -5, 5, ylim = c(0, 1))

for(i in -1000:1000){
    curve(snorm(x, i/100), add = TRUE, col = rgb((i + 1000)/2000, 0, 0))
}


colSums(test2[test[, 1] == attr(test, "ui"), drop = FALSE])
