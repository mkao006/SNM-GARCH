# ================ #
# Discrete density #
# ================ #

# pt    Points
# pr    Probabilities at the points

dden = function(pt, pr=1) {
  if (is.null(pt) ) d = list(pt=NULL, pr=NULL)
  else {
    k = max(length(pt), length(pr), na.rm=TRUE)
    pt = rep(pt, len=k)
    if(is.null(pr)) pr = rep(1, length(pt))
    else pr = rep(pr, len=k)
    d = list(pt=pt, pr=pr/sum(pr))
  }
  class(d) = "dden"
  d
}

is.null.dden = function (d) is.null(d$pt)

is.dden = function (d) any(class(d) == "dden")

print.dden = function (d, ...) {
  if (is.null(d)) b = matrix(nrow=0, ncol=2)
  else b = cbind(d$pt, d$pr)
  dimnames(b) = list(NULL, c("pt", "pr"))
  print(b, ...)
}

plot.dden = function (d, ...) {
  ylim = c(0, max(d$pr))
  plot(d$pt, d$pr, type="h", col="blue", lwd=2, ylim=ylim,
       xlab="", ylab="Probability", ...)
}

# Sort

sort.dden = function(d) {
  if( is.null(d) ) return(d)
  index = order(d$pt)
  d$pt = d$pt[index]
  d$pr = d$pr[index]
  d
}

is.unsorted.dden = function(d) is.unsorted(d$pt)

# Unique

unique.dden = function(d, prec=0) {
  if( length(d$pt) == 1 ) return(d)
  if( is.unsorted.dden(d) ) d = sort.dden(d)
  prec = rep(prec, len=2)
  if ( all(prec < 0) ) return(d)
  pt2 = pt = d$pt
  pr2 = pr = d$pr
  j  = pr <= prec[2]
  pt = pt[!j]
  pr = pr[!j]
  index = 0
  repeat {
    if( length(pt) == 0 ) break
    j = abs(pt[1] - pt) <=  prec[1]
    index = index + 1
    pt2[index] = weighted.mean(pt[j], pr[j])
    pr2[index] = sum( pr[j] )
    pt = pt[!j]
    pr = pr[!j]
  }
  dden(pt=pt2[1:index], pr=pr2[1:index])
}

