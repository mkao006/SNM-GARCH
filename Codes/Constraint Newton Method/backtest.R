library(tseries)
spData <- ts(get.hist.quote("SPY", quote = c("Close"), start = "2006-01-01", end = "2008-12-31"), frequency = 365,
             start = c(2006, 3))

mean(diff(log(spData)))


cr <- cReturns(spData)
mean(cr[[1]])

spa <- dataSnoop(spData,bSamples=3,test="SPA")

z <- rep(c(1,1,0,0,0,0), 100)
res <- benchmark(fun1(z), fun2(z),
                  fun1c(z), fun2c(z),
                  funRcpp(z),
                  columns=c("test", "replications", "elapsed",
                            "relative", "user.self", "sys.self"),
                  order="relative",
                  replications=1000)
print(res)
