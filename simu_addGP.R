setwd('~/Downloads/test')
source('gp.R')

###

# fsimu <- function(x) return(2 * sum(x[1:10]))

# fsimu <- function(x) return(10 * sin(pi * x[1]) + 10 * sin(pi * x[2]) + 20 * (x[3] -0.5)**2 + 10 * x[7])

fsimu <- function(x) return(10 * sin(pi * x[1]*x[2]) + 10 * cos(pi * (x[3]*x[4] + x[5])) + 20 * (x[6] -0.5)**2 + 10 * x[7])

n <- 100
p <- 50

set.seed(100)
x <- matrix(runif(n*p), n, p)
y <- apply(x, 1, fsimu) + rnorm(n, sd = 1)
x.train <- x
y.train <- y
x.test  <- matrix(runif(2*n*p), 2*n, p)
ytst <- apply(x.test, 1, fsimu)


###
ncomp <- 10
max.rank <- min(n, max(40, 2 * log(n)^2))
prox.lamsq <- c(0.99, 0.94, 0.88, 0.81, 0.73, 0.67)
Rsq.rhosq <- c(0.00, 0.25, 0.50, 0.70, 0.85, 0.90, 0.99)
# Rsq.rhosq <- c(0.00, 0.05, 0.25, 0.50, 0.70, 0.85, 0.90, 0.99)

## addGP - using Toggle variable selection
o1 <- add.gp(x.train, y.train, nsweep=1000, ncomp=ncomp, max.rank=max.rank,
         prox=prox.lamsq, Rsq=Rsq.rhosq)

## addGP - using pMTM variable selection
o2 <- add.gp(x.train, y.train, nsweep=200, ncomp=ncomp, max.rank=max.rank,
         prox=prox.lamsq, Rsq=Rsq.rhosq, pmtm.budget=5*ncomp, varp.update=0)
# pmtm.budget=floor(log(p)**2)


###
out <- o2
summary(out, burn=0.3)

agp.pred <- predict.addGP(out, x.test=x.test, burn=0.3, nsamp=200, fsamp=TRUE)
par(mfrow = c(1,3), mar = rep(5,5,4,2))
plot(agp.pred$sig2, ty = "l", xlab = "", ylab = "")
title(xlab = "Thinned iterations", ylab = "sigSq")
plot(ytst, rowMeans(agp.pred$f.samp), xlab = "", ylab = "")
abline(0,1, col = "red", lwd = 2)
title(xlab = "y.pred", ylab = "y.true", main = paste("rmse = ", round(rmse(rowMeans(agp.pred$f.samp), ytst), 2)))

plot(out$var.p, ty = "h", xlab = "predictor", ylab = "propensity score")
plot(rowSums(o2$active.store), ty = "l")