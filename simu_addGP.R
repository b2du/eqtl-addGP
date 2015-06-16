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
agp.rmse <- rmse(ytst, rowMeans(agp.pred$f.samp))
plot(ytst, rowMeans(agp.pred$f.samp), xlab = "y.pred", ylab = "y.true", main = paste("rmse = ", round(agp.rmse, 2)))
abline(0,1, col = "red", lwd = 2)
# plot(out$var.p, ty = "h", xlab = "predictor", ylab = "propensity score")

par(mfrow = c(1,2), mar = rep(5,5,4,2))
plot(agp.pred$sig2, ty = "l", xlab = "Thinned iterations", ylab = "sigSq")
plot(rowSums(o2$active.store), ty = "l", xlab = "MCMC iteration", ylab = "# active components")


###
require(BayesTree)
bt <- bart(x.train, as.numeric(y.train), x.test, sigest = 1, verbose = FALSE)
y.bt.test.hat.mean <- bt$yhat.test.mean
y.bt.test.hat.lowr <- apply(bt$yhat.test, 2, quantile, p = .025)
y.bt.test.hat.uppr <- apply(bt$yhat.test, 2, quantile, p = .975)

require(randomForest)
rfo <- randomForest(x.train, y.train)
y.rfo.test.hat.mean <- predict(rfo, newdata = x.test)

require(glmnet)
lass <- cv.glmnet(x.train, y.train)
lass <- glmnet(x.train, y.train, lambda = lass$lambda.min)
y.test.lass.hat <- predict(lass, newx = x.test)

cat("addGP RMSE =", round(agp.rmse, 4),"\n")
cat("BART RMSE =", bart.rmse <- round(rmse(ytst, y.bt.test.hat.mean), 4),"\n")
cat("RF RMSE =", rf.rmse <- round(rmse(ytst, y.rfo.test.hat.mean), 4),"\n")
cat("LASSO RMSE =", lm.rmse <- round(rmse(ytst, y.test.lass.hat), 4),"\n")
cat("NULL RMSE =", null.rmse <- round(rmse(ytst, mean(y.train)), 4),"\n")
