add.gp <- function(x.train, y.train, nsweep = 1e3, prox = c(0.99, 0.94, 0.88, 0.81, 0.73), Rsq = c(0.01, 0.25, 0.5, 0.7, 0.86, 0.99), lp.lam.f = 0, lp.rho.f = 0, ncomp, p.incl, p.ref = 0.1, a = 1, b = 1, max.rank = nrow(x.train), tol = c(1e-4,1e-6), burn = 0.1, pmtm.budget, varp.update, lpenalty){
    
    xnames <- dimnames(x.train)[[2]]
    
    n <- length(y.train)
    x.train <- matrix(x.train, nrow = n)
    p <- ncol(x.train)
    
    if(missing(ncomp)) ncomp <- ceiling(p / 10)
    if(missing(p.incl)) p.incl <- 1 / p
    if(missing(pmtm.budget)) pmtm.budget <- 0
    if(missing(varp.update)) varp.update <- 99
    
    sx <- apply(x.train, 2, function(z) diff(range(z)))
    x <- scale(x.train, center = FALSE, scale = sx)
    dimnames(x)[[2]] <- xnames
    if(missing(lpenalty)) lpenalty <- dt(scale(x), df = 10, log = TRUE)
    
    my <- mean(y.train)
    sy <- sd(y.train)
    obs <- as.numeric(scale(y.train))
    
    lamsq <- -log(prox) / 0.01
    rhosq <- Rsq / (1 - Rsq)
    lplam <- -lp.lam.f * sqrt(lamsq)
    lprho <- -lp.rho.f * log(ncomp) * log1p(rhosq)
    
    p.swap <- c(0, (1-p.ref) * dpois(1:p, 4)^(0.6))
    p.add <- c(1, (1-p.ref) * ppois(1:p, 1, lower = FALSE)^(0.3))
    p.ref <- c(0, rep(p.ref, p))
    p.rem <- 1 - p.add - p.swap - p.ref
    
    p.move <- rbind(p.add, p.rem, p.swap, p.ref)
    
    ## MAIN FUNCTION FOR FITTING THE ADDITIVE GP MODEL.
    ## --ARGUMENTS--
    ##   xvar     : the predictor matrix as a column major vector
    ##   yvar     : the response vector
    ##   pincl    : inclusion probability
    ##   dim      : integet vector of problem dimensions: n = #obs, p = pred, ncomp = #comp, max_rank = #max low rank, nlam = #lambda, nrho = #rho, nsweep = #mcmc iterations, nburn = #mcmc burn-in, budget = pMTM neighborhood budget
    ##   lamsqR   : grid of lam^2 values, length = nlam
    ##   rhosqR   : grid of rho^2 vales, length = nrho
    ##   hpar     : (a, b) of the gamma(a,b) prior on 1 / sigma^2
    ##   lplamR   : log-prior over lam^2 grid
    ##   lprhoR   : log-prior over rho^2 grid
    ##   pmove    : move probabilities, 4*(p + 1) vector of (p.add, p.remove, p.swap, p.refresh) for comp size = 0, 1, ..., p
    ##   lpenalty : log-penalty score for knots selection
    ##   tolchol  : tolerance levels for Cholesky factorizations
    ## --OUTPUTS--
    ##   nprop    : tally of different moves proposed
    ##   nacpt    : acceptance counts by move types
    ##   varp     : variable importance vector
    ##   ix_store : Markov chain sample of inclusion
    ##   par_store: Markov chain sample of covariance parameters
    ##   active_store : Flags for active components
    dyn.load("gp.so")
    if(pmtm.budget > 0){
        st <- system.time(
            out.c <- .C("addGP_pMTM", xvar = as.double(x), yvar = as.double(obs),
                        p.incl = as.double(p.incl), dim = as.integer(c(n, p, 
                                                                       ncomp, max.rank, length(lamsq), length(rhosq), nsweep, 
                                                                       nburn = burn*nsweep, pmtm.budget)), lamsqR = as.double(lamsq), 
                        rhosqR = as.double(rhosq), hpar = as.double(c(a,b)), 
                        lplamR = as.double(lplam), lprhoR = as.double(lprho), 
                        pmove = as.double(p.move), lpen = as.double(lpenalty), 
                        nprop = double(4), nacpt = double(4), varp = double(p), 
                        varp_update = as.integer(varp.update), 
                        ix.store = character(nsweep * ncomp), 
                        par.store = double(2*ncomp*nsweep), 
                        active.store = integer(nsweep * ncomp), 
                        tolchol = as.double(tol)))
    } else{
        st <- system.time(
            out.c <- .C("addGP", xvar = as.double(x), yvar = as.double(obs), 
                        p.incl = as.double(p.incl), dim = as.integer(c(n, p, 
                                                                       ncomp, max.rank, length(lamsq), length(rhosq), nsweep)),
                        lamsqR = as.double(lamsq), rhosqR = as.double(rhosq), 
                        hpar = as.double(c(a,b)), lplamR = as.double(lplam), 
                        lprhoR = as.double(lprho), pmove = as.double(p.move), 
                        lpen = as.double(lpenalty), nprop = double(4), 
                        nacpt = double(4), varp = double(p), 
                        ix.store = character(nsweep * ncomp), 
                        par.store = double(2*ncomp*nsweep), 
                        active.store = integer(nsweep * ncomp), 
                        tolchol = as.double(tol)))
    }
    dyn.unload("gp.so")
    out <- list(ix.store = t(matrix(out.c$ix.store, nrow = ncomp)), 
                pars.store = t(matrix(out.c$par.store, nrow = 2 * ncomp)),
                active.store = t(matrix(out.c$active.store, nrow = ncomp)),
                nprop = out.c$nprop, nacpt = out.c$nacpt, var.p = out.c$varp, 
                my = my, sy = sy, sx = sx, x = x, obs = obs, a = a, b = b, 
                lamsq.grid = lamsq, rhosq.grid = rhosq, runtime = st[3])
    class(out) <- "addGP"
    return(out)
    
}



predict.addGP <- function(out, x.test, burn = 0.1, nsamp = 200, fsamp = FALSE){
    
    ix.store <- out$ix.store
    pars.store <- out$pars.store
    
    nsweep <- nrow(pars.store)
    ncomp <- ncol(ix.store)
    
    incr <- (1 - burn) / nsamp
    ss <- unique(ceiling(nsweep * seq(burn + incr, 1, length = nsamp)))
    
    x <- out$x
    obs <- out$obs
    
    n <- nrow(x)
    p <- ncol(x)
    xnew <- scale(x.test, center = FALSE, scale = out$sx)
    nnew <- nrow(xnew)
    
    x.diff <- apply(x, 2, function(z) (rep(z, n) - rep(z, each = n))^2)
    xnew.diff <- sapply(1:p, function(j) (rep(xnew[,j], nnew) - rep(xnew[,j], each = nnew))^2)
    xcross.diff <- sapply(1:p, function(j) (rep(x[,j], nnew) - rep(xnew[,j], each = n))^2)
    
    a <- out$a
    b <- out$b
    
    get.f <- function(jj){
        K <- 0
        Knew <- 0
        Kcross <- 0
        
        for(k in 1:ncomp){
            par <- pars.store[jj,2*(k-1) + 1:2]
            ix <- uncollapse(ix.store[jj,k], collapse = "\\.", mode = "numeric")
            
            d2x <- matrix(rowSums(x.diff[,ix,drop=FALSE]), n, n)
            d2x.new <- matrix(rowSums(xnew.diff[,ix,drop=FALSE]), nnew, nnew)
            d2x.cross <- matrix(rowSums(xcross.diff[,ix,drop=FALSE]), n, nnew)
            K <- K + par[1] * exp(-par[2] * d2x)
            Knew <- Knew + par[1] * exp(-par[2] * d2x.new)
            Kcross <- Kcross + par[1] * exp(-par[2] * d2x.cross)
        }
        
        R <- chol(K + diag(1, n))
        z <- backsolve(R, obs, transpose = TRUE)
        sigSq <- 1/rgamma(1, (a + n)/2, (a * b + sum(z^2)) / 2)
        
        ZZ <- backsolve(R, Kcross, transpose = TRUE)
        fmean <- c(crossprod(ZZ, z))
        w <- fmean
        if(fsamp){
            fvar <- sigSq * (Knew - crossprod(ZZ)) 
            hvar <- chol(fvar, pivot = TRUE)
            pivot <- attr(hvar,"pivot")
            w <- rep(NA, nnew)
            w[pivot] <- fmean[pivot] + c(crossprod(hvar, rnorm(nnew)))
        }
        return(c(sigSq, w))
    }
    require(multicore)
    # f.samp <- out$my + out$sy * sapply(ss, get.f)
    temp <- matrix(unlist(mclapply(ss, get.f)), ncol = length(ss))
    f.samp <- out$my + out$sy * temp[-1,]
    sig2 <- out$sy^2 * temp[1,]
    return(list(f.samp=f.samp, sig2 = sig2))
}

summary.addGP <- function(out, burn = 0.1, nsamp = 200, ngraph = 10, use.names = FALSE){
    ix.store <- out$ix.store
    pars.store <- out$pars.store
    rhosq.grid <- out$rhosq.grid
    lamsq.grid <- out$lamsq.grid
    
    nsweep <- nrow(pars.store)
    ncomp <- ncol(ix.store)
    
    incr <- (1 - burn) / nsamp
    ss <- unique(ceiling(nsweep * seq(burn + incr, 1, length = nsamp)))
    
    plotpar <- par()
    par(mfrow = c(2,2), font.main = 1)
    
    boxplot(sqrt(pars.store[ss,2*(1:ncomp)]), outline = FALSE, col = 3, ylim = sqrt(range(lamsq.grid)), ann = FALSE, bty = "n")
    title(xlab = "Component", ylab = "Posterior spread", main = expression(lambda))
    abline(h = sqrt(lamsq.grid), lty = 1, col = 2)
    boxplot(sqrt(pars.store[ss,2*(1:ncomp) - 1]), outline = FALSE, col = 3, ann = FALSE, bty = "n", ylim = sqrt(range(rhosq.grid)))
    title(xlab = "Component", ylab = "Posterior spread", main = expression(rho))
    abline(h = sqrt(rhosq.grid), lty = 1, col = 2)
    
    
    p <- ncol(out$x)
    gp.cocc <- matrix(0,p,p)
    v.incl <- rep(0,p)
    for(jj in ss){
        jj.incl <- rep(FALSE, p)
        for(k in 1:ncomp){
            ixn <- uncollapse(ix.store[jj, k], collapse = "\\.", mode = "numeric")
            gp.cocc[ixn,ixn] <- gp.cocc[ixn,ixn] + 1
            jj.incl[ixn] <- TRUE
        }
        v.incl <- v.incl + jj.incl
    }
    v.incl <- v.incl / length(ss)
    if(use.names){
        dimnames(gp.cocc)[[1]] <- dimnames(gp.cocc)[[2]] <- dimnames(out$x)[[2]]
    } else {
        dimnames(gp.cocc)[[1]] <- dimnames(gp.cocc)[[2]] <- 1:p    
    }
    gp.cocc <- gp.cocc / length(ss)
    gp.cocc <- 20 *  gp.cocc / max(gp.cocc)
    plot(v.incl, ty = "h", bty = "n", xlab = "Predictor index", ylab = "Inclusion probability", main = "Predictor importance", ylim = c(0,1))
    
    
    #  iz <- which(diag(gp.cocc) > quantile(diag(gp.cocc), 1 - (ngraph - 0.5) / p))  
    iz <- order(diag(gp.cocc), decreasing = TRUE)[1:ngraph]
    require(igraph)
    graph <- graph.adjacency(gp.cocc[iz,iz,drop = FALSE], mode = "undirected", weighted = TRUE, diag = FALSE)
    plot(graph, vertex.size=diag(gp.cocc[iz,iz,drop = FALSE]), edge.width=E(graph)$weight, vertex.color = "orange", vertex.frame.color = "orange", vertex.label.cex = (diag(gp.cocc[iz,iz,drop = FALSE])/20)^.3)
    title(main = "Predictor interaction")
    
    suppressWarnings(par(plotpar))
    invisible(list(graph = gp.cocc, inclusion = v.incl))
}

## Utility functions

logsum <- function(lx) return(max(lx) + log(sum(exp(lx - max(lx)))))

expand.grid <- function(x, y){
    nx <- ifelse(is.null(dim(x)), length(x), nrow(x))
    ny <- ifelse(is.null(dim(y)), length(y), nrow(y))
    cbind(kronecker(rep(1, ny), x), kronecker(y, rep(1, nx)))
}

uncollapse <- function(str, collapse = "", mode = "character"){
    a <- unlist(strsplit(str, collapse))
    mode(a) <- mode
    return(a)
}

rmse <- function(a, b) return(sqrt(mean((a - b)^2)))

