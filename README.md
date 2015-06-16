# addGP
Bayesian Additive GP regression

GP is self-regulatory, i.e., it imposes some penalty on large rho.  The discrete grid choice of lambda also imposes regularization.  Additional penalty could be imposed by setting lp.lam.f and/or lp.rho.f to positive values (not used in existing simulated or real-data experiments).

Note: Variaible selection for addGP can be performed by either "Toggle" and "pMTM".  These correspond to specific stochastic search variable selection schemes via Metropolis Hastings.  Both allow for adaptive MCMC via the updating of predictor propensity scores (i.e., "varp"); pMTM is a multiple-try analogue of Toggle, wherein multiple (i) add-remove, (ii) remove-add, and (iii) swap-swap paired moves are considered at each MCMC iteration.


1. gp.c gives an implementation written completely in C.  This implementation allows for an optional low-rank approximation.  The default "max-rank" value is set at n = #training samples.  Consider reducing this for larger training data.  A potential rule of thumb is (2 * log(n))^2.

2. gp.R gives an implementation written mostly in C, wrapped by R.  The C code draws helper functions from mydefs.c and spchol.c.  Additional description of input arguments to addGP are provided within.  To compile you must run:
        ```
        R CMD SHLIB gp.c mydefs.c spchol.c
        ```

3. simu_gp.R runs addGP for simulated regression examples.  Specific options for pMTM include (a) pmtm.budget: neighborhood budget across "ncomp" components.  When pmtm.budget = 0, the implementation defaults to "Toggle" moves.  A reasonable default for the maximum number of components is sqrt(p), with a small per-component pMTM budget (e.g., 5); and (b) varp.update = {0, 99} (default: 99, adaptation is disabled; '0' will update using a basic diminishing adaptation scheme, see [arXiv-preprint](http://arxiv.org/abs/1411.7009) for details; additional "similarity" based updating may be implemented).


Input / Ouput summary for addGP                                                                                                          
====
1. ARGUMENTS
    * xvar - the predictor matrix as a column major vectors
    * yvar - the response vector
    * pincl - inclusion probability
    * dim - integet vector of problem dimensions:
        * n = # training samples
        2. p = # predictors
        3. ncomp = # compomponents
        4. max_rank = # max low rank
        5. nlam = length of discrete lambda grid
        6. nrho = length of discrete rho grid
        7. nsweep = #MCMC iterations
        8. nburn = #mcmc burn-in
        9. budget = pMTM neighborhood budget
    * lamsqR - grid of lam^2 values, length = nlam
    * rhosqR - grid of rho^2 vales, length = nrho
    * hpar - (a, b) of the gamma(a,b) prior on 1 / sigma^2
    * lplamR - log-prior over lam^2 grid
    * lprhoR - log-prior over rho^2 grid                                                                                                              
    * pmove -  move probabilities, 4*(p + 1) vector of (p.add, p.remove, p.swap, p.refresh) for comp size = 0, 1, ..., p
    * lpenalty - log-penalty score for knots selection
    * tolchol - tolerance levels for Cholesky factorizations                                                                                                 

2. OUTPUTS
    * nprop - tally of different moves proposed                                                                                                            
    * nacpt - acceptance counts by move types
    * varp - variable importance vector
    * ix_store - Markov chain sample of inclusion
    * par_store - Markov chain sample of covariance parameters
    * active_store - Flags for active components
    