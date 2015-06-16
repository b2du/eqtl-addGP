# addGP
Bayesian Additive GP regression

1. gp.c gives an implementation written completely in C.

2. gp.R gives an implementation written mostly in C, wrapped by R.  The C code draws helper functions from mydefs.c and spchol.c.  Additional description of input arguments to addGP are provided within.  To compile you must run:
```
R CMD SHLIB gp.c mydefs.c spchol.c
```

3. The C implementation allows for an optional low-rank approximation.  The default max-rank value is set at n = #training samples.  Consider reducing this for large training data.  A potential rule of thumb is (2 * log(n))^2.

4. GP is self-regulatory, i.e., it imposes some penalty on large rho.  The discrete grid choice of lambda also imposes regularization.  Additional penalty could be imposed by setting lp.lam.f and/or lp.rho.f to positive values.  Not currently used in simulated or real-data experiments.

5. simu_gp.R runs addGP for simulated regression examples.  Specific options for pMTM include (a) pmtm.budget: neighborhood budget across "ncomp" components.  When pmtm.budget = 0, the implementation defaults to "Toggle" moves.  A reasonable default for the maximum number of components is sqrt(p), with a small per-component pMTM budget (e.g., 5); and (b) varp.update = {0, 99} (default: 99, adaptation is disabled; '0' will update using a basic diminishing adaptation scheme, see [arXiv-preprint](http://arxiv.org/abs/1411.7009) for details; additional "similarity" based updating may be implemented).


Note: Variaible selection for addGP can be performed by either "Toggle" and "pMTM".  These correspond to specific stochastic search variable selection schemes via Metropolis Hastings.  Both allow for adaptive MCMC via the updating of predictor propensity scores (i.e., "varp"); pMTM is a multiple-try analogue of Toggle, wherein multiple (i) add-remove, (ii) remove-add, and (iii) swap-swap paired moves are considered at each MCMC iteration.
