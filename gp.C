#include "R.h"
#include "Rmath.h"
#include "R_ext/Applic.h"

#include "mydefs.h"
#include "spchol.h"


// GLOBAL VARIABLES
int n, p, max_rank, toprint, dopivot, nlam, nrho, ngrid;
double **x, *y, **lpenMat, priorVar, priorSD, sig_shape, sig_rate, *lsgrid, *lamsq, *rhosq, *lplam, *lprho, lpincl;

// VARIABLE REPEATEDLY USED BY AUX FUNCTIONS
int *rank, *pivot, *grank, *gpvt;
double **R, *d, *d2, *a, **F, **G, **H, *Rad, *gd, *yp, *z, *ryp, *pryp, *gptol, *cpar, *lpen, **S;


// VARIABLE SELECTION RADIAL BASIS COVARIANCE FUNCTION
double sevscov(int i, int j, double *cpar, int *include, double **K0, int addK){
	
	int l;
	double dist = 0.0, diff;
	for(l = 0; l < p; l++){
		if(include[l]){		
			diff = x[i][l] - x[j][l];
			dist += diff * diff;
		}
	}
	return cpar[0] * exp(-cpar[1] * dist) + addK * K0[i][j];
}	

// LOG LIKELIHOOD CALCULATION WITH A FIXED LAM & A VECTOR OF RHO
// --ARGUMENTS--
//   include : a logical vector of inclusion
//   K0      : covariance matrix contribution from other components
//   irho    : index for rho grid element
//   ilam    : index for lambda grid element
// --OUTPUT--
//   lsgrid  : a vector of log posterior scores for rho grid
void lsFn(int *include, double **K0, int irho, int ilam, double *lsgrid){
	
	int i, j, k, l;
	cpar[0] = rhosq[irho]; cpar[1] = lamsq[ilam];
	spchol(R, n, gptol[0], pivot, rank, max_rank, sevscov, cpar, include, K0, 1, n, d, dopivot, 1, 0, lpen, d2);

	for(i = 0; i < rank[0]; i++) a[pivot[i]] = 1.0;
	for(i = rank[0]; i < n; i++) a[pivot[i]] = 1.0 + d[i];
		
	for(i = 0; i < rank[0]; i++){				
		for(j = 0; j < n; j++)
			Rad[j] = R[i][j] / a[pivot[j]];
		triprod(R, rank[0], n, Rad, F[i], 0);
	}
	for(i = 0; i < rank[0]; i++) F[i][i] += 1.0;
		
	chol(G, rank[0], gptol[1], gpvt, grank, rank[0], gd, F, 0, 0);
	
	double logdetA = 0.0, logdetG = 0.0;
	
	for(i = 0; i < n; i++) logdetA += log(a[i]);
	for(i = 0; i < rank[0]; i++) logdetG += log(G[i][i]);
	
	for(i = 0; i < n; i++) yp[i] = y[pivot[i]] / a[pivot[i]];
	triprod(R, rank[0], n, yp, ryp, 0);
	
	for(i = 0; i < rank[0]; i++) pryp[gpvt[i]] = ryp[i];
	trisolve(G, rank[0], pryp, z, 1);  
	
	for(i = 0; i < n; i++) yp[i] = y[i] / sqrt(a[i]);
	
	lsgrid[0] = (- 0.5 * logdetA - logdetG 
				 - 0.5 * ((double)n + sig_shape) * log1p((sumsquares(yp,n) - sumsquares(z,rank[0])) / (sig_shape * sig_rate)) 
				 + lprho[irho] + lplam[ilam]);
}

// THE WORKS. DOES THREE THINGS:
// 1) LOG POSTERIOR CALCULATION WITH DICRETE INTEGRATION OVER BOTH LAMBDA AND RHO
// 2) GIBBS DRAW OF (LAM, RHO)
// 3) CALCULATES THE CORRESPONDING nxn COVARIANCE MATRIX
// --ARGUMENTS--
//   include   : a logical vector of inclusion
//   csize     : size of the component = sum(include)
//   K0        : covariance matrix contribution from other components
//   draw_cpar : a logical bit to flag update for (rhosq, lamsq) 
// --OUTPUT--
//   RETURN    : returns the log-posterior values
//   par       : Gibbs draw of (rhosq, lamsq)
//   K         : covariance matrix of the current component obtained with above par
double works(int *include, int csize, double **K0, int draw_cpar, double *par, double **K){
	int i, j, skip = 0;
	for(i = 0; i < n; i++){
		lpen[i] = 0.0;
		for(j = 0; j < p; j++) if(include[j]) lpen[i] += lpenMat[j][i];
	}
	for(i = 0; i < nlam; i++){
		for(j = 0; j < nrho; j++){
			lsFn(include, K0, j, i, lsgrid + skip);
			skip++;
		}
	}
	
	double ls = logsum(lsgrid, ngrid) + (double)csize * lpincl;

	if(draw_cpar){
	  int ii = rdraw(ngrid, lsgrid, 1);
	  par[0] = rhosq[ii % nrho]; par[1] = lamsq[ii / nrho];
	  for(i = 0; i < n; i++){
	    for(j = 0; j < i; j++) K[i][j] = K[j][i] = sevscov(i, j, par, include, K0, 0);
	    K[i][i] = par[0];
	  }
	}
	return ls;
}

// MAIN FUNCTION FOR FITTING THE ADDITIVE GP MODEL.
// --ARGUMENTS--
//   xvar     : the predictor matrix as a column major vector
//   yvar     : the response vector
//   pincl    : inclusion probability
//   dim      : integet vector of problem dimensions: n = #obs, p = #pred, ncomp = #comp, max_rank = #max low rank, nlam = #lambda, nrho = #rho
//   lamsqR   : grid of lam^2 values, length = nlam
//   rhosqR   : grid of rho^2 vales, length = nrho
//   hpar     : (a, b) of the gamma(a,b) prior on 1 / sigma^2
//   lplamR   : log-prior over lam^2 grid
//   lprhoR   : log-prior over rho^2 grid
//   pmove    : move probabilities, 4*(p + 1) vector of (p.add, p.remove, p.swap, p.refresh) for comp size = 0, 1, ..., p
//   lpenalty : log-penalty score for knots selection
//   tolchol  : tolerance levels for Cholesky factorizations
// --OUTPUTS--
//   nprop    : tally of different moves proposed
//   nacpt    : acceptance counts by move types
//   varp     : variable importance vector
//   ix_store : Markov chain sample of inclusion
//   par_store: Markov chain sample of covariance parameters
void addGP(double *xvar, double *yvar, double *pincl, int *dim, double *lamsqR, double *rhosqR, 
		   double *hpar, double *lplamR, double *lprhoR, double *pmove, double *lpenalty, 
		   double *nprop, double *nacpt, double *varp, char **ix_store, double *par_store, 
		   double *tolchol){
	
	n = dim[0];
	p = dim[1];
	int ncomp = dim[2];
	max_rank = dim[3];
	nlam = dim[4];
	nrho = dim[5];
	int nsweep = dim[6];
	toprint = 1; 
	dopivot = 1;
	
	ngrid = nlam * nrho;
	lamsq = lamsqR;
	rhosq = rhosqR;
	lplam = lplamR;
	lprho = lprhoR;
	
	sig_shape = hpar[0];
	sig_rate = hpar[1];
	lpincl = log(pincl[0]);
	
	int i, j, k, l;
	
	x = mymatrix(n, p);
	y = yvar;
	R = mymatrix(max_rank, n);
	d = vect(n);
	d2 = vect(n);
	a = vect(n);
	rank = ivect(1); rank[0] = max_rank;
	pivot = ivect(n);
	F = mymatrix(max_rank, max_rank);
	G = mymatrix(max_rank, max_rank);
	H = mymatrix(max_rank, max_rank);
	Rad = vect(n);
	gd = vect(max_rank);
	grank = ivect(1);
	gpvt = ivect(max_rank);
	yp = vect(n);
	z = vect(max_rank);
	ryp = vect(max_rank);
	pryp = vect(max_rank);
	gptol = tolchol;
	cpar = vect(2);
	lpenMat = mymatrix(p, n);
	lpen = vect(n);
	lsgrid = vect(ngrid);
	
	int pos = 0;
	for(j = 0; j < p; j++)
		for(i = 0; i < n; i++)
			x[i][j] = xvar[pos++];

	pos = 0;
	for(j = 0; j < p; j++)
		for(i = 0; i < n; i++)
			lpenMat[j][i] = lpenalty[pos++];
	
	GetRNGstate();

	Rprintf("n = %d, p = %d, ncomp = %d, nsweep = %d, max_rank = %d\n", n, p, ncomp, nsweep, max_rank); 
	
	double **K0 = mymatrix(n, n);
	for(i = 0; i < n; i++) for(j = 0; j < n; j++) K0[i][j] = 0.0;

	int **ix = (int **)R_alloc(ncomp, sizeof(int *)), *comp_size = ivect(ncomp);
	for(k = 0; k < ncomp; k++){
		ix[k] = ivect(p);
		for(j = 0; j < p; j++) ix[k][j] = 0;
		comp_size[k] = 0;
	}
	double ***K = (double ***)R_alloc(ncomp, sizeof(double **));
	for(k = 0; k < ncomp; k++) K[k] = mymatrix(n, n);
	
	double *ls = vect(ncomp);
	double *par = vect(2 * ncomp);
	for(k = 0; k < ncomp; k++){
	  ls[k] = works(ix[k], comp_size[k], K0, 1, par + 2*k, K[k]);

		if(ncomp > 1){
			for(i = 0; i < n; i++) for(j = 0; j < n; j++) K0[i][j] += K[k][i][j];			
		}
	}

	int sweep, move_type, nempty, ticker = nsweep / 10;
	nprop[0] = nprop[1] = nprop[2] = nprop[3] = 0.0;
	nacpt[0] = nacpt[1] = nacpt[2] = nacpt[3] = 0.0;
	double *varpc = vect(p), *varp_k = vect(p), *varpc_k = vect(p);
	for(j = 0; j < p; j++){
		varp[j] = 1.0;
		varpc[j] = 1.0;
	}
	int jj, jj2, par_pos = 0, ix_pos = 0, *ixnew = ivect(p);
	double u, log_alpha, pfwd, pbwd, lsnew, incr, *parnew = vect(2), **Knew = mymatrix(n, n);
		
	for(sweep = 0; sweep < nsweep; sweep++){
		for(k = 0; k < ncomp; k++){
			if(ncomp > 1){
				for(i = 0; i < n; i++) for(j = 0; j < n; j++) K0[i][j] -= K[k][i][j];
				ls[k] = works(ix[k], comp_size[k], K0, 1, par + 2*k, K[k]);
			}
			move_type = rdraw(4, pmove + 4 * comp_size[k], 0);

			nprop[move_type] += 1.0;
			
			switch(move_type){
				
				case 0 : // add
					for(j = 0; j < p; j++){
						varp_k[j] = ix[k][j] ? 0.0 : varp[j]; 
						ixnew[j] = ix[k][j];
					}
					jj = rdraw(p, varp_k, 0);
					ixnew[jj] = 1;
					lsnew = works(ixnew, comp_size[k] + 1, K0, 1, parnew, Knew);
					pfwd = pmove[4 * comp_size[k]] * varp_k[jj] / sum(varp_k, p);
					
					for(j = 0; j < p; j++) varpc_k[j] = ixnew[j] ? varpc[j] : 0.0;
					pbwd = pmove[4 * comp_size[k] + 4 + 1] * varpc_k[jj] / sum(varpc_k, p);
					
					log_alpha = lsnew - ls[k] + log(pbwd) - log(pfwd);
					u = runif(0.0, 1.0);
					if(log(u) < log_alpha){
						nacpt[0] += 1.0;
						ls[k] = lsnew;
						ix[k][jj] = 1;
						for(i = 0; i < n; i++) for(j = 0; j < n; j++) K[k][i][j] = Knew[i][j];
						par[2*k] = parnew[0]; par[2*k + 1] = parnew[1];
						comp_size[k]++;
					}
					break;
					
				case 1 : //remove
					for(j = 0; j < p; j++){
						varpc_k[j] = ix[k][j] ? varpc[j] : 0.0;
						ixnew[j] = ix[k][j];
					}
					jj = rdraw(p, varpc_k, 0);
					ixnew[jj] = 0;
					lsnew = works(ixnew, comp_size[k] - 1, K0, 1, parnew, Knew);
					
					pfwd = pmove[4 * comp_size[k] + 1] * varpc_k[jj] / sum(varpc_k, p);
					
					for(j = 0; j < p; j++) varp_k[j] = ixnew[j] ?  0.0 : varp[j];
					pbwd = pmove[4 * comp_size[k] - 4] * varp_k[jj] / sum(varp_k, p);
					
					log_alpha = lsnew - ls[k] + log(pbwd) - log(pfwd);
					u = runif(0.0, 1.0);
					if(log(u) < log_alpha){
						nacpt[1] += 1.0;
						ls[k] = lsnew;
						ix[k][jj] = 0;
						for(i = 0; i < n; i++) for(j = 0; j < n; j++) K[k][i][j] = Knew[i][j];
						par[2*k] = parnew[0]; par[2*k + 1] = parnew[1];
						comp_size[k]--;
					}
					break;
					
				case 2 : //swap
					for(j = 0; j < p; j++){
						varp_k[j] = ix[k][j] ? 0.0 : varp[j]; 
						varpc_k[j] = ix[k][j] ? varpc[j] : 0.0;
						ixnew[j] = ix[k][j];
					}
					jj = rdraw(p, varp_k, 0);
					jj2 = rdraw(p, varpc_k, 0);
					ixnew[jj] = 1; ixnew[jj2] = 0;
					lsnew = works(ixnew, comp_size[k], K0, 1, parnew, Knew);
					
					pfwd = (varp_k[jj] / sum(varp_k, p)) * (varpc_k[jj2] / sum(varpc_k, p));
					varp_k[jj] = 0.0; varp_k[jj2] = varp[jj2];
					varpc_k[jj] = varpc[jj]; varpc_k[jj2] = 0.0;
					pbwd = (varpc_k[jj] / sum(varpc_k, p)) * (varp_k[jj2] / sum(varp_k, p));
					
					log_alpha = lsnew - ls[k] + log(pbwd) - log(pfwd);
					u = runif(0.0, 1.0);
					if(log(u) < log_alpha){
						nacpt[2] += 1.0;
						ls[k] = lsnew;
						ix[k][jj] = 1; ix[k][jj2] = 0;
						for(i = 0; i < n; i++) for(j = 0; j < n; j++) K[k][i][j] = Knew[i][j];
						par[2*k] = parnew[0]; par[2*k + 1] = parnew[1];
					}
					break;
					
				default : //refresh
					nacpt[3] += 1.0;
					break;
			}
			if(ncomp > 1)
				for(i = 0; i < n; i++) for(j = 0; j < n; j++) K0[i][j] += K[k][i][j];

			incr = R_pow((double)(ncomp * (1+sweep)), -0.67);
			for(j = 0; j < p; j++){
				if(ix[k][j]){
					varp[j] += incr;
					varpc[j] = 1.0 / varp[j];
				}
			}
			
		}
		
		for(k = 0; k < ncomp; k++){
			par_store[par_pos++] = par[2*k];
			par_store[par_pos++] = par[2*k + 1];
			locator_string(ix[k], p, ix_store[ix_pos++]);

			if(par[2*k] == 0) 
			  for(j=0; j<p; j++) ix[k][j] = 0;
		}
		if((sweep + 1) % ticker == 0){
			Rprintf("sweep = %d. Active: ", sweep + 1);
			nempty = 0;
			for(k = 0; k < ncomp; k++){
				if(comp_size[k] == 0) {
					nempty++;	
				} else {
					Rprintf(" %s ", ix_store[ix_pos - ncomp + k]);
				}
			}
			Rprintf(" Empty: %d\n", nempty);
		}
	}
	PutRNGstate();
	
}

// UPDATE INCLUSION
// --ARGUMENTS--
//   accept    : 0 / 1 flag for a predictor inclusion proposal
//   nburn     : #burn-in mcmc iterations
//   sweep     : current mcmc iteration
//   nactive   : #active components; ncomp_min <= nactive <= ncomp
//   option    : update option: {0: basic, 1: similarity (depends on p x p matrix "S"); else do nothing}
// --OUTPUTS--
//   varp      : inclusion probability for predictor toggle
//   varpc     : inclusion probability for predictor toggle
void update_varp(int accept, int nburn, int sweep, int nactive, int option, int *ix, double *varp, double *varpc){
  if(1){
    int i, j, nsim;
    double csim, incr = R_pow(nactive, -0.67);
    incr *= sweep > nburn ? R_pow(sweep-nburn, -0.67) : ((double)sweep) / nburn;
    
    switch(option){
    case 0 : // basic
      for(i = 0; i < p; i++){
	varp[i] += incr * ix[i];
	varpc[i] = 1.0 / varp[i];
      }
      break;

    case 1 : // similarity	
      for(i = 0; i < p; i++){	  
	if(ix[i])
	  varp[i] += incr * ix[i];
	else{
	  nsim = 0; csim = 0.0;
	  for(j = 0; j < p; j++){
	    csim += S[i][j] * ix[j];
	    nsim += S[i][j] * ix[j] > 0 ? 1 : 0;
	  }
	  varp[i] += incr * (csim / nsim);
	  varpc[i] = 1.0 / varp[i];
	}
      }
      break;

    default : 
      break;
    }
  }
}

// WEIGHT FUNCTION
// --ARGUMENTS--
//   ix        : vector of predictor inclusion
//   varp      : vector of predictor propensity scores
//   j         : predictor index, (0, 2, ..., p-1)
//   move_type : 0 "add"; 1 "remove"; 2 "swap-add"
//   m         : neighborhood budget / size
// --OUTPUTS--
//   togg_pr   : inclusion probability for predictor toggle
double weightFn(int *ix, double *varp, int j, int move_type, int m){
  double togg_pr;
  int csize;
  
  switch(move_type){
  case 0 : // add
    togg_pr = ix[j]==0 ? m*varp[j] / (m*varp[j] + p) : 0.0;
    break;
    
  case 1 : // remove
    togg_pr = ix[j]==1 ? 1.0 : 0.0;
    break;

  case 2 : // swap (add)
    csize = isum(ix, p);
    togg_pr = ix[j]==0 ? (m*varp[j] / (m*varp[j] + p)) / csize : 0.0;
    break;

  default :
    break;
  }

  return togg_pr;
}

// MAIN FUNCTION FOR FITTING THE ADDITIVE GP MODEL.
// --ARGUMENTS--
//   xvar     : the predictor matrix as a column major vector
//   yvar     : the response vector
//   pincl    : inclusion probability
//   dim      : integet vector of problem dimensions: n = #obs, p = pred, ncomp = #comp, max_rank = #max low rank, nlam = #lambda, nrho = #rho, nsweep = #mcmc iterations, nburn = #mcmc burn-in, budget = pMTM neighborhood budget
//   lamsqR   : grid of lam^2 values, length = nlam
//   rhosqR   : grid of rho^2 vales, length = nrho
//   hpar     : (a, b) of the gamma(a,b) prior on 1 / sigma^2
//   lplamR   : log-prior over lam^2 grid
//   lprhoR   : log-prior over rho^2 grid
//   pmove    : move probabilities, 4*(p + 1) vector of (p.add, p.remove, p.swap, p.refresh) for comp size = 0, 1, ..., p
//   lpenalty : log-penalty score for knots selection
//   tolchol  : tolerance levels for Cholesky factorizations
// --OUTPUTS--
//   nprop    : tally of different moves proposed
//   nacpt    : acceptance counts by move types
//   varp     : variable importance vector
//   ix_store : Markov chain sample of inclusion
//   par_store: Markov chain sample of covariance parameters
//   active_store : Flags for active components
void addGP_pMTM(double *xvar, double *yvar, double *pincl, int *dim, double *lamsqR, double *rhosqR, 
		double *hpar, double *lplamR, double *lprhoR, double *pmove, double *lpenalty, 
		double *nprop, double *nacpt, double *varp, int *varp_update, char **ix_store, 
		double *par_store, int *active_store, double *tolchol){
	
  n = dim[0];
  p = dim[1];
  int ncomp = dim[2];
  max_rank = dim[3];
  nlam = dim[4];
  nrho = dim[5];
  int nsweep = dim[6];
  int nburn = dim[7];
  int budget = dim[8];
  int varp_option = varp_update[0];
  toprint = 1; 
  dopivot = 1;
	
  ngrid = nlam * nrho;
  lamsq = lamsqR;
  rhosq = rhosqR;
  lplam = lplamR;
  lprho = lprhoR;
	
  sig_shape = hpar[0];
  sig_rate = hpar[1];
  lpincl = log(pincl[0]) - log(1-pincl[0]);
	
  int i, j, k, l;
	
  x = mymatrix(n, p);
  y = yvar;
  R = mymatrix(max_rank, n);
  d = vect(n);
  d2 = vect(n);
  a = vect(n);
  rank = ivect(1); rank[0] = max_rank;
  pivot = ivect(n);
  F = mymatrix(max_rank, max_rank);
  G = mymatrix(max_rank, max_rank);
  H = mymatrix(max_rank, max_rank);
  Rad = vect(n);
  gd = vect(max_rank);
  grank = ivect(1);
  gpvt = ivect(max_rank);
  yp = vect(n);
  z = vect(max_rank);
  ryp = vect(max_rank);
  pryp = vect(max_rank);
  gptol = tolchol;
  cpar = vect(2);
  lpenMat = mymatrix(p, n);
  lpen = vect(n);
  lsgrid = vect(ngrid);  

  int pos = 0;
  for(j = 0; j < p; j++)
    for(i = 0; i < n; i++)
      x[i][j] = xvar[pos++];

  pos = 0;
  for(j = 0; j < p; j++)
    for(i = 0; i < n; i++)
      lpenMat[j][i] = lpenalty[pos++];

  int *active = ivect(ncomp);
  for(k = 0; k < ncomp; k++) active[k] = 1;
	
  GetRNGstate();

  Rprintf("n = %d, p = %d, ncomp = %d, nsweep = %d, nburn = %d, max_rank = %d\n", n, p, ncomp, nsweep, nburn, max_rank); 
  
  double **K0 = mymatrix(n, n);
  for(i = 0; i < n; i++) for(j = 0; j < n; j++) K0[i][j] = 0.0;
  
  int **ix = (int **)R_alloc(ncomp, sizeof(int *)), *comp_size = ivect(ncomp);
  for(k = 0; k < ncomp; k++){
    ix[k] = ivect(p);
    for(j = 0; j < p; j++) ix[k][j] = 0;
    comp_size[k] = 0;
  }
  double ***K = (double ***)R_alloc(ncomp, sizeof(double **));
  for(k = 0; k < ncomp; k++) K[k] = mymatrix(n, n);
  
  double *ls = vect(ncomp);
  double *par = vect(2 * ncomp);
  for(k = 0; k < ncomp; k++){
    ls[k] = works(ix[k], comp_size[k], K0, 1, par + 2*k, K[k]);
    
    if(ncomp > 1){
      for(i = 0; i < n; i++) for(j = 0; j < n; j++) K0[i][j] += K[k][i][j];			
    }
  }

  int sweep, move_type, nempty, ncomp_min = fmin(log(ncomp), ncomp), ticker = nsweep / 10;
  nprop[0] = nprop[1] = nprop[2] = nprop[3] = 0.0;
  nacpt[0] = nacpt[1] = nacpt[2] = nacpt[3] = 0.0;
  double *varpc = vect(p), *varp_k = vect(p), *varpc_k = vect(p);
  for(j = 0; j < p; j++){
    varp[j] = 1.0;
    varpc[j] = 1.0;
  }
  
  int jj, jj2, par_pos = 0, ix_pos = 0, ix2_pos = 0, *ixnew = ivect(p);
  double u, log_alpha, pfwd, pbwd, lsnew, incr, *parnew = vect(2), **Knew = mymatrix(n, n);

  int nactive = ncomp, cbudget = budget / nactive, togg_sum = 0, togg = 0, accept = 0, *ixnewc = ivect(p);
  double *lsfwd = vect(p), *lsbwd = vect(p), *lsswp = vect(p*p), *lsswpc = vect(p*p), *probs = vect(p*p), *probsc = vect(p*p);
  double a_pincl, b_pincl, p_update; 

  // default initialization for similarity matrix "S"
  S = mymatrix(p, p);
  for(i = 0; i < p; i++){
    for(j = 0; j < i; j++) S[i][j] = S[j][i] = 0.0;
    S[i][i] = 1.0;
  }
  
  for(sweep = 0; sweep < nsweep; sweep++){
    // update active components
    if(sweep+1 >= nburn && 1){      
      if(nactive < ncomp){
	p_update = sqrt(ncomp_min) / (sqrt(nactive) * (ncomp - nactive));
	nactive = 0;
	for(k = 0; k < ncomp; k++){	
	  if(nactive + ncomp-k-1 < ncomp_min)
	    active[k] = 1;
	  else
	    active[k] = runif(0.0, 1.0) < fmax(p_update, active[k]) ? 1 : 0;

	  nactive += active[k];
	}       	
      }      
      cbudget = budget / nactive;
    }
    
    // update pincl
    a_pincl = 1.0 + (double)isum(comp_size, ncomp);
    b_pincl = (double)(p-1 + p * nactive - isum(comp_size, ncomp));
    lpincl = rbeta(a_pincl, b_pincl);
    lpincl = log(lpincl) - log(1-lpincl);
  
    for(k = 0; k < ncomp; k++){
      if(active[k]){
	if(ncomp > 1){
	  for(i = 0; i < n; i++) for(j = 0; j < n; j++) K0[i][j] -= K[k][i][j];
	  ls[k] = works(ix[k], comp_size[k], K0, 1, par + 2*k, K[k]);
	}
	move_type = rdraw(4, pmove + 4 * comp_size[k], 0);

	nprop[move_type] += 1.0;

	for(j = 0; j < p; j++) ixnew[j] = ixnewc[j] = ix[k][j];

	accept = togg = togg_sum = 0;
	log_alpha = -INFINITY;

	switch(move_type){
	case 0 : // add
	  for(j = 0; j < p; j++){	  
	    togg = runif(0.0, 1.0) < weightFn(ix[k], varp, j, 0, cbudget) ? 1 : 0;
	    ixnew[j] = togg ? 1 : ix[k][j];
	    lsfwd[j] = togg ? works(ixnew, comp_size[k] + 1, K0, 0, parnew, Knew) : -INFINITY;

	    ixnew[j] = ix[k][j];
	    togg_sum += togg; 
	  }

	  if(togg_sum > 0){
	    jj = rdraw(p, lsfwd, 1);
	    ixnew[jj] = ixnewc[jj] = 1;
	    lsnew = lsfwd[jj];
	    pfwd = pmove[4 * comp_size[k]] * weightFn(ix[k], varp, jj, 0, cbudget);
	  
	    //rev
	    for(j = 0; j < p; j++){	    
	      togg = runif(0.0, 1.0) < weightFn(ixnew, varp, j, 1, cbudget) ? 1 : 0;
	      togg = (togg || j==jj);
	      ixnewc[j] = togg ? 0 : ixnew[j];
	      lsbwd[j] = togg ? works(ixnewc, comp_size[k], K0, 0, parnew, Knew) : -INFINITY;
	    
	      ixnewc[j] = ixnew[j];
	    }
	    pbwd = pmove[4 * comp_size[k] + 4 + 1] * weightFn(ixnew, varp, jj, 1, cbudget);
	    
	    log_alpha = lsnew - ls[k] + log(pbwd) - log(pfwd) + (lsbwd[jj]-logsum(lsbwd, p)) - (lsfwd[jj]-logsum(lsfwd, p));
	    if((sweep + 1) % ticker == 0) 
	      Rprintf("Component = %d, pMTM move = A/R (%d / %d), pred = %d, accept pr = %2.3f\n",k+1,togg_sum,comp_size[k]+1,jj+1,exp(log_alpha));
	  } else{
	    // do single A/R toggle
	    for(j = 0; j < p; j++){
	      varp_k[j] = ix[k][j] ? 0.0 : varp[j]; 
	      ixnew[j] = ix[k][j];
	    }
	    jj = rdraw(p, varp_k, 0);
	    ixnew[jj] = 1;
	    lsnew = works(ixnew, comp_size[k] + 1, K0, 1, parnew, Knew);
	    pfwd = pmove[4 * comp_size[k]] * varp_k[jj] / sum(varp_k, p);
					
	    for(j = 0; j < p; j++) varpc_k[j] = ixnew[j] ? varpc[j] : 0.0;
	    pbwd = pmove[4 * comp_size[k] + 4 + 1] * varpc_k[jj] / sum(varpc_k, p);
					
	    log_alpha = lsnew - ls[k] + log(pbwd) - log(pfwd);
	  }

	  u = runif(0.0, 1.0);
	  if(log(u) < log_alpha){
	    accept = 1;
	    nacpt[0] += 1.0;
	    ls[k] = lsnew;
	    ix[k][jj] = 1;
	    if(togg_sum > 0) works(ix[k], comp_size[k] + 1, K0, 1, parnew, Knew);

	    for(i = 0; i < n; i++) for(j = 0; j < n; j++) K[k][i][j] = Knew[i][j];
	    par[2*k] = parnew[0]; par[2*k + 1] = parnew[1];
	    comp_size[k]++;
	  }

	  break;
					
	case 1 : //remove
	  for(j = 0; j < p; j++){	  
	    togg = runif(0.0, 1.0) < weightFn(ix[k], varp, j, 1, cbudget) ? 1 : 0;
	    ixnew[j] = togg ? 0 : ix[k][j];
	    lsfwd[j] = togg ? works(ixnew, comp_size[k] - 1, K0, 0, parnew, Knew) : -INFINITY;
	  
	    ixnew[j] = ix[k][j];
	  }	
	
	  jj = rdraw(p, lsfwd, 1);
	  ixnew[jj] = ixnewc[jj] = 0;
	  lsnew = lsfwd[jj];
	  pfwd = pmove[4 * comp_size[k] + 1] * weightFn(ix[k], varp, jj, 1, cbudget);
	  
	  //rev
	  for(j = 0; j < p; j++){
	    togg = runif(0.0, 1.0) < weightFn(ixnew, varp, j, 0, cbudget) ? 1 : 0;
	    togg = (togg || j==jj);
	    ixnewc[j] = togg ? 1 : ixnew[j];
	    lsbwd[j] = togg ? works(ixnewc, comp_size[k], K0, 0, parnew, Knew) : -INFINITY;
	    	    
	    ixnewc[j] = ixnew[j];
	    togg_sum += togg;
	  }
	  pbwd = pmove[4 * comp_size[k] - 4] * weightFn(ixnew, varp, jj, 0, cbudget);
	
	  log_alpha = lsnew - ls[k] + log(pbwd) - log(pfwd) + (lsbwd[jj]-logsum(lsbwd, p)) - (lsfwd[jj]-logsum(lsfwd, p));	  
	  if((sweep + 1) % ticker == 0) 
	    Rprintf("Component = %d, pMTM move = R/A (%d / %d), pred = %d, accept pr = %2.3f\n",k+1,comp_size[k],togg_sum,jj+1,exp(log_alpha));

	  u = runif(0.0, 1.0);
	  if(log(u) < log_alpha){
	    accept = 1;
	    nacpt[1] += 1.0;
	    ls[k] = lsnew;
	    ix[k][jj] = 0;
	    works(ix[k], comp_size[k] - 1, K0, 1, parnew, Knew);
	    for(i = 0; i < n; i++) for(j = 0; j < n; j++) K[k][i][j] = Knew[i][j];
	    par[2*k] = parnew[0]; par[2*k + 1] = parnew[1];
	    comp_size[k]--;
	  }
	  break;
					
	case 2 : //swap
	  for(int r = 0; r < p; r++){
	    for(int a = 0; a < p; a++){
	      probs[r*p + a] = weightFn(ix[k], varp, r, 1, cbudget) * weightFn(ix[k], varp, a, 2, cbudget);
	      lsswp[r*p + a] = -INFINITY;
	      if(runif(0.0, 1.0) < probs[r*p + a]){
		ixnew[r] = 0; ixnew[a] = 1;
		lsswp[r*p + a] = works(ixnew, comp_size[k], K0, 0, parnew, Knew);

		ixnew[r] = 1; ixnew[a] = 0;
		togg_sum++;
	      }
	    }
	  }
	    
	  if(togg_sum > 0){
	    int swpdraw = rdraw(p*p, lsswp, 1);
	    jj = swpdraw / p;
	    jj2 = swpdraw - jj*p;
	  
	    ixnew[jj] = ixnewc[jj] = 0;
	    ixnew[jj2] = ixnewc[jj2] = 1;
	    lsnew = lsswp[swpdraw];
	    pfwd = probs[swpdraw];
	  
	    //rev
	    for(int r = 0; r < p; r++){
	      for(int a = 0; a < p; a++){
		probsc[r*p + a] = weightFn(ixnew, varp, r, 1, cbudget) * weightFn(ixnew, varp, a, 2, cbudget);
		lsswpc[r*p + a] = -INFINITY;
	      
		if(runif(0.0, 1.0) < probsc[r*p + a] || (a*p + r == swpdraw)){
		  ixnewc[r] = 0; ixnewc[a] = 1;
		  lsswpc[r*p + a] = works(ixnewc, comp_size[k], K0, 0, parnew, Knew);
		
		  ixnewc[r] = 1; ixnewc[a] = 0;
		  togg_sum++;
		}
	      }
	    }
	    pbwd = probsc[jj2*p + jj];
	  
	    log_alpha = lsnew - ls[k] + log(pbwd) - log(pfwd) + (lsswpc[jj2*p+jj]-logsum(lsswpc, p*p)) - (lsswp[swpdraw]-logsum(lsswp, p*p));
	    if((sweep + 1) % ticker == 0) 
	      Rprintf("Component = %d, pMTM move = S/S (%d), pred[R/A] = %d/%d, accept pr = %2.3f\n",k+1,2*comp_size[k]+togg_sum,jj+1,jj2+1,exp(log_alpha));
	  } else{
	    // do single S/S toggle
	    for(j = 0; j < p; j++){
	      varp_k[j] = ix[k][j] ? 0.0 : varp[j]; 
	      varpc_k[j] = ix[k][j] ? varpc[j] : 0.0;
	      ixnew[j] = ix[k][j];
	    }
	    jj = rdraw(p, varpc_k, 0);
	    jj2 = rdraw(p, varp_k, 0);
	    ixnew[jj] = 0; ixnew[jj2] = 1;
	    lsnew = works(ixnew, comp_size[k], K0, 1, parnew, Knew);
					
	    pfwd = (varp_k[jj] / sum(varp_k, p)) * (varpc_k[jj2] / sum(varpc_k, p));
	    varp_k[jj] = 0.0; varp_k[jj2] = varp[jj2];
	    varpc_k[jj] = varpc[jj]; varpc_k[jj2] = 0.0;
	    pbwd = (varpc_k[jj] / sum(varpc_k, p)) * (varp_k[jj2] / sum(varp_k, p));
					
	    log_alpha = lsnew - ls[k] + log(pbwd) - log(pfwd);
	  }	

	  u = runif(0.0, 1.0);
	  if(log(u) < log_alpha){
	    accept = 1;
	    nacpt[2] += 1.0;
	    ls[k] = lsnew;
	    ix[k][jj] = 0; ix[k][jj2] = 1;
	    if(togg_sum > 0) works(ix[k], comp_size[k], K0, 1, parnew, Knew);

	    for(i = 0; i < n; i++) for(j = 0; j < n; j++) K[k][i][j] = Knew[i][j];
	    par[2*k] = parnew[0]; par[2*k + 1] = parnew[1];
	  }
	
	  break;
					
	default : //refresh
	  nacpt[3] += 1.0;
	  break;
	}	
	if(ncomp > 1)
	  for(i = 0; i < n; i++) for(j = 0; j < n; j++) K0[i][j] += K[k][i][j];
	
	update_varp(accept, nburn, sweep+1, nactive, varp_option, ix[k], varp, varpc);
      }
    }
		
    for(k = 0; k < ncomp; k++){
      /*if(par[2*k] == 0.0){
	for(j = 0; j < p; j++) ix[k][j] = 0;
	comp_size[k] = 0;
	}*/
      if(sweep+1 >= nburn && 1) active[k] = (comp_size[k]==0 && par[2*k] <= rhosq[0]) ? 0 : 1;

      nactive = isum(active, ncomp);
      par_store[par_pos++] = par[2*k];
      par_store[par_pos++] = par[2*k + 1];
      active_store[ix2_pos++] = active[k];
      locator_string(ix[k], p, ix_store[ix_pos++]);
    }
    if((sweep + 1) % ticker == 0){
      Rprintf("\nsweep = %d. Active: %d. ", sweep + 1, nactive);
      Rprintveci("Active-flags: ", "%d", active, ncomp);
      Rprintf("Components: ");
      nempty = 0;
      for(k = 0; k < ncomp; k++){
	if(comp_size[k] == 0) {
	  nempty++;	
	} else {
	  Rprintf(" %s ", ix_store[ix_pos - ncomp + k]);
	}
      }      
      
      Rprintf(" Empty: %d\n----\n", nempty);
    }
  }
  PutRNGstate();
	
}


