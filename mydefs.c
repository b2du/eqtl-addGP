#include "R.h"
#include "Rmath.h"
#include "R_ext/Applic.h"
#include "mydefs.h"

double * vect(int n){
	return (double *)R_alloc(n, sizeof(double));
}

int * ivect(int n){
	return (int *)R_alloc(n, sizeof(int));
}

double ** mymatrix(int nr, int nc){
	int   i;
	double **m;	
	m = (double **) R_alloc(nr, sizeof(double *));
	for (i = 0; i < nr; i++)
		m[i] = (double *) R_alloc(nc, sizeof(double));
	return m;
}


void Rprintvec(char *a, char *format, double *x, int n){
  int i;
	Rprintf("%s", a);
  for(i = 0; i < n; i++)
    Rprintf(format, x[i]);
  Rprintf("\n");
}

void Rprintmat(char *a, char *format, double **x, int m, int n, int flip){
	int i, j;
	Rprintf("%s\n", a);
  for(i = 0; i < m; i++){
		for(j = 0; j < n; j++)
			Rprintf(format, x[i][j]);
		Rprintf("\n");
	}
}


void Rprintveci(char *a, char *format, int *x, int n){
	int i;
	Rprintf("%s", a);
  for(i = 0; i < n; i++)
    Rprintf(format, x[i]);
  Rprintf("\n");
}


double sumsquares(double *x, int n){
	double ss = 0.0;
	int i;
	for(i = 0; i < n; i++)
		ss += x[i] * x[i];
	return ss;
}

double inprod(double *x, double *y, int n){
	double ip = 0.0;
	int i;
	for(i = 0; i < n; i++)
		ip += x[i] * y[i];
	return ip;
}


double rnormtrunc(double mu, double sigma, double lo, double hi){
	double u = runif(0.0, 1.0);
	double p = u * pnorm(hi, mu, sigma, 1, 0) + (1.0 - u) * pnorm(lo, mu, sigma, 1, 0);
	if(p <= 0.0) p = 1.0e-10;
	if(p >= 1.0) p = 1.0 - 1.0e-10;
	return qnorm(p, mu, sigma, 1, 0);
}



double matern(double x, double phi, int kappa){ 
	/*Returns the Matern function for x, phi and kappa.*/ 
	
	/* Variables */ 
	double ans, cte; 
	double uphi=x/phi; 
		
	/* Matern */ 
	
	if (uphi==0) return 1; 
	else{ 
		if (kappa==0.5) 
			ans = exp(-uphi); 
		else { 
			cte = R_pow(2, (-(kappa-1)))/gammafn(kappa); 
			ans = cte * R_pow(uphi, kappa) * bessel_k(uphi,kappa,1); 
		} 
	} 
	
	return ans; 
} 

double vmax(double *x, int n){
	int i;
	double xmax = x[0];
	for(i = 1; i < n; i++) if(x[i] > xmax) xmax = x[i];
	return xmax;
}

double logsum(double *lx, int n){
	double lxmax = vmax(lx, n), a = 0.0;
	int i;
	for(i = 0; i < n; i++) a += exp(lx[i] - lxmax);
	return lxmax + log(a);
}

double sum(double *x, int n){
	double a = 0.0;
	int i;
	for(i = 0; i < n; i++) a += x[i];
	return a;
}

int isum(int *x, int n){
        int a = 0;
	for(int i = 0; i < n; i++) a += x[i];
	return a;
}

int rdraw(int n, double *prob, int inlog){
	double psum, u = runif(0.0, 1.0), cprob;
	int j = 0;

	if(inlog){
		psum = logsum(prob, n);	
		cprob = exp(prob[0] - psum);
		while(u > cprob && j < n - 1){
			j++;
			cprob += exp(prob[j] - psum);
		}
	} else {
		psum = sum(prob, n);
		cprob = prob[0] / psum;
		while(u > cprob && j < n - 1){
			j++;
			if(prob[j] > 0.0) cprob += prob[j] / psum;
		}
	}
	return j;
}


void locator_string(int *ix, int n, char *a){
	const char *fmt[2];	fmt[0] = "%d";	fmt[1] = ".%d";	
	int i, skip = 0;
	for(i = 0; i < n; i++){
		if(ix[i]){
			sprintf(a + skip, fmt[skip > 0], i + 1);
			skip = strlen(a);
		}
	}
}
