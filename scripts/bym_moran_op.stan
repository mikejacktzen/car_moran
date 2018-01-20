data {
int<lower=0> N; // number of areas
int<lower=0> p; // number of covariates
int<lower=0> q; // number of eigenvalues in dim red
matrix[q,q] Q_s;  // q by q matrix for sparse precision random effect
matrix[N,q] M; // matrix of q<<n eigenvectors
vector[N] E; // expected cases
matrix[N,p] X; // covariate matrix

int<lower=0> Z[N]; // observed cases
}
 
transformed data {
matrix[q,q] Sig; // inverse to be cholesky decomposed
matrix[q,q] L_sre; //lower cholesky factor for sre
 
// enforce symmetricity to feed into cholesky
 
// inverse to be cholesky decomposed
// note: Q_s sparse replaces precision matrix of of vanilla ICAR
Sig = (inverse(Q_s));
 
for (i in 1:q)
for (j in 1:i)
Sig[i,j] = Sig[j,i];
 
// L cholesky decompose of spatial re cov
L_sre = cholesky_decompose(Sig);
}
 
parameters {
//////////////////////////////////
// Non-varying Fixed Effects
//////////////////////////////////

real beta0; // intercept
vector[p] Beta; // vector of covariate params
 
//////////////////////////////////
// spatial random effects
////////////////////////////////// 

vector[q] gamma_z; // N(0,1) to be affinely transformed into CAR sre
 
//////////////////////////////////
// spatial effect precision prior param
////////////////////////////////// 

real<lower=1e-5> prec_sre; // precision of CAR
}
 
transformed parameters {

//////////////////////////////////////////////////////
// Turn precision hyper param to sd hyper param
// MAT transform unit normal to car 
//////////////////////////////////////////////////////

real tau_sq;
vector[q] gamma; // spatial re after affine trans
 
// turn precisions to variances
tau_sq = 1.0 / (prec_sre * prec_sre);
 
// car sre after affine transform
gamma = sqrt(tau_sq)*(L_sre*gamma_z);
}
 
model {

//////////////////////////////////
//overall poisson risk by N
//////////////////////////////////

vector[N] log_mu;  //overall poisson risk

//////////////////////////////////
// spatial hyper prior 
//////////////////////////////////

prec_sre ~ gamma(1,1); // precision of spatial re

////////////////////////////////// 
// spatial re unit normal (before MAT affine transform)
//////////////////////////////////

gamma_z ~ normal(0, 1);


//////////////////////////////////
// Nonvarying fixed effect priors
//////////////////////////////////

// beta0: no explicit prior / implicit flat
Beta ~ normal(0.0, sqrt(1e5)); // each Beta_p has exchangeable normal priors
 
 
//////////////////////////////////
// likelihood
//////////////////////////////////

// Recall:
// Z is N by 1 vector
// log_mu is N by 1 vector
// E is N by 1 vector 
// beta0 is scalar but implies N by 1 vector of same scalars
// Beta is p by 1 vector
// gamma is q by 1 vector
// M is N by q matrix
 
// process model
log_mu = log(E) + 
M*gamma +
beta0 + 
X*Beta
;

// data model
Z ~ poisson_log(log_mu);
}

