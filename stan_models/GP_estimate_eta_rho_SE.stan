// Predict from Gaussian Process
// estimate sigma_sq and rho_sq
// All data parameters must be passed as a list to the Stan call
// Based on original file from https://code.google.com/p/stan/source/browse/src/models/misc/gaussian-process/

data {
int<lower=1> N1;
vector[N1] x1;
vector[N1] y1;
int<lower=1> N2;
vector[N2] x2;
real sigma_sq;

}

transformed data {
int<lower=1> N;
vector[N1+N2] x;
vector[N1+N2] mu;
N <- N1 + N2;
for (n in 1:N1) x[n] <- x1[n];
for (n in 1:N2) x[N1 + n] <- x2[n];
for (i in 1:N) mu[i] <- 0;
}
parameters {
  vector[N2] y2;
  real<lower=0> rho_sq;
  real<lower=0> eta_sq;
}
transformed parameters {
}
model {
vector[N] y;
matrix[N, N] Sigma;
  // off-diagonal elements
  for (i in 1:(N-1)){
    for (j in (i+1):N){
  
  //// Squared Exponential (RBF) Kerenl
    Sigma[i,j] <- eta_sq * exp(-rho_sq * pow(x[i] - x[j],2));
    Sigma[j,i] <- Sigma[i,j];
    }
  }
  // diagonal elements
  for (k in 1:N){
    Sigma[k,k] <- 1 + sigma_sq; // + jitter
  }
  rho_sq ~ cauchy(0,5);
  eta_sq ~ cauchy(0,5);
  for (n in 1:N1) y[n] <- y1[n];
  for (n in 1:N2) y[N1 + n] <- y2[n];
  y ~ multi_normal(mu,Sigma);
}
