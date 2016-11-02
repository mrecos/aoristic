// Sample from Gaussian process
// All data parameters must be passed as a list to the Stan call
// Based on original file from https://code.google.com/p/stan/source/browse/src/models/misc/gaussian-process/
  
  data {
    int<lower=1> N;
    real x[N];
    real eta_sq;
    real rho_sq;
    real sigma_sq;
    real a;
    real p;
    real l;
    real pi;
    real sigma_b;
    real sigma_v;
    real offset_1;
  }
transformed data {
  vector[N] mu;
  cov_matrix[N] Sigma;
  for (i in 1:N) 
    mu[i] <- 0;
  for (i in 1:N) 
    for (j in 1:N)
    
  //// RBF, aka squared exp kernel, aka Gaussian KERNEL
     Sigma[i,j] <- eta_sq * exp(-rho_sq*pow(x[i] - x[j],2)) +
     if_else(i==j, sigma_sq, 0.0);
        
    // Sigma[i,j] <- eta_sq * exp(-(pow((x[i] - x[j]),2)/(1/rho_sq))) + 
    // if_else(i==j, sigma_sq, 0.0);
     
  //// PERIODIC KERNEL
    // Sigma[i,j] <- eta_sq * exp(-((2*pow(sin(pi*fabs((x[i] - x[j]))/p),2))/pow(l,2))) +           if_else(i==j, sigma_sq, 0.0);
    
  //// LINEAR KERNEL
    // Sigma[i,j] <- (sigma_b + sigma_v * ((x[i] - offset_1) * (x[j] - offset_1))) +
    // if_else(i==j, sigma_sq, 0.0);
     
  //// RBF, aka squared exp * periodic kernel
    // Sigma[i,j] <- eta_sq * exp(-rho_sq*pow(x[i] - x[j],2)) *
    //     (eta_sq * exp(-((2*pow(sin(pi*fabs((x[i] - x[j]))/p),2))/pow(l,2)))) +
    //     if_else(i==j, sigma_sq, 0.0);
         
  //// linear * periodic kernel
     // Sigma[i,j] <-  (sigma_b + sigma_v * ((x[i] - offset_1) * (x[j] - offset_1))) *
     // (eta_sq * exp(-((2*pow(sin(pi*fabs((x[i] - x[j]))/p),2))/pow(l,2)))) +
     // if_else(i==j, sigma_sq, 0.0);
     
  //// linear + periodic kernel
     // Sigma[i,j] <-  (sigma_b + sigma_v * ((x[i] - offset_1) * (x[j] - offset_1))) +
     // (eta_sq * exp(-((2*pow(sin(pi*fabs((x[i] - x[j]))/p),2))/pow(l,2)))) +
     // if_else(i==j, sigma_sq, 0.0);

  //// linear * linear kernel
     // Sigma[i,j] <- (sigma_b + sigma_v * ((x[i] - offset_1) * (x[j] - offset_1))) *
     // (sigma_b + sigma_v * ((x[i] - offset_1) * (x[j] - offset_1))) +
     // if_else(i==j, sigma_sq, 0.0);
     
}
parameters {
  vector[N] y;
}
model {
  y ~ multi_normal(mu,Sigma);
}

