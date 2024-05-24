data {
  int<lower=1> N;   // Number of individuals
  int<lower=1> J;   // Number of Betas: aka n weeks + 1 is_weekend
  matrix[N, J] dum; // matrix of indicator values
  int Y[N];         // observed reporting delays
}

parameters {
  vector[J] betas;         // one param for each week + 1 is_weekend
  real<lower=0.01> phi;    // a single dispersion param
}

transformed parameters {
  
  vector[N] mu;          // each person has their own mu
  mu = exp(dum * betas); // dot-product, gives mu_vector
  
}

model {
  
  // prior for beta and size
  betas ~ normal(0, 1);
  phi ~ normal(0, 1);
  
  // likelihood
  Y ~ neg_binomial_2(mu, phi);
  
}

generated quantities {
 vector[N] y_rep;
 for(n in 1:N){
  y_rep[n] = neg_binomial_2_rng(mu[n], phi); //posterior draws to get posterior predictive checks
 }
}

