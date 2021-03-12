// two-component mixture model

data {
  int<lower=1> N; // number of data points
  real y[N]; // observations
}

parameters {
  ordered[2] mu; // mean of mixture components, ordered for identifiability
  real<lower=0> sigma[2]; // standard deviation of mixture components
  simplex[2] mixing_p; // mixing proportions
}


model {
  vector[2] log_component;
  mu ~ normal(0, 2);
  sigma ~ normal(0, 2);
  mixing_p ~ dirichlet(rep_vector(1,2));
  for (n in 1:N){
    for (k in 1:2)
      log_component[k] = log(mixing_p[k]) + normal_lpdf(y[n] | mu[k], sigma[k]);
    target += log_sum_exp(log_component);
  }
}

generated quantities{
  vector[N] z; // discrete latent variable
  vector[2] log_component;
  for (n in 1:N){
    for (k in 1:2)
      log_component[k] = log(mixing_p[k]) + normal_lpdf(y[n] | mu[k], sigma[k]);
    z[n] = -exp(log_component[1] - log_sum_exp(log_component))+2;
  }

}



