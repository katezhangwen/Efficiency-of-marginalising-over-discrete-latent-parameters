// three-component mixture model

data {
  int<lower=1> N; // number of data points
  real y[N]; // observations
}

parameters {
  ordered[3] mu; // mean of mixture components, ordered for identifiability
  real<lower=0> sigma[3]; // standard deviation of mixture components
  simplex[3] mixing_p; // mixing proportions
}


model {
  vector[3] log_component;
  mu ~ normal(0, 7);
  sigma ~ lognormal(0, 1);
  mixing_p ~ dirichlet(rep_vector(1,3));
  for (n in 1:N){
    for (k in 1:3)
      log_component[k] = log(mixing_p[k]) + normal_lpdf(y[n] | mu[k], sigma[k]);
    target += log_sum_exp(log_component);
  }
}




generated quantities{
  vector[N] z; // discrete latent variable
  vector[3] log_component;
  for (n in 1:N){
    for (k in 1:3)
      log_component[k] = log(mixing_p[k]) + normal_lpdf(y[n] | mu[k], sigma[k]);
    z[n] = -2*exp(log_component[1] - log_sum_exp(log_component))
            -exp(log_component[2] - log_sum_exp(log_component))+3;
  }

}


