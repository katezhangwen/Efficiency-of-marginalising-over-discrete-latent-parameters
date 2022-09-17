## Efficiency-of-marginalising-over-discrete-latent-parameters

For Bayesian model involving discrete variables, it is possible to sample a value of the discrete parameter at each iteration using methods such as Gibbs Sampling. With algorithms that do not support sampling discrete parameters, an alternative is to marginalise out the discrete parameters from the likelihood. This repo contains scripts which explore which approach is more computationally efficient. 

Preprint of our analysis is available on arXiv: https://arxiv.org/abs/2204.06313 

## About this repo
- **Models**: All stan models and marginalised and non-marginalised versions of JAGS models
- **MCMC**: R scripts which simulate data and measure the performance of each model
- **Results**: Visualisation of results
