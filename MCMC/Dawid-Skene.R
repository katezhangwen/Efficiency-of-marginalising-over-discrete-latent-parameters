# loading all dependencies
library("rater")
require(rjags)
require(coda)
require(here)
require(rstan)
require(posterior)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
NT <- 5 # number of trials
set.seed(123)


# run jags and stan and record the quantities of interest

# quantities of interest:
# *time/min ess
# * min ess
# * time
# * Rhat

# the quantities will be stored in a dataframe with the following columns
# * model
# * quantity
# * trial 
# * value
result <- data.frame()

# apart from the quantities of interest, create lists 
# for model objects
stan_mod_list <- list()
jags_full_mod_list <- list()
jags_marg_mod_list <- list()
jags_full_restricted_mod_list <- list()


# functions to restrict the set of samplers that jags can use
all_sampler <- function(){
  set.factory("bugs::BinomSlice", TRUE, type="sampler")
  set.factory("bugs::RW1", TRUE, type="sampler")
  set.factory("bugs::Censored", TRUE, type="sampler")
  set.factory("bugs::Sum", TRUE, type="sampler")
  set.factory("bugs::DSum", TRUE, type="sampler")
  set.factory("bugs::Conjugate", TRUE, type="sampler")
  set.factory("bugs::Dirichlet", TRUE, type="sampler")
  set.factory("bugs::MNormal", TRUE, type="sampler")
  set.factory("base::Finite", TRUE, type="sampler")
  set.factory("base::Slice", TRUE, type="sampler")
}

res_sampler <- function(){
  set.factory("bugs::BinomSlice", FALSE, type="sampler")
  set.factory("bugs::RW1", FALSE, type="sampler")
  set.factory("bugs::Censored", FALSE, type="sampler")
  set.factory("bugs::Sum", FALSE, type="sampler")
  set.factory("bugs::DSum", FALSE, type="sampler")
  set.factory("bugs::Conjugate", FALSE, type="sampler")
  set.factory("bugs::Dirichlet", TRUE, type="sampler")
  set.factory("bugs::MNormal", FALSE, type="sampler")
  set.factory("base::Finite", FALSE, type="sampler")
  set.factory("base::Slice", TRUE, type="sampler")
}


for(t in 1:NT){
  # simulating data
  #' Simulate incomplete data from the Dawid-Skene model.
  #'
  #' @param I The number of ratings.
  #' @param J The number of raters.
  #' @param K The number of categories.
  #' @param n_ratings The total number of ratings in the simulated data.
  #' @param pi The pi parameter
  #' @param theta The theta parameter
  #'
  #' @return A list with two components:
  #'   * z: The latent class of each item;
  #'   * data: The simulated data in a long format.
  #'
  simulate_ds_incomplete <- function(I, J, K,
                                     n_ratings,
                                     pi, theta) {
    
    # We need there to be at least one rating per item.
    stopifnot(n_ratings > I)
    
    n_extra_ratings <- n_ratings - I
    
    z <- sample(1:K, size = I, prob = pi, replace = TRUE)
    
    p <- rep(1 / I, I)
    n_extra_ratings_per_item <- rmultinom(1, n_extra_ratings, p)
    
    # Add in the minimal number of ratings which are distributed perfectly
    # uniformly.
    n_ratings_per_item <- n_extra_ratings_per_item + 1
    
    # NB: This is in ragged list form
    rater_indexes <- lapply(n_ratings_per_item, function(x) {
      # We currently sample from the raters uniformly - this could be extended.
      sample(1:J, size = x, replace = TRUE)
    })
    
    # Use rater_indexes as a placeholder.
    item_indexes <- rep(1:I, lengths(rater_indexes))
    out <- cbind(item_indexes, unlist(rater_indexes), 0)
    for (i in 1:nrow(out)) {
      out[i, 3] <- sample(1:K, 1, prob = theta[out[i, 2], z[[out[i, 1]]], ])
    }
    colnames(out) <- c("item", "rater", "rating")
    
    list(
      z = z,
      data = out
    )
  }
  
  #' Build the theta parameter of the Dawid-Skene model
  #'
  #' @param on_diag The on diagonal entries of each matrix
  #' @param J The number of raters (number of matrices)
  #' @param K The number of categories (the dimensions of each matrix)
  #'
  #' @return A c(J, K, K) array; the theta parameter
  #'
  build_theta <- function(on_diag, J, K) {
    
    stopifnot(length(on_diag) == 1)
    stopifnot(is.numeric(on_diag))
    
    theta <- array(0, dim = c(J, K, K))
    for (j in 1:J) {
      theta[j, ,] <- array((1 - on_diag) / (K - 1))
      for (k in 1:K) {
        theta[j, k, k] <- on_diag
      }
    }
    
    theta
  }
  
  
  I <- 100
  K <- 5
  J <- 5
  n_ratings <- N <- 500
  pi <- rep(1 / K, K)
  on_diag <- 0.7
  theta <- build_theta(0.7, J, K)
  
  sim <- simulate_ds_incomplete(I, J, K, n_ratings, pi, theta)
  
  
  # setup for models

  ii <- sim$data[,1] # item index for annotation n
  jj <- sim$data[,2] # annotator for annotation n
  y <- sim$data[,3] # annotation for observation n
  alpha <- rep(3,K)
  beta <- matrix(ncol = K, nrow = K)
  p <- 0.6
  for (i in 1:K){
    for (j in 1:K){
      if(i ==j){
        beta[i,j] <- N*p
      }
      else {
        beta[i,j] <- N*(1-p)/(K-1)
      }
    }
  }
  
  beta_slice <- beta
  beta <- array(dim = c(J, K, K))
  for (j in 1:J) {
    beta[j, , ] <- beta_slice
  }
  
  data_jags <-list(N=N, y=y, K=K, I=I, J=J, ii=ii, jj=jj)
  data_stan <- list(N=N, y=y, K=K, I=I, J=J, ii=ii, jj=jj, beta=beta, alpha=alpha)
  iterations <- 6000
  burnin <- floor(iterations/2)
  chains <- 3
  parameters = c("pi", "theta","z")
  
  pi_init <- rep(1/K, K)
  theta_init <- array(0.2 / (K - 1), c(J, K, K))
  for (j in 1:J) {
    diag(theta_init[j, ,]) <- 0.8
  }
  
  
  # stan
  # record stan computation time
  time <- system.time(model_fit <-
                        suppressMessages(
                          stan(file = here("Models","DawidSkene.stan"),
                               data = data_stan, iter=iterations,
                               init = function(n) list(theta = theta_init, pi = pi_init),
                               chain=chains, warmup=burnin)))["elapsed"]
  
  #time <- system.time(model.rater <- rater(sim$data, "dawid_skene"))["elapsed"]
  result <- rbind(result, c("stan", "Computation Time", t, time))
  
  # record stan model
  stan_mod_list[[t]] <- model_fit
  stan_summary <- summarise_draws(as_draws_array(model_fit))
  
  #model_fit <- get_stanfit(model.rater)
  #stan_mod_list[[t]] <- model_fit
  
  # record stan min ess(130 is the number of continuous pars)
  ess <- min(stan_summary$ess_bulk[0:130])
  result <- rbind(result, c("stan", "Min Effective Sample Size", t, ess))
  
  # record stan timeper min ess
  timeperess <-  time/ess
  result <- rbind(result, c("stan", "Time per min Effective Sample", t, timeperess))
  
  # record stan rhat
  rhat <- max(stan_summary$rhat[0:130])
  result <- rbind(result, c("stan", "Rhat", t, rhat))
  
  # unrestricted set of samplers for jags
  all_sampler()
  # jags-full
  # record jags-full computation time
  
  t1 <- system.time(model.fit <- jags.model(file =here("Models","Dawid-Skene.txt"), 
                                            data=data_jags, n.chains=chains))["elapsed"] 
  t2 <-  system.time(model.samples <- coda.samples(model.fit, parameters, 
                                                   n.iter=iterations))["elapsed"]
  time <- t1 + t2
  result <- rbind(result, c("jags-full", "Computation Time", t, time))
  
  # record jags-full model
  jags_full_mod_list[[t]] <- model.samples
  jags_summary <- summarise_draws(as_draws_array(model.samples))
  
  # record jags-full min ess
  ess <- min(jags_summary$ess_bulk[0:130])
  result <- rbind(result, c("jags-full", "Min Effective Sample Size", t, ess))
  
  # record jags-full timeper min ess
  timeperess <-  time/ess
  result <- rbind(result, c("jags-full", "Time per min Effective Sample", t, timeperess))
  
  # record jags-full rhat
  rhat <- max(jags_summary$rhat[0:130])
  result <- rbind(result, c("jags-full", "Rhat", t, rhat))
  
  # jags-marg
  # record jags-marg computation time
  t1 <- system.time(model.fit <- jags.model(file =here("Models","Dawid-Skene-marginalised.txt"), 
                                            data=data_jags, n.chains=chains))["elapsed"] 
  t2 <- system.time(model.samples <- coda.samples(model.fit, parameters, 
                                                  n.iter=iterations))["elapsed"] 
  time <- t1 + t2
  result <- rbind(result, c("jags-marg", "Computation Time", t, time))
  
  # record jags-marg model
  jags_marg_mod_list[[t]] <- model.samples
  jags_summary <- summarise_draws(as_draws_array(model.samples))
  
  # record jags-marg min ess
  ess <- min(jags_summary$ess_bulk[0:130])
  result <- rbind(result, c("jags-marg", "Min Effective Sample Size", t, ess))
  
  # record jags-marg timeper min ess
  timeperess <-  time/ess
  result <- rbind(result, c("jags-marg", "Time per min Effective Sample", t, timeperess))
  
  # record jags-marg rhat
  rhat <- max(jags_summary$rhat[0:130])
  result <- rbind(result, c("jags-marg", "Rhat", t, rhat))
  
  # restructed set of samplers for jags
  res_sampler()
  
  # jags-full
  # record jags-full computation time
  t1 <- system.time(model.fit <- jags.model(file =here("Models","Dawid-Skene.txt"), 
                                            data=data_jags, n.chains=chains))["elapsed"] 
  t2 <-  system.time(model.samples <- coda.samples(model.fit, parameters, 
                                                   n.iter=iterations))["elapsed"]
  time <- t1 + t2
  result <- rbind(result, c("jags-full-restricted", "Computation Time", t, time))
  
  # record jags-full model
  jags_full_restricted_mod_list[[t]] <- model.samples
  jags_summary <- summarise_draws(as_draws_array(model.samples))
  
  # record jags-full min ess
  ess <- min(jags_summary$ess_bulk[0:130])
  result <- rbind(result, c("jags-full-restricted", "Min Effective Sample Size", t, ess))
  
  # record jags-full timeper min ess
  timeperess <-  time/ess
  result <- rbind(result, c("jags-full-restricted", "Time per min Effective Sample", t, timeperess))
  
  # record jags-full rhat
  rhat <- max(jags_summary$rhat[0:130])
  result <- rbind(result, c("jags-full-restricted", "Rhat", t, rhat))
}

colnames(result) <- c("model", "quantity", "trial", "value")
# save results
saveRDS(result, file = paste(here("Results"), "/Dawid-Skene-result-new-ess-bulk.rds", sep=""))







