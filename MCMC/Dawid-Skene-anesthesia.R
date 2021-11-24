# loading all dependencies
library("rater")
require(rjags)
require(coda)
require(here)
require(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
load(file = here("data","anesthesia.rda"))
NT <- 2 # number of trials
set.seed(1)


#rater(anesthesia, dawid_skene())

# using the rater package to run stan
data("anesthesia", package = "rater")

# some setup
N <- nrow(anesthesia) # total number of annotations
K <- 4 # number of annotation categories
I <- 45 # number of items
J <- 5 # number of annotators
ii <- anesthesia$item # item index for annotation n
jj <- anesthesia$rater # annotator for annotation n
y <- anesthesia$rating # annotation for observation n

data_jags <- data_stan <- list(N=N, y=y, K=K, I=I, J=J, ii=ii, jj=jj)
iterations <- 6000
burnin <- floor(iterations/2)
chains <- 3
parameters = c("pi", "theta","z")


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


for(t in 1:NT){
  # stan
  # record stan computation time
  time <- system.time(model.rater <- rater(anesthesia, dawid_skene()))["elapsed"]
  result <- rbind(result, c("stan", "Computation Time", t, time))
  
  # record stan model
  model_fit <- get_stanfit(model.rater)
  stan_mod_list[[t]] <- model_fit
  
  # record stan min ess(84 is the number of continuous pars)
  ess <- min(summary(model_fit)$summary[,"n_eff"][0:84])
  result <- rbind(result, c("stan", "Min Effective Sample Size", t, ess))
  
  # record stan timeper min ess
  timeperess <-  time/ess
  result <- rbind(result, c("stan", "Time per min Effective Sample", t, timeperess))
  
  # record stan rhat
  rhat <- max(mcmc_diagnostics(model.rater)[,"Rhat"])
  result <- rbind(result, c("stan", "Rhat", t, rhat))
  
  # jags-full
  # record jags-full computation time
  time <- 
    system.time(model.fit <- jags.model(file =here("Models","Dawid-Skene.txt"), 
                                        data=data_jags, n.chains=chains))["elapsed"] 
  + system.time(model.samples <- coda.samples(model.fit, parameters, 
                                              n.iter=iterations))["elapsed"] 
  result <- rbind(result, c("jags-full", "Computation Time", t, time))
  
  # record jags-full model
  jags_full_mod_list[[t]] <- model.samples
  
  # record jags-full min ess
  ess <- min(effectiveSize(model.samples)[0:84])
  result <- rbind(result, c("jags-full", "Min Effective Sample Size", t, ess))
  
  # record jags-full timeper min ess
  timeperess <-  time/ess
  result <- rbind(result, c("jags-full", "Time per min Effective Sample", t, timeperess))
  
  # record jags-full rhat
  disc <- gelman.diag(model.samples, multivariate = FALSE)
  rhat <- max(disc$psrf[1:84,"Upper C.I."])
  result <- rbind(result, c("jags-full", "Rhat", t, rhat))
  
  # jags-marg
  # record jags-marg computation time
  time <- 
    system.time(model.fit <- jags.model(file =here("Models","Dawid-Skene-marginalised.txt"), 
                                        data=data_jags, n.chains=chains))["elapsed"] 
  + system.time(model.samples <- coda.samples(model.fit, parameters, 
                                              n.iter=iterations))["elapsed"] 
  result <- rbind(result, c("jags-marg", "Computation Time", t, time))
  
  # record jags-marg model
  jags_marg_mod_list[[t]] <- model.samples
  
  # record jags-marg min ess
  ess <- min(effectiveSize(model.samples)[0:84])
  result <- rbind(result, c("jags-marg", "Min Effective Sample Size", t, ess))
  
  # record jags-marg timeper min ess
  timeperess <-  time/ess
  result <- rbind(result, c("jags-marg", "Time per min Effective Sample", t, timeperess))
  
  # record jags-marg rhat
  disc <- gelman.diag(model.samples, multivariate = FALSE)
  rhat <- max(disc$psrf[1:84,"Upper C.I."])
  result <- rbind(result, c("jags-marg", "Rhat", t, rhat))
}

colnames(result) <- c("model", "quantity", "trial", "value")
# save results
saveRDS(result, file = paste(here("Results"), "/Dawid-Skene-result.rds", sep=""))


    
        
        
        