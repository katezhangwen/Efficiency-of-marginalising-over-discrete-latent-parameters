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
set.seed(123)


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

# the quantities will be stored in separate dataframes
# with the following rows and columns
#             trial number
# stan
# jags-full
# jags-marg
time_df <- data.frame(row.names = c("stan", "jags-full", "jags-marg"))
timeperess_df <- data.frame(row.names = c("stan", "jags-full", "jags-marg"))
ess_df <- data.frame(row.names = c("stan", "jags-full", "jags-marg"))
rhat_df <- data.frame(row.names = c("stan", "jags-full", "jags-marg"))

# apart from the quantities of interest, create lists 
# for model objects
stan_mod_list <- list()
jags_full_mod_list <- list()
jags_marg_mod_list <- list()

for(t in 1:NT){
  # stan
  # record stan computation time
  time_df["stan",t] <- system.time(model.rater <- rater(anesthesia, dawid_skene()))["elapsed"]
  
  # record stan model
  model_fit <- get_stanfit(model.rater)
  stan_mod_list[[t]] <- model_fit
  
  # record stan min ess(84 is the number of continuous pars)
  ess_df["stan",t] <- min(summary(model_fit)$summary[,"n_eff"][0:84])
  
  # record stan timeper min ess
  timeperess_df["stan",t] <-  time_df["stan",t]/ess_df["stan",t]
  
  # record stan rhat
  rhat_df["stan",t] <- max(mcmc_diagnostics(model.rater)[,"Rhat"])
  
  # jags-full
  # record jags-full computation time
  time_df["jags-full",t] <- 
    system.time(model.fit <- jags.model(file =here("Models","Dawid-Skene.txt"), 
                                        data=data_jags, n.chains=chains))["elapsed"] 
  + system.time(model.samples <- coda.samples(model.fit, parameters, 
                                              n.iter=iterations))["elapsed"] 
  # record jags-full model
  jags_full_mod_list[[t]] <- model.samples
  
  # record jags-full min ess
  ess_df["jags-full",t] <- min(effectiveSize(model.samples)[0:84])
  
  # record jags-full timeper min ess
  timeperess_df["jags-full",t] <- time_df["jags-full",t]/ess_df["jags-full",t]
  
  # record jags-full rhat
  disc <- gelman.diag(model.samples, multivariate = FALSE)
  rhat_df["jags-full",t] <- max(disc$psrf[1:84,"Upper C.I."])
  
  
  # jags-marg
  # record jags-marg computation time
  time_df["jags-marg",t] <- 
    system.time(model.fit <- jags.model(file =here("Models","Dawid-Skene-marginalised.txt"), 
                                        data=data_jags, n.chains=chains))["elapsed"] 
  + system.time(model.samples <- coda.samples(model.fit, parameters, 
                                              n.iter=iterations))["elapsed"] 
  # record jags-marg model
  jags_marg_mod_list[[t]] <- model.samples
  
  # record jags-marg min ess
  ess_df["jags-marg",t] <- min(effectiveSize(model.samples)[0:84])
  
  # record jags-marg timeper min ess
  timeperess_df["jags-marg",t] <- time_df["jags-marg",t]/ess_df["jags-marg",t]
  
  # record jags-marg rhat
  disc <- gelman.diag(model.samples, multivariate = FALSE)
  rhat_df["jags-marg",t] <- max(disc$psrf[1:84,"Upper C.I."])
  
}

# save results
saveRDS(time_df, file = paste(here("Results"), "Dawid-Skene-time.rds"))
saveRDS(rhat_df, file = paste(here("Results"), "Dawid-Skene-rhat.rds"))
saveRDS(timeperess_df, file = paste(here("Results"), "Dawid-Skene-timeperess.rds"))
saveRDS(ess_df, file = paste(here("Results"), "Dawid-Skene-ess.rds"))

    
        
        
        