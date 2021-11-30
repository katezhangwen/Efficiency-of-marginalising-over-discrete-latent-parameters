# loading all dependencies
require(rjags)
require(coda)
require(here)
require(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
load.module("mix")
NT <- 5 # number of trials
set.seed(123)

# some setup
iterations <- 3000
burnin <- floor(iterations/2)
chains <- 3
parameters = c("mu", "sigma", "mixing_p","z")

# functions to restrict the set of samplers that jags can use
all_sampler <- function(){
  set.factory("mix::LDA", TRUE, type="sampler")
  set.factory("mix::DirichletCat", TRUE, type="sampler")
  set.factory("bugs::BinomSlice", TRUE, type="sampler")
  set.factory("bugs::RW1", TRUE, type="sampler")
  set.factory("bugs::Censored", TRUE, type="sampler")
  set.factory("bugs::Sum", TRUE, type="sampler")
  set.factory("bugs::DSum", TRUE, type="sampler")
  set.factory("bugs::Conjugate", TRUE, type="sampler")
  set.factory("bugs::Dirichlet", TRUE, type="sampler")
  set.factory("zbugs::MNormal", TRUE, type="sampler")
  set.factory("base::Finite", TRUE, type="sampler")
  set.factory("base::Slice", TRUE, type="sampler")
  set.factory("mix::TemperedMix", TRUE, type="sampler")
}

res_sampler <- function(){
  set.factory("mix::LDA", FALSE, type="sampler")
  set.factory("mix::DirichletCat", FALSE, type="sampler")
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
  set.factory("mix::TemperedMix", FALSE, type="sampler")
}

# parameters and functions for generating data
mixing_equal <- c(1/3, 1/3, 1/3)
mixing_diff <- c(0.5, 0.3, 0.2)
mixing_pair <- c(0.5, 0.3, 0.2)
prec = c(1/4, 1/4, 1/4)
sigma = sqrt(1/prec)

mu_equal_big <- c(-10.5, 0, 10.5)
mu_equal_small <- c(-6, 0 , 6)

mu_diff_big <- c(-9, -3, 12 )
mu_diff_small <- c(-4.5, 0 , 6.5)


mixing_p1 <- mixing_p2 <- mixing_p3 <- mixing_p4 <- mixing_equal

mixing_p5 <- mixing_p6 <- mixing_p7 <- mixing_p8 <- mixing_diff

mu1 <- mu5 <- mu_equal_big
mu2 <- mu6 <- mu_equal_small
mu3 <- mu7 <- mu_diff_big
mu4 <- mu8 <- mu_diff_small



N <- 200



library(LaplacesDemon)
data_generation <- function(datanum){
  I = numeric(N)
  x = numeric(N)
  I = rcat(n=N, eval(as.name(paste("mixing_p", datanum, sep=""))))
  for(i in 1:N){
    x[i] = rnorm(n=1, 
                 mean=eval(as.name(paste("mu", datanum, sep="")))[I[i]], 
                 sd=sigma)
  }
  return(x)
}



# running jags and stan for data 1,2,3,4 and record the quantities of
# interest 
# quantities of interest:
# *time/min ess
# * min ess
# * time
# * Rhat

# the quantities will be stored in a dataframe with the following columns
# * model
# * model type
# * data
# * quantity
# * trial 
# * value
result <- data.frame()
stan_mod_list <- list()
jags_mod_list <- list()
jags_full_restricted_mod_list <- list()

# start running jags and stan
for (datanum in 1:8){
  
  for (t in 1:NT){
    data <- data_generation(datanum)
    y <- data
    N <- length(data)
    data_jags <- data_stan <- list(N=N, y=y)
    
    # for stan
    # record stan computation time
    time <- system.time(model_fit <-
                          suppressMessages(
                            stan(file = here("Models","three-componentModellnorm.stan"),
                                 data = data_stan, iter=iterations,
                                 chain=chains, warmup=burnin)))["elapsed"]
    result <- rbind(result, c("stan", NA, datanum, "Computation Time", t, time))
    
    # record stan model
    stan_mod_list[[t]] <- model_fit
    
    # record stan min ess
    ess <- min(summary(model_fit)$summary[,"n_eff"][0:9])
    result <- rbind(result, c("stan",NA,  datanum, "Min Effective Sample Size", t, ess))
    
    # record stan timeper min ess
    timeperess <-  time/ess
    result <- rbind(result, c("stan", NA,datanum,  "Time per min Effective Sample", t, timeperess))
    
    # record stan rhat
    temp <- summary(model_fit)
    rhat <- max(temp$summary[,"Rhat"][1:9])
    result <- rbind(result, c("stan",NA, datanum,  "Rhat", t, rhat))
    
    # unrestricted set of samplers for jags
    all_sampler()
    for(j in 1:3){
      model <- paste("three-componentModel",j,"lnorm.txt",sep="")
      
      # for jags-full 
      # computation time
      t1 <- system.time(model.fit <- 
                            jags.model(file =here("Models",model), 
                                       data=data_jags, 
                                       n.chains=chains))["elapsed"]
      t2 <- system.time(model.samples <- 
                      coda.samples(model.fit, parameters, 
                                   n.iter=iterations))["elapsed"]
      time <- t1 + t2
      
      result <- rbind(result, c("jags", j, datanum,"Computation Time", t, time))
      
      # model
      jags_mod_list[[t]] <- model.samples
      
      # ess
      ess <- min(effectiveSize(model.samples)[0:9])
      result <- rbind(result, c("jags",j, datanum,"Min Effective Sample Size", t, ess))
      
      # time per ess
      timeperess <- time/ess
      result <- rbind(result, c("jags",j, datanum,"Time per min Effective Sample", t, timeperess))
      
      #rhat
      disc <- gelman.diag(model.samples, multivariate = FALSE)
      rhat <- max(disc$psrf[1:9,"Upper C.I."])
      result <- rbind(result, c("jags",j,datanum, "Rhat", t, rhat))
    }
    # restricted set of samplers for jags
    res_sampler()
    model <- paste("three-componentModel1lnorm.txt")
    # computation time
    t1 <- system.time(model.fit <- 
                          jags.model(file =here("Models",model), 
                                     data=data_jags, 
                                     n.chains=chains))["elapsed"]
    t2 <-  system.time(model.samples <- 
                    coda.samples(model.fit, parameters, 
                                 n.iter=iterations))["elapsed"]
    time <- t1 + t2
    result <- rbind(result, c("jags-full-restricted", NA,datanum, "Computation Time", t, time))
    
    # model
    jags_full_restricted_mod_list[[t]] <- model.samples
    
    # ess
    ess <- min(effectiveSize(model.samples)[0:9])
    result <- rbind(result, c("jags-full-restricted",NA, datanum,"Min Effective Sample Size", t, ess))
    
    # time per ess
    timeperess <- time/ess
    result <- rbind(result, c("jags-full-restricted",NA, datanum,"Time per min Effective Sample", t, timeperess))
    
    #rhat
    disc <- gelman.diag(model.samples, multivariate = FALSE)
    rhat <- max(disc$psrf[1:9,"Upper C.I."])
    result <- rbind(result, c("jags-full-restricted",NA, datanum,"Rhat", t, rhat))
    
  }
}
# cleaning up model name
colnames(result) <- c("model","type","data","quantity","trial", "value")
mod_name <- result$model
mod_name[result$model == "jags" & result$type == 1] <- "jags-full"
mod_name[result$model == "jags" & result$type == 2] <- "jags-marg-inbuilt"
mod_name[result$model == "jags" & result$type == 3] <- "jags-marg"
result$model <- mod_name
result <- subset(result, select = -c(type))

saveRDS(result, file = paste(here("Results"), "/three-component-result.rds", sep=""))
