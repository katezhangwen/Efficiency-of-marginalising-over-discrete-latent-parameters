# loading all dependencies
require(rjags)
require(coda)
require(here)
require(rstan)
require(ggplot2)
require(ggpubr)
library(ggmcmc)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
load.module("mix")
NT <- 2 # number of trials
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
mixing_p1 <- c(0.50, 0.50)
mu1 = c(5, -5)
prec1 = c(1/4, 1/4)
sigma1 = sqrt(1/prec1)

mixing_p2 <- c(0.50, 0.50)
mu2 = c(2.5, -2.5)
prec2 = c(1/4, 1/4)
sigma2 = sqrt(1/prec2)


mixing_p3 <- c(0.70, 0.30)
mu3 = c(5, -5)
prec3 = c(1/4, 1/4)
sigma3 = sqrt(1/prec3)


mixing_p4 <- c(0.70, 0.30)
mu4 = c(2.5, -2.5)
prec4 = c(1/4, 1/4)
sigma4 = sqrt(1/prec4)

N <- 200



library(LaplacesDemon)
data_generation <- function(datanum){
  I = numeric(N)
  x = numeric(N)
  I = rcat(n=N, eval(as.name(paste("mixing_p", datanum, sep=""))))
  for(i in 1:N){
    x[i] = rnorm(n=1, 
                 mean=eval(as.name(paste("mu", datanum, sep="")))[I[i]], 
                 sd=eval(as.name(paste("sigma", datanum, sep="")))[I[i]])
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
for (datanum in 1:4){

  for (t in 1:NT){
    data <- data_generation(datanum)
    y <- data
    N <- length(data)
    data_jags <- data_stan <- list(N=N, y=y)
    
    # for stan
    # record stan computation time
    time <- system.time(model_fit <-
                          suppressMessages(
                            stan(file = here("Models","two-componentModellnorm.stan"),
                                 data = data_stan, iter=iterations,
                                 chain=chains, warmup=burnin)))["elapsed"]
    result <- rbind(result, c("stan", NA, datanum, "Computation Time", t, time))
    
    # record stan model
    stan_mod_list[[t]] <- model_fit
    
    # record stan min ess
    ess <- min(summary(model_fit)$summary[,"n_eff"][0:6])
    result <- rbind(result, c("stan",NA,  datanum, "Min Effective Sample Size", t, ess))
    
    # record stan timeper min ess
    timeperess <-  time/ess
    result <- rbind(result, c("stan", NA,datanum,  "Time per min Effective Sample", t, timeperess))
    
    # record stan rhat
    temp <- summary(model_fit)
    rhat <- max(temp$summary[,"Rhat"][1:6])
    result <- rbind(result, c("stan",NA, datanum,  "Rhat", t, rhat))
    
    # unrestricted set of samplers for jags
    all_sampler()
    for(j in 1:3){
      model <- paste("two-componentModel",j,"lnorm.txt",sep="")
      
      # for jags-full 
      # computation time
      time <- system.time(model.fit <- 
                            jags.model(file =here("Models",model), 
                                       data=data_jags, 
                                       n.chains=chains))["elapsed"]
      + system.time(model.samples <- 
                      coda.samples(model.fit, parameters, 
                                   n.iter=iterations))["elapsed"]
      
      result <- rbind(result, c("jags", j, datanum,"Computation Time", t, time))
      
      # model
      jags_mod_list[[t]] <- model.samples
      
      # ess
      ess <- min(effectiveSize(model.samples)[0:6])
      result <- rbind(result, c("jags",j, datanum,"Min Effective Sample Size", t, ess))
      
      # time per ess
      timeperess <- time/ess
      result <- rbind(result, c("jags",j, datanum,"Time per min Effective Sample", t, timeperess))
      
      #rhat
      disc <- gelman.diag(model.samples, multivariate = FALSE)
      rhat <- max(disc$psrf[1:6,"Upper C.I."])
      result <- rbind(result, c("jags",j,datanum, "Rhat", t, rhat))
    }
    # restricted set of samplers for jags
    res_sampler()
    model <- paste("two-componentModel1lnorm.txt")
    # computation time
    time <- system.time(model.fit <- 
                          jags.model(file =here("Models",model), 
                                     data=data_jags, 
                                     n.chains=chains))["elapsed"]
    + system.time(model.samples <- 
                    coda.samples(model.fit, parameters, 
                                 n.iter=iterations))["elapsed"]
    
    result <- rbind(result, c("jags-full-restricted", NA,datanum, "Computation Time", t, time))
    
    # model
    jags_full_restricted_mod_list[[t]] <- model.samples
    
    # ess
    ess <- min(effectiveSize(model.samples)[0:6])
    result <- rbind(result, c("jags-full-restricted",NA, datanum,"Min Effective Sample Size", t, ess))
    
    # time per ess
    timeperess <- time/ess
    result <- rbind(result, c("jags-full-restricted",NA, datanum,"Time per min Effective Sample", t, timeperess))
    
    #rhat
    disc <- gelman.diag(model.samples, multivariate = FALSE)
    rhat <- max(disc$psrf[1:6,"Upper C.I."])
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

saveRDS(result, file = paste(here("Results"), "/two-component-result.rds", sep=""))
