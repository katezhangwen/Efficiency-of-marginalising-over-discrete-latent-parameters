library(rstan)
require(rjags)
require(coda)
require(here)

load(file = "data/anesthesia.rda")

# data
N <- nrow(anesthesia) # total number of annotations
K <- 4 # number of annotation categories
I <- 45 # number of items
J <- 5 # number of annotators
ii <- anesthesia$item # item index for annotation n
jj <- anesthesia$rater # annotator for annotation n
y <- anesthesia$rating # annotation for observation n


# parameters

alpha <- rep(3,K)
n <- 8 # pseudocount of hypothetical prior observations
p <- 0.6 # approximate probability of a correct rating
beta <- array((n*(1-p)/(K-1)),c(J,K,K))
beta
for(j in 1:J){
  for (k in 1:K){
    beta[j,k,k] <- n*p
  }
}


# setting
iterations <- 3000
burnin <- floor(iterations/2)
chains <- 3

data_stan <- list(N=N, y=y, K=K, I=I, J=J, ii=ii, jj=jj,beta=beta, alpha=alpha)
data_jags <- list(N=N, y=y, K=K, I=I, J=J, ii=ii, jj=jj)

# with stan
model_fit <- suppressMessages(stan(file = here("Models","DawidSkene.stan"), data = data_stan
                  ,iter=iterations, chain=chains, warmup=burnin))
summary(model_fit)

# with jags (discrete model)
parameters = c("pi", "theta","z")
model.fit <- jags.model(file = here("Models","Dawid-Skene.txt"), data=data_jags, n.chains=chains)
model.samples <- coda.samples(model.fit, parameters,n.iter=iterations)
summary(window(model.samples, start = burnin))

# with jags (marginalised model)
parameters = c("pi", "theta","z")
model.fit <- jags.model(file = here("Models","Dawid-Skene-marginalised.txt"), data=data_jags, n.chains=chains)
model.samples <- coda.samples(model.fit, parameters,n.iter=iterations)
summary(window(model.samples, start = burnin))

