set.seed(1056915)
library(LaplacesDemon)

# this data set will generate node inconsistency error, not used
# generate data from a two-component normal mixture model
N = 80
mixing_p = c(0.35, 0.65)
mu = c(34, 62)
prec = c(1/81, 1/16)
sigma = sqrt(1/prec)

I = numeric(N)
x = numeric(N)

I = rcat(n=N, mixing_p)
for(i in 1:N){
  x[i] = rnorm(n=1, mean=mu[I[i]], sd=sigma[I[i]])
}

x
par(mar = c(5,5,3,1))
hist(x, breaks = 15)
plot(density(x))


write.table(x, file = "data1.txt", 
            row.names = FALSE, col.names = FALSE)

