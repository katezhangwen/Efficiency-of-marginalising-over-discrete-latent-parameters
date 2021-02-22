set.seed(201229)
library(LaplacesDemon)

# generate data from a two-component normal mixture model
N = 800
mixing_p = c(0.60, 0.40)
mu = c(0, -1)
prec = c(9, 4)
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



write.table(x, file = "data3.txt", 
            row.names = FALSE, col.names = FALSE)

