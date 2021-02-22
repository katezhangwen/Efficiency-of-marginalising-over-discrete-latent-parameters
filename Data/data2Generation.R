set.seed(201229)
library(LaplacesDemon)

# generate data from a two-component normal mixture model
N = 350
mixing_p = c(0.20, 0.80)
mu = c(10, -2)
prec = c(1/4, 1/25)
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


write.table(x, file = "data2.txt", 
            row.names = FALSE, col.names = FALSE)

