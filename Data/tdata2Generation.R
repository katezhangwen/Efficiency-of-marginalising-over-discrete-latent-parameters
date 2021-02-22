set.seed(201229)
library(LaplacesDemon)

# generate data from a three-component normal mixture model
N = 200
mixing_p = c(0.40, 0.30, 0.3)
mu = c(7, -5, 1)
prec = c(1, 1/3, 1/4)
sigma = sqrt(1/prec)

I = numeric(N)
x = numeric(N)

I = rcat(n=N, mixing_p)
I
for(i in 1:N){
  x[i] = rnorm(n=1, mean=mu[I[i]], sd=sigma[I[i]])
}

x
par(mar = c(5,5,3,1))
hist(x, breaks = 15)
plot(density(x))


write.table(x, file = "tdata2.txt", 
            row.names = FALSE, col.names = FALSE)



