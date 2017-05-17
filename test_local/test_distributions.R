library(exactmaxsel)
library(ggplot2)
library(microbenchmark)

## k = 2 ------------------------------------------------------------------

## TODO: Compare for other n0/n1/A -> other data, simulate? 
## Compare all methods for k = 2, small b
dat <- iris[1:100, ]
y <- factor(dat$Species)
x <- dat$Petal.Length
b <- seq(0, 10, 0.01)
betensky <- maxchisq:::pmaxchisq_betensky(b = b, k = nlevels(y))
miller <- maxchisq:::pmaxchisq_miller(b = b)
exact <- maxchisq:::pmaxchisq_exact(b = b, y = y, x = x)
perm <- maxchisq:::pmaxchisq_permutation(b = b, y = y, x = x)
##perm <- sapply(b, maxchisq:::pmaxchisq_permutation, y = y, x = x)
plot(b, betensky, col = "blue", ylim = c(0, 1))
points(b, miller, col = "red")
points(b, exact, col = "green")
points(b, perm, col = "cyan")

## TODO: Compare for other n0/n1/A -> other data, simulate? 
## Compare all methods for k = 2, large b
dat <- iris[1:100, ]
y <- factor(dat$Species)
x <- dat$Petal.Length
b <- seq(10, 20, 0.01)
betensky <- maxchisq:::pmaxchisq_betensky(b = b, k = nlevels(y))
miller <- maxchisq:::pmaxchisq_miller(b = b)
exact <- sapply(b, exactmaxsel::Ford, n0 = table(y)[1], n1 = table(y)[2], 
                A = table(x), statistic = "chi2")
perm <- maxchisq:::pmaxchisq_permutation(b = b, y = y, x = x)
plot(b, betensky, col = "blue", xlim = c(min(b), max(b)))
points(b, miller, col = "red")
points(b, 1-exact, col = "green")
points(b, perm, col = "cyan")

## k > 2 ------------------------------------------------------------------

## Compare Betensky and Permutation for different values of k
n <- 100
b <- seq(0, 20, 0.5)
k <- seq(2, 10)
permute <- c(TRUE, FALSE)
params <- expand.grid(b = b, k = k, permute = permute)
cdf <- apply(params, 1, function(p) {
  if (p["permute"]) {
    x <- rnorm(n)
    y <- factor(round(runif(n, 0.5, p["k"]+0.5)))
    levels(y) <- 1:p["k"]
    return(maxchisq:::pmaxchisq_permutation(b = p["b"], y = y, x = x))
  } else {
    return(maxchisq:::pmaxchisq_betensky(b = p["b"], k = p["k"]))
  }
})
dat <- cbind(params, cdf)
ggplot(dat, aes(x = b, y = cdf, colour = permute)) + 
  geom_point() + 
  facet_wrap(~k)

## Check approximations ------------------------------------------------------------------

maxchisq:::pmaxchisq_permutation(b = 5, y = iris$Species, x = iris$Petal.Length)

microbenchmark(NEW = maxchisq:::pmaxchisq_permutation(b = 5, y = iris$Species, x = iris$Petal.Length),
               times = 1)

b <- 0:50
k <- 3
maxchisq:::pmaxchisq_betensky(b = b, k = k)

b <- 0:50
maxchisq:::pmaxchisq_miller(b = b)

maxchisq:::pmaxchisq_permutation(b = b, y = iris$Species, x = iris$Petal.Length)

