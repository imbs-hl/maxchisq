library(maxchisq)
context("distributions")

n <- 50
x <- rnorm(n)
y <- factor(round(rbinom(n, 1, 0.5)))
b <- runif(20, 0, 50)

test_that("Betensky pvalue is between 0 and 1", {
  p <- maxchisq:::pmaxchisq_betensky(b = b, k = nlevels(y))
  expect_true(all(p <= 1 & p >= 0))
})

test_that("Miller pvalue is between 0 and 1", {
  p <- maxchisq:::pmaxchisq_miller(b = b)
  expect_true(all(p <= 1 & p >= 0))
})

test_that("Permutation pvalue is between 0 and 1", {
  p <- maxchisq:::pmaxchisq_permutation(b = b, y = y, x = x)
  expect_true(all(p <= 1 & p >= 0))
})

test_that("Exact pvalue is between 0 and 1", {
  p <- maxchisq:::pmaxchisq_exact(b = b, y = y, x = x)
  expect_true(all(p <= 1 & p >= 0))
})
