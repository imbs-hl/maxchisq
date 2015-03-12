library(maxchisq)
context("maxstat_chisq")

n <- 50
x <- rnorm(n)
y <- factor(round(rbinom(n, 1, 0.5)))

test_that("maxstat_chisq object has desired structure", {
  betensky <- maxstat_chisq(y = y, x = x, pval_method = "betensky")
  expect_is(betensky, "numeric")
  expect_equal(length(betensky), 3)
  expect_named(betensky, c("cutpoint", "teststat", "pvalue"))
})

test_that("maxstat_chisq_test object has desired structure", {
  betensky <- maxstat_chisq_test(Species ~ ., iris, pval_method = "betensky")
  expect_is(betensky, "list")
  expect_equal(length(betensky), 3)
  expect_named(betensky, c("covariate", "cutpoint", "pvalue"))
})
