
library(maxchisq)

x <- iris$Sepal.Length
y <- iris$Species

## Possible splits
minprop <- 0.1
maxprop <- 0.9
n <- length(y)
k <- nlevels(y)
x_sorted <- sort(x)
y_sorted <- y[order(x)]
all_values <- unique(x_sorted)
quantiles <- quantile(x, c(minprop, maxprop))
possible_splits <- all_values[all_values >= quantiles[1] & all_values < quantiles[2]]

## Compute test statistics
teststats <- sapply(possible_splits, function(split) {
  idx <- x <= split
  observed <- cbind(table(y[idx]), table(y[!idx]))
  chisq.test(observed)$statistic
})

## Get maximum
possible_splits[which.max(teststats)]
max(teststats)

## Use new function
maxstat_chisq(y = y, x = x, pval_method = "betensky")

## Apply to data frame
maxstat_chisq_test(Species ~ ., iris, pval_method = "betensky")

## Only one covar
maxstat_chisq_test(Species ~ Petal.Length, iris, pval_method = "betensky")
