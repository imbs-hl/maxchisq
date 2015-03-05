
pmaxchisq_betensky <- function(b, k, minprop = 0.1, maxprop = 1-minprop) {
  
  ## Compute p value
  r <- maxprop * (1 - minprop) / (minprop * (1 - maxprop))
  p <- dgamma(b, shape = (k - 1)/2, scale = 2) * ((b - k + 1)*log(r) + 4)
  
  ## Return 1 for b==0
  p[b == 0] <- 1
  
  ## Negative pvalues occur for small b, approximation is wrong, pvalue should be 1
  p[p < 0] <- 1
  
  ## If 1-F is rising, approximation is wrong, pvalue should be 1
  next_b <- b + 1e-5 ##.Machine$double.eps
  next_p <- dgamma(next_b, shape = (k - 1)/2, scale = 2) * ((next_b - k + 1)*log(r) + 4)
  p[p < next_p] <- 1
  
  ## Restrict to [0,1]
  pmin(p, 1)
}

pmaxchisq_miller <- function(b, minprop = 0.1, maxprop = 1-minprop) {
  
  ## Compute p value
  r <- maxprop * (1 - minprop) / (minprop * (1 - maxprop))
  db <- dnorm(b)
  p <- 4 * db/b + db * (b - 1/b) * log(r)
  
  ## Return 1 for b==0
  p[b == 0] <- 1
  
  ## Negative pvalues occur for small b, approximation is wrong, pvalue should be 1
  p[p < 0] <- 1
  
  ## If 1-F is rising, approximation is wrong, pvalue should be 1
  next_b <- b + 1e-5 ##.Machine$double.eps
  next_db <- dnorm(next_b)
  next_p <- 4 * next_db/next_b + next_db * (next_b - 1/next_b) * log(r)
  p[p < next_p] <- 1
  
  ## Restrict to [0,1]
  pmin(p, 1)
}

## TODO: Compare speed if only C++ function for test statistic used
## TODO: Possible splits and num_left are computed here and in main function.. But should this work without main function, too?
##' @import Rcpp
pmaxchisq_permutation <- function(b, y, x,
                                  minprop = 0.1, maxprop = 1-minprop, 
                                  num_permutations = 10000) {
  n <- length(y)
  k <- nlevels(y)
  class_counts <- tabulate(y, nbins = k)
  
  x_sorted <- sort(x)
  all_values <- unique(x_sorted)
  quantiles <- quantile(x, c(minprop, maxprop))
  possible_splits <- all_values[all_values > quantiles[1] & all_values < quantiles[2]]
  
  ## For each cutpoint determine last position in x for which x <= cutpoint, 
  ## i.e. the number of observations left of the cutpoint
  num_left <- sapply(possible_splits, function(value) {
    max(which(x_sorted <= value))
  })
  
  ## Run Rcpp permutation function
  pmaxchisq_permutation_internal(b, k, n, num_permutations, as.numeric(y), num_left, class_counts)
}