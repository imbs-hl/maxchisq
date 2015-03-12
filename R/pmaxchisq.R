
##' P-value estimation for a maximally selected chi squared statistic using the approximation by Betensky & Rabinowitz (1999). 
##'
##' @title Betensky & Rabinowitz p-value estimation
##' @param b Quantile, i.e., the value of the test statistic
##' @param k Number of rows in contingency table
##' @param minprop Lower quantile for cutpoint 
##' @param maxprop Upper quantile for cutpoint 
##' @return P-value estimation
##' @author Marvin N. Wright
##' @references
##'   Betensky, R. A., & Rabinowitz, D. (1999). Maximally Selected \eqn{\chi2} Statistics for k×2 Tables. Biometrics 55, 317-320.
pmaxchisq_betensky <- function(b, k, minprop = 0.1, maxprop = 1-minprop) {
  ## Internal function, no argument checks!
  
  ## Compute p value
  r <- maxprop * (1 - minprop) / (minprop * (1 - maxprop))
  p <- dgamma(b, shape = (k - 1)/2, scale = 2) * ((b - k + 1)*log(r) + 4)

  ## Return 1 for b==0
  p[b == 0] <- 1

  ## Negative pvalues occur for small b, approximation is wrong, pvalue should be 1
  p[p < 0] <- 1

  ## If 1-F is rising, approximation is wrong, pvalue should be 1
  next_b <- b + 1e-5
  next_p <- dgamma(next_b, shape = (k - 1)/2, scale = 2) * ((next_b - k + 1)*log(r) + 4)
  p[p < next_p] <- 1

  ## Restrict to [0,1]
  pmin(p, 1)
}

##' P-value estimation for a maximally selected chi squared statistic using the approximation by Miller & Siegmund (1982). 
##'
##' @title Miller & Siegmund p-value estimation
##' @param b Quantile, i.e., the value of the test statistic
##' @param minprop Lower quantile for cutpoint 
##' @param maxprop Upper quantile for cutpoint 
##' @return P-value estimation
##' @author Marvin N. Wright
##' @references
##'   Miller, R., & Siegmund, D. (1982). Maximally selected chi square statistics. Biometrics 38, 1011-1016.
pmaxchisq_miller <- function(b, minprop = 0.1, maxprop = 1-minprop) {
  ## Internal function, no argument checks!
  
  ## Compute p value
  r <- maxprop * (1 - minprop) / (minprop * (1 - maxprop))
  db <- dnorm(b)
  p <- 4 * db/b + db * (b - 1/b) * log(r)

  ## Return 1 for b==0
  p[b == 0] <- 1

  ## Negative pvalues occur for small b, approximation is wrong, pvalue should be 1
  p[p < 0] <- 1

  ## If 1-F is rising, approximation is wrong, pvalue should be 1
  next_b <- b + 1e-5
  next_db <- dnorm(next_b)
  next_p <- 4 * next_db/next_b + next_db * (next_b - 1/next_b) * log(r)
  p[p < next_p] <- 1

  ## Restrict to [0,1]
  pmin(p, 1)
}

## TODO: Compare speed if only C++ function for test statistic used
## TODO: Possible splits and num_left are computed here and in main function.. But should this work without main function, too?
##' P-value estimation for a maximally selected chi squared statistic using permutations.
##' 
##' @title Permutation p-value estimation
##' @param b Quantile, i.e., the value of the test statistic
##' @param y Response
##' @param x Covariate
##' @param minprop Lower quantile for cutpoint 
##' @param maxprop Upper quantile for cutpoint 
##' @param num_permutations Number of permutations used
##' @return P-value estimation
##' @author Marvin N. Wright
#' @import Rcpp
pmaxchisq_permutation <- function(b, y, x,
                                  minprop = 0.1, maxprop = 1-minprop,
                                  num_permutations = 10000) {
  ## Internal function, no argument checks!
  
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
