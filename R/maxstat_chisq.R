
## TODO: Not needed to compute all test statistics?
##' @importFrom Rcpp evalCpp
##' @useDynLib maxchisq
##' @export
maxstat_chisq <- function(y, x, minprop = 0.1, maxprop = 1-minprop, pval_method = "betensky", ...) {
  n <- length(y)
  k <- nlevels(y)
  
  ## TODO: Check parameters
  ## If y is a factor, etc ..
  if (pval_method == "miller" & k != 2) {
    stop("Error: Miller & Siegmund approximation only applicable to 2-class problems.")
  }
  
  ## Possible splits
  x_sorted <- sort(x)
  y_sorted <- y[order(x)]
  all_values <- unique(x_sorted)
  quantiles <- quantile(x, c(minprop, maxprop))
  possible_splits <- all_values[all_values > quantiles[1] & all_values < quantiles[2]]
  
  ## Observed observations left and right of split per class 
  num_left <- sapply(possible_splits, function(cutpoint) {
    tabulate(y_sorted[x_sorted <= cutpoint], nbins = k)
  })
  num_right <- tabulate(y_sorted) - num_left
  
  ##browser()
  ## Test statistics for all cutpoints
  teststats <- sapply(1:length(possible_splits), function(i) {
    chisq_statistic(num_left[, i], num_right[, i], k = k, n = n)
  })
  
  ## P-values for all cutpoints
  if (pval_method == "betensky") {
    pvalues <- pmaxchisq_betensky(b = possible_splits, k = k, minprop = minprop, maxprop = maxprop)
  } else if (pval_method == "miller") {
    pvalues <- pmaxchisq_miller(b = possible_splits, minprop = minprop, maxprop = maxprop)
  } else if (pval_method == "permutation") {
    pvalues <- pmaxchisq_permutation(b = possible_splits, y = y, x = x, minprop = minprop, maxprop = maxprop, ...)
  } else if (pval_method == "none") {
    pvalues <- NA
  } else {
    stop("Error: Unknown pval_method.")
  }
  
  ## Return best cutpoint, test statistic and p-value
  if (pval_method == "none") {
    best_idx <- which.max(teststats)  
  } else {
    best_idx <- which.min(pvalues)  
  }
  list(cutpoint = possible_splits[best_idx], teststat = teststats[best_idx], pvalue = pvalues[best_idx])
}


