
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
  
  ## Test statistics for all cutpoints
  teststats <- sapply(1:length(possible_splits), function(i) {
    chisq_statistic(num_left[, i], num_right[, i], k = k, n = n)
  })
  
  best_idx <- which.max(teststats)  
  best_cutpoint <- possible_splits[best_idx]
  best_teststat <- teststats[best_idx]
  
  ## P-value for best cutpoint
  if (pval_method == "betensky") {
    pvalue <- pmaxchisq_betensky(b = best_cutpoint, k = k, minprop = minprop, maxprop = maxprop)
  } else if (pval_method == "miller") {
    pvalue <- pmaxchisq_miller(b = best_cutpoint, minprop = minprop, maxprop = maxprop)
  } else if (pval_method == "permutation") {
    pvalue <- pmaxchisq_permutation(b = best_cutpoint, y = y, x = x, minprop = minprop, maxprop = maxprop, ...)
  } else if (pval_method == "none") {
    pvalue <- NA
  } else {
    stop("Error: Unknown pval_method.")
  }
  
  ## Return best cutpoint, test statistic and p-value
  c(cutpoint = best_cutpoint, teststat = best_teststat, pvalue = pvalue)
}

##' @export
maxstat_chisq_test <- function(formula, data, na.action, ...) {
  ## TODO: Check parameters
  
  formula <- formula(formula)
  if (class(formula) != "formula") {
    stop("Error: Invalid formula.")
  }
  
  ## Apply formula
  data_model <- model.frame(formula, data, na.action = na.action)
  
  ## Apply maximally selected chi2 statistic to each covariate
  maxstats <- data.frame(t(sapply(data_model[, -1], maxstat_chisq, y = data_model[, 1], ...)))
  
  ## Select covariate with minimal p value
  best_idx <- max.col(t(-maxstats$pvalue))
  
  ## Return best covariate, cutpoint and p value
  list(covariate = rownames(maxstats)[best_idx], 
    cutpoint = maxstats[best_idx, ]$cutpoint, 
    pvalue = maxstats[best_idx, ]$pvalue)
}


