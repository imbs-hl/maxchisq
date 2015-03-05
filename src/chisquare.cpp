#include <Rcpp.h>

// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// TODO: Use matrix instead of 2 vectors?
// [[Rcpp::export]]
double chisq_statistic(std::vector<double> observed_left, std::vector<double> observed_right, 
                       size_t k, size_t n) {

  // Expected counts
  double colsum_left = std::accumulate(observed_left.begin(), observed_left.end(), 0);
  double colsum_right = std::accumulate(observed_right.begin(), observed_right.end(), 0);
  std::vector<double> expected_left;
  std::vector<double> expected_right;
  expected_left.resize(k, 0);
  expected_right.resize(k, 0);
  for (size_t m = 0; m < k; ++m) {
    size_t rowsum = observed_left[m] + observed_right[m];        
    expected_left[m] = rowsum  * colsum_left / n;
    expected_right[m] = rowsum  * colsum_right / n;
  }
  
  // Compute chi2 test statistic
  double teststat = 0;
  for (size_t m = 0; m < k; ++m) {        
    teststat += (observed_left[m] - expected_left[m]) * (observed_left[m] - expected_left[m]) / expected_left[m];
    teststat += (observed_right[m] - expected_right[m]) * (observed_right[m] - expected_right[m]) / expected_right[m];
  }
  
  return teststat;
}

// [[Rcpp::export]]
std::vector<double> pmaxchisq_permutation_internal(std::vector<double> b, 
                                size_t k, size_t n, size_t num_permutations, 
                                std::vector<size_t> y, 
                                std::vector<size_t> num_left,
                                std::vector<size_t> class_counts) {
  // Use R seed
  Rcpp::RNGScope scope;    
  
  // Initialize pvalue vector
  std::vector<double> pvalues;
  pvalues.resize(b.size(), 0);
    
  // Do permutations
  for (size_t i = 0; i < num_permutations; ++i) {
            
    // Permute response
    std::random_shuffle(y.begin(), y.end(), randWrapper);
    
    // Compute chi2 statistic for each cutpoint and return maximum
    double max_teststat = 0;
    for (size_t j = 0; j < num_left.size(); ++j) {
            
      // Observed counts
      std::vector<double> class_counts_left;
      class_counts_left.resize(k, 0);
      for (size_t l = 0; l < num_left[j]; ++l) {
        ++class_counts_left[y[l]-1];
      }            
      std::vector<double> class_counts_right;
      class_counts_right.resize(k, 0);
      for (size_t m = 0; m < k; ++m) {
        class_counts_right[m] = class_counts[m] - class_counts_left[m];
      }
      
      // Test statistic
      double teststat = chisq_statistic(class_counts_left, class_counts_right, k, n);
      
      // Save if new maximum
      if (teststat > max_teststat) {
        max_teststat = teststat;
      }
    }
    
    // If teststat is at least as extreme as b, count for pvalue
    for (size_t j = 0; j < b.size(); ++j) {
      if (max_teststat >= b[j]) {
        ++pvalues[j];
      } 
    }
  }
  
  // P value is proportion of tests at least as extreme as observed
  for (size_t j = 0; j < pvalues.size(); ++j) {
    pvalues[j] /= (double) num_permutations;
  }
  return pvalues;
}


