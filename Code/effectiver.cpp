#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame calculate_R_effective(List List, int time, NumericMatrix b_matrix, NumericVector params) {
  
  double mu = params["mu"];
  double K = params["K"];
  double gamma = params["gamma"];
  
  List S = List[0];
  std::vector<int> times;
  std::vector<double> RE_values;
  std::vector<double> total_sus_values;
  
  for (int j = 0; j < time; j++) {
    NumericVector S_j = S[j];
    int n = S_j.size();
    NumericMatrix F_mat(n, n);
    NumericMatrix V_mat(n, n);
    
    for (int i = 0; i < n; i++) {
      for (int k = 0; k < n; k++) {
        F_mat(i, k) = (i == k) ? S_j[i] * b_matrix(i, k) : S_j[i] * b_matrix(i, k);
      }
      V_mat(i, i) = 1.0 / (gamma + mu);
    }
    
    NumericMatrix FV = F_mat * V_mat;
    NumericVector eigenvalues = eigen(FV);
    
    double max_RE = 0.0;
    for (double val : eigenvalues) {
      if (std::abs(val) > max_RE) {
        max_RE = std::abs(val);
      }
    }
    
    double total_sus = sum(S_j);
    
    times.push_back(j + 1);
    RE_values.push_back(max_RE);
    total_sus_values.push_back(total_sus);
  }
  
  return DataFrame::create(_["time"] = times,
                           _["RE"] = RE_values,
                           _["total_sus"] = total_sus_values);
}