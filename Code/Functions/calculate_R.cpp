#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat Calculate_R(arma::mat S_mat,
                      int time,
                      arma::mat b_matrix,
                      Rcpp::NumericVector param) {

  
  //Parameters
  double mu = param["mu"];
  double K = param["K"];
  double gamma = param["gamma"];
  int species_n = S_mat.n_cols;
  
  arma::mat RE_mat(time, 1);
  
  for (int j = 0; j < time - 1; j++) {
    
    arma::mat F_mat = diagmat(S_mat.row(j)) * b_matrix;
    
    Rcpp::NumericVector V_vec= Rcpp::rep(1/(gamma + mu), species_n );
    arma::mat V_mat = arma::diagmat(Rcpp::as<arma::vec>(V_vec));  // Convert NumericVector to arma::vec and create diagonal matrix
    
    arma::mat FV = F_mat * V_mat;    
    
    
    arma::cx_vec eigval = eig_gen(FV);
    arma::vec real_eigval = arma::real(eigval);  // Extract real parts
    double eigenvalue = real_eigval.max();  // Get the largest real eigenvalue
    
    
    RE_mat(j,0) = eigenvalue;
  }

  
  return(RE_mat);

}



