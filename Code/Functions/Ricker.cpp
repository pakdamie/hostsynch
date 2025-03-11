#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List Ricker_Model(int n_species,
                 int times,
                 double initial_values,
                 arma::mat rmatrix,
                 arma::mat bmatrix,
                 double K,
                 double mu,
                 double gamma,
                 double delta_T) {
  
  
  //SIR matrices
  arma::mat S_mat(times, n_species, arma::fill::zeros); 
  arma::mat I_mat(times, n_species, arma::fill::zeros);
  arma::mat R_mat(times, n_species, arma::fill::zeros);
  
  // Initialize S_mat with initial values
  S_mat.row(0).fill(initial_values);
  
  arma::rowvec N_vec; 
  arma::rowvec newbirths; 
  arma::rowvec newinfection; 
  
  arma::rowvec new_deaths_S;
  arma::rowvec new_deaths_I;
  arma::rowvec new_deaths_R;
  
  arma::rowvec new_recoveries;
  // Total individuals in each species at time j
  for (int j = 0; j < times - 1; j++) {
    
    N_vec = S_mat.row(j) + I_mat.row(j) + R_mat.row(j);
    newbirths = rmatrix.row(j) % N_vec  % exp(-1.0 / K * N_vec);
    newinfection = S_mat.row(j) % (bmatrix * I_mat.row(j).t()).t();
    
    
    arma::rowvec new_deaths_S = mu * S_mat.row(j);
    arma::rowvec new_deaths_I = mu * I_mat.row(j);
    arma::rowvec new_deaths_R = mu * R_mat.row(j);
    arma::rowvec new_recoveries = gamma * I_mat.row(j);  
    
    
    arma::rowvec S_change= newbirths - new_deaths_S - newinfection;
    arma::rowvec I_change = newinfection - new_recoveries - new_deaths_I;
    arma::rowvec R_change= new_recoveries - new_deaths_R;
    
    arma::rowvec total_change_S = S_mat.row(j) + (S_change * delta_T);
    arma::rowvec total_change_I = I_mat.row(j) + (I_change* delta_T);
    arma::rowvec total_change_R = R_mat.row(j) + (R_change* delta_T);
    
    total_change_S = arma::clamp(total_change_S, 0, 1e30);
    total_change_I = arma::clamp(total_change_I, 0, 1e30);
    total_change_R = arma::clamp(total_change_R, 0, 1e30);
    
    
    S_mat.row(j + 1) = total_change_S;
    I_mat.row(j + 1) = total_change_I;
    R_mat.row(j + 1) = total_change_R;
  }
  
  
  return Rcpp::List::create(
    Rcpp::Named("HS") = S_mat,
    Rcpp::Named("HI") = I_mat,
    Rcpp::Named("HR") = R_mat
  );
}

