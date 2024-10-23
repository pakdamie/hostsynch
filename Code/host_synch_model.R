#' Calculate the r for all species over time
#'
#' @param n_species The number of species 
#' @param times How long the model will run for 
#' @param base_r What is the baseline growth rate for all species 
#' @param sd_shift The SD of the normal distribution (how far does the
#' species deviate)
#' @param sd_env The SD of the environmental noise
#' 
#' @return A matrix where the columns are the species and the rows are time
#' @export
#'
#' @examples r_matrix(100, 100, 2.5, 1,1)
r_matrix <- function(n_species, times, base_r, sd_shift = 1, sd_env = 1){
  
  #Initialize an empty matrix to be filled with growth rates
  rmat <- matrix(0,nrow = times, ncol = n_species )
  
  for (k in seq(1, nspecies)){
    shift = rnorm(1, mean = 0, sd = sd_shift )
    
    rmat[,k]<- (base_r + sin(seq(1,times)/365 + shift)) + 
               rnorm(times, mean = 0, sd =sd_env )
  }
  return(rmat)
}


#' Simulate the discrete SIR/Ricker model
#'
#' @param n_species The number of species 
#' @param times  How long the model will run for 
#' @param rmatrix The matrix from `r_matrix `
#' @param beta The transmission coefficient 
#' @param mu Mortality coefficient
#' @param K The carrying capacity
#' @param gamma Recovery rate
#' @param initial_values The number of susceptible individuals to start out with
#' @param delta_T The time step
#'
#' @return A list of three matrixes (S, I and R)
#' @export
#'
#' @examples
#' 
ricker_SIR_model <- function(
  n_species, times, rmatrix, 
  beta, mu, d ,K ,
  gamma, initial_values, delta_T){

  #These are the matrix where are tracking the abundance
  compartment_label <- c(
    "S_mat", "I_mat", "R_mat"
  )
  
 #Create the matrix where columns are the number of species
 # and the rows are the number of times
  for (i in 1:length(compartment_label)) {
    assign(
      compartment_label[i],
      matrix(0, nrow = times, ncol = n_sp))
  }

  S_mat[1,] <- initial_values #Initialize everyone the same
  I_mat[1,] <- 1 #initial seed of infection
  
  for (j in seq(1, times - 1)){

  #Total individuals in each species at time j
  N <- S_mat[j,] + I_mat[j,] + R_mat[j,]
  
  new_births <- sapply(N * exp(rmatrix[j,] * (1-N)/K), function(x) rpois(lambda = x, n= 1))

  #WHow many individuals get infected by other individuals both within and
  #between species
  new_infection_sp <- matrix(0, ncol = 1, nrow = n_sp)
  
  for (sp in 1:(n_sp)){
    new_infection_sp[sp] <- sum(beta * S_mat[j,sp] * I_mat[j,])
  }
  
  new_deaths_S <- mu * S_mat[j,]
  new_deaths_I <- mu * I_mat[j,]
  new_deaths_R <- mu * R_mat[j,]
  new_recoveries <- gamma* I_mat[j,]
  
  S_change <- new_births + new_deaths_S - new_infections_sp
  I_change <- new_infections_sp - new_recoveries - new_deaths_I
  R_change <- new_recoveries - new_deaths_R
  

  S_mat[j +1, ] <- S_mat[j, ] +   (S_CHANGE * delta_T)
  I_mat[j +1, ] <- I_mat[j, ] +   (I_CHANGE * delta_T)
  R_mat[j +1, ] <- R_mat[j, ] +   (R_CHANGE * delta_T)
  
  S_mat [S_mat < 0] <- 0
  I_mat [I_mat < 0] <- 0
  R_mat[ R_mat < 0] <- 0

  
  
  }
  return(list(S_mat, I_mat, R_mat))  
}
