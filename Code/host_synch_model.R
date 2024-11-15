library(reshape2)
library(ggplot2)
library(viridis)
library (synchrony)

beta_creator <- function(n_species = 100){
  beta_interspecific = 0.01
  beta_intraspecific = 0.05

  Beta_Mat <- matrix(beta_interspecific, nrow = n_species, ncol = n_species)  
  diag(Beta_Mat) <- beta_intraspecific
  
  return(Beta_Mat)
}



#' Calculate the growth rate (r) for all species over time
#'
#'
#'
#' @param n_species The number of species (default = 100)
#' @param times How long the model will run for 
#' @param base_r What is the baseline growth rate for all species 
#' @param sd_shift The SD of the normal distribution (how much does the
#' species deviate)
#' @param sd_env The SD of the environmental noise (how much does the species
#' deviate)

#' @return A matrix with columns as species and the rows as time
#' @export
#'
#' @examples r_matrix(100, 100, 2.5, 1,1)
r_matrix <- function(n_species = 100, times, base_r, sd_shift = 1, sd_env = 1,seasonal=1){
  
  #Initialize an empty matrix to be filled with growth rates
  rmat <- matrix(0, nrow = times, ncol = n_species )
  
  # For each species, we give them an intrinsic 'shift' 
  
  for (k in seq(1, n_species)){
    shift = rnorm(1, mean = 0, sd = sd_shift)
    
    rmat[ , k] <- base_r + sin((seq(1, times) -  shift)/seasonal) + 
               rnorm(times, mean = 0, sd = sd_env)
  }
  return(rmat)
}


#' Simulate the discrete SIR/Ricker model
#'
#'#Note to damie, recode to get n_species and times variable
#'from the rmatrix itself; also put all parameters into a dataframe
#'for less verbose function inputs
#'
#' @param n_species The number of species (default = 100)
#' @param times  How long the model will run for 
#' @param rmatrix The matrix from `r_matrix `
#' @param beta The transmission coefficient 
#' @param mu Mortality coefficient
#' @param K The carrying capacity
#' @param gamma Recovery rate
#' @param initial_values The number of susceptible individuals to start out with
#' @param delta_T The time step
#'
#' @return A list of three matrices (S, I and R)
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
      matrix(0, nrow = times, ncol =   n_species))
  }

  S_mat[1, ] <- initial_values # Initialize everyone the same
  I_mat[1, ] <- 1 # Initial seed of infection
  
  for (j in seq(1, times - 1)){

  #Total individuals in each species at time j
  N <- S_mat[j, ] + I_mat[j, ] + R_mat[j, ]
  
  
  #Ricker births 
  
  new_births <- (N * exp(rmatrix[j, ] * (1-(N/K))))
  
  
  #How many individuals get infected by other individuals both within and
  #between species
  
  # TD: I think this is odd, as it's treating cross-species transmission with
  # the same weight as within-species transmission. Perhaps beta should be a 
  # matrix of values, with diag being within-species transmission and off-diag
  # being cross-species transmission (often much much lower values
  #How many individuals get infected by other individuals both within and
  #between species
  

  new_infections_sp <-colSums(diag(I_mat[j,]) %*% (diag(S_mat[j,]) %*% beta))
  
  new_deaths_S <- mu * S_mat[j, ]
  new_deaths_I <- mu * I_mat[j, ]
  new_deaths_R <- mu * R_mat[j, ]
  new_recoveries <- gamma * I_mat[j, ]
  
  S_change <- new_births - new_deaths_S - new_infections_sp
  I_change <- new_infections_sp - new_recoveries - new_deaths_I
  R_change <- new_recoveries - new_deaths_R
  

  S_mat[j +1, ] <- S_mat[j, ] +   (S_change * delta_T)
  I_mat[j +1, ] <- I_mat[j, ] +   (I_change * delta_T)
  R_mat[j +1, ] <- R_mat[j, ] +   (R_change * delta_T)
  
  S_mat [S_mat < 0] <- 0
  I_mat [I_mat < 0] <- 0
  R_mat[ R_mat < 0] <- 0


  }
  return(list(S_mat, I_mat, R_mat))  
}

###Example code 

beta_matrix <- beta_creator (n_species = 10)

Rmat_test <- r_matrix(n_species = 10, 
                      times = 365, 
                      base_r = 1.2, 
                      sd_shift = 0.01, 
                      sd_env = 0.1)
image(Rmat_test)



###Rewrite in rcpp 
model_sim <- ricker_SIR_model(n_species = 
                              10, times = 100, rmatrix = Rmat_test, 
                              beta = beta_matrix , mu =1.2,K =1000,
                              gamma = 1e-2, 
                              initial_values = 10, 
                              delta_T = 1)


full_SIR_DF <- data.frame(model_sim[[1]] )
full_SIR_DF$time <- seq(1, nrow(full_SIR_DF))
full_SIR_DF_melt <- melt(full_SIR_DF, id.vars = 'time')


ggplot(full_SIR_DF_melt, aes(x = time, y = (value ), color= variable)) + 
  geom_line(size = 0.5, alpha = 0.5) + 
  scale_color_viridis(discrete= TRUE, option = 'turbo') +
  xlab("Time") + 
  ylab("Total Abundance") + ylim(600, 1200) + 
  theme_classic() + 
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.position = 'none') 

#Thoughts about using the community synchrony metric


time_series <- subset(full_SIR_DF, select = -c(time))
community.sync(time_series)


