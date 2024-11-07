library(reshape2)
library(ggplot2)
library(viridis)
library (synchrony)

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
r_matrix <- function(n_species = 100, times, base_r, sd_shift = 1, sd_env = 1){
  
  #Initialize an empty matrix to be filled with growth rates
  rmat <- matrix(0, nrow = times, ncol = n_species )
  
  # For each species, we give them an intrinsic 'shift' 
  
  for (k in seq(1, n_species)){
    shift = rnorm(1, mean = 0, sd = sd_shift )
    
    rmat[ , k]<- base_r + sin(seq(1, times)/365 + shift) + 
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
  new_births <- N * exp(rmatrix[j, ] * (1-(N/K)))
  

  #How many individuals get infected by other individuals both within and
  #between species
  new_infections_sp <- matrix(0, ncol = 1, nrow =   n_species)
  
  for (sp in 1:(n_species)){
    
    new_infections_sp[sp] <- sum(beta * S_mat[j, sp] * I_mat[j, ])
  }
  
  new_deaths_S <- mu * S_mat[j, ]
  new_deaths_I <- mu * I_mat[j, ]
  new_deaths_R <- mu * R_mat[j, ]
  new_recoveries <- gamma * I_mat[j, ]
  
  S_change <- new_births + new_deaths_S - new_infections_sp
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

Rmat_test <- r_matrix(n_species = 10, 
                      times = 100, 
                      base_r = 3.5, 
                      sd_shift = 10, 
                      sd_env = 0.5)

###Rewrite in rcpp 
model_sim <- ricker_SIR_model(n_species = 10, times = 100, rmatrix = Rmat_test, 
                             beta = 0.1, mu = 0.5,K = 50,
                             gamma = 0.05, 
                             initial_values = 100, 
                             delta_T = 1)


full_SIR_DF <- data.frame(model_sim[[1]]+ model_sim[[2]] + model_sim[[3]])
full_SIR_DF$time <- seq(1, nrow(full_SIR_DF))
full_SIR_DF_melt <- melt(full_SIR_DF, id.vars = 'time')


ggplot(full_SIR_DF_melt, aes(x = time, y = log10(value), color= variable)) + 
  geom_line(size = 0.9) + 
  scale_color_viridis(discrete= TRUE, option = 'turbo') +
  xlab("Time") + 
  ylab("Total Abundance") + 
  theme_classic() + 
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.position = 'none') 

#Thoughts about using the community synchrony metric


time_series <- subset(full_SIR_DF, select = -c(time))
community.sync(time_series)


