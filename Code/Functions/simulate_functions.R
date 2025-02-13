#' Retrieve parameter data.frame
#'
#' Sets of parameter values for the modeling
#'
#'
#' @param type (character) If you want the standard type
#'
#' @returns
#' @export
#'
#' @examples
create_parameters <- function(timesteps, type = "standard") {
  
  param_df <-
    c(infection_time = 50,
      mu = 0.01, # Mortality rate
      K = 100,  # Carrying capacity
      gamma = 1/7, #Recovery rate 
      initial_values = 100,
      delta_T = 1,
      time = timesteps
    )
  
  
  #there is no infection
  param_no_inf <-param_df 
  param_no_inf["infection_time"] <- 1e10
  
  switch_param<- switch(type,
         "standard" = param_df,
         "no_inf" = param_no_inf)
  
  
  return(switch_param)
}

#' Simulate the different parameters for each species' performance curve
#'
#' To create different performance curves (growth rate versus the environment)
#' for each species, we must change the parameters associated with the gaussian
#' function.
#'
#'
#' @param (numeric) sigma What is the standard deviation for all of the normal distribution?
#' @param (numeric) r0 What is the baseline growth rate for the community?
#' @param (numeric) E_0 What is the baseline optimum environment for the community?
#' @param (numeric) r_sigma_0 What is the standard deviation for the response in the community?
#' @param (numeric) n_species How many species are there?
#'
#' @returns a matrix object
#' @export
#'
#' @examples
simulate_variablity_species <- function(sigma_i, r0, E_0, r_sigma_0, n_species) {
  rnorm_rmax <- rnorm(n_species, mean = r0, sd = 2) # growth rate
  rnorm_Eopt <- rnorm(n_species, mean = E_0, sd = sigma_i) # optimum Environment
  rnorm_r_sigma <- rnorm(n_species, mean = r_sigma_0, sd = 1) # std.

  rmatrix <- cbind(r0, rnorm_Eopt,r_sigma_0)

  return(rmatrix)
}

#' Generate the growth rate at the given environment value
#'
#' Assuming a Gaussian function, given the parameters of the species,
#' what is the growth rate of the species at Environment condition (time t)
#'
#' @param E (numeric) environment
#' @param rmat_row
#'
#' @returns
#' @export
#'
#' @examples
simulate_gaussian_curve <- function(E, rmat_row) {
  r0 <- rmat_row[1]
  E_opt <- rmat_row[2]
  r_sigma <- rmat_row[3]

  gaus <- r0 * exp(-((E - E_opt)^2) / (2 * r_sigma^2))
  return(gaus)
}

#' Simulate the environmental fluctuations independently
#'
#' @param sd_envir (numeric)
#' @param timestep (numeric)
#' @param seasonal (numeric)
#'
#' @returns
#' @export
#'
#' @examples
simulate_env_flucs <- function(sd_envir, timestep = 365 * 5, seasonal = 10) {
  
  environ <- 10 * sin((2*pi*seq(1, timestep)) / seasonal) + rnorm(timestep, mean = 0, sd = sd_envir)
  return(cbind(seq(1:timestep), environ))
}


#' Calculate the growth rate (r) for all species over time
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
#' @examples r_matrix(100, 100, 2.5, 1, 1)
simulate_r_matrix <- function(
    n_species = 10, sigma_i = 1, r0 = 2, E_0 = 2, r_sigma_0 = 5,
    sd_envir = 0.15, timestep = 365 * 10, seasonal = 10) {
  
  environmental_factor <- simulate_env_flucs(sd_envir, timestep, seasonal)
  species_trait <- simulate_variablity_species(sigma_i, r0, E_0, r_sigma_0, n_species)

  
  # For each time-step (col)
  r_mat <- matrix(0, nrow = timestep, ncol = n_species)

  for (e in seq(1, nrow(environmental_factor))) {
    environ_at_t <- environmental_factor[e, "environ"]

    species_r <- apply(species_trait, 1, function(x) simulate_gaussian_curve(environ_at_t, x))

    r_mat[e,] <- species_r
  }
  return(r_mat)
}

simulate_beta_gamma <- function(n, mean_beta, CV_desired){
  
  std_wanted <- CV_desired * mean_beta
  
  varianced_wanted <- std_wanted^2
  
  # Calculate scale to achieve the target mean
  shape <- (mean_beta^2) / (varianced_wanted)
  rate <- (mean_beta) / (varianced_wanted)
  
  return(dgamma(seq(0, 9e-4, length.out = n), shape = shape, rate = rate))
  
}



#' Simulate the beta matrix for disease transmission using a gamma
#' distribution
#'
#' @param n_species Number of species needed
#' @param mean_beta The community base-line disease transmission
#' @param sd_beta The standard deviation.
#'
#' @returns
#' @export
#'
#' @examples
simulate_betas <- function(n_species = 10, mean_beta, CV_desired) {
  
  
  std_wanted <- CV_desired * mean_beta
  varianced_wanted <- std_wanted^2
  
  # Calculate scale to achieve the target mean
  shape <- (mean_beta^2) / (varianced_wanted)
  rate <- (mean_beta) / (varianced_wanted)

  beta_intra <- rgamma(n_species, shape = shape, rate = rate)
  
  beta_inter <- beta_intra -
    runif(n_species, min = 0, max = beta_intra)


  Beta_Mat <- matrix(beta_inter,
    nrow = n_species,
    ncol = n_species
  )

  diag(Beta_Mat) <- beta_intra

  
  return(Beta_Mat)
}



#' Simulate the discrete SIR/Ricker model
#'
#' #Note to damie, recode to get n_species and times variable
#' from the rmatrix itself; also put all parameters into a dataframe
#' for less verbose function inputs
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
ricker_SIR_model <- function(bmatrix, mu, K, gamma, 
                             infection_time, rmatrix, n_species, times,
                             initial_values, delta_T) {
  # These are the matrix where are tracking the abundance
  compartment_label <- c(
    "S_mat", "I_mat", "R_mat"
  )

  # Create the matrix where columns are the number of species
  # and the rows are the number of times
  for (i in 1:length(compartment_label)) {
    assign(
      compartment_label[i],
      matrix(0, nrow = times, ncol = n_species)
    )
  }

  S_mat[1, ] <- initial_values # Initialize everyone the same
  

  for (j in seq(1, times - 1)) {
    # Total individuals in each species at time j
    N <- S_mat[j, ] + I_mat[j, ] + R_mat[j, ]

    # Ricker births
    new_births <- (N * exp(rmatrix[j, ] * (1 - (N/K)))) - N

    # How many individuals get infected by other individuals both within and
    # between species
    
    if(j == infection_time){
    I_mat[j, sample(1:n_species,1) ] <- 10 # Initial seed of infection
    }
    
    new_infections_sp <-  S_mat[j, ] * (bmatrix %*% I_mat[j, ])

    new_deaths_S <- mu * S_mat[j, ]
    new_deaths_I <- mu * I_mat[j, ]
    new_deaths_R <- mu * R_mat[j, ]
    new_recoveries <- gamma * I_mat[j, ]

    S_change <- new_births - new_deaths_S - new_infections_sp
    I_change <- new_infections_sp - new_recoveries - new_deaths_I
    R_change <- new_recoveries - new_deaths_R


    S_mat[j + 1, ] <- S_mat[j, ] + (S_change * delta_T)
    I_mat[j + 1, ] <- I_mat[j, ] + (I_change * delta_T)
    R_mat[j + 1, ] <- R_mat[j, ] + (R_change * delta_T)

    S_mat[j + 1, ][S_mat[j + 1, ] < 0] <- 0
    I_mat[j + 1, ][I_mat[j + 1, ] < 0] <- 0
    R_mat[j + 1, ][R_mat[j + 1, ] < 0] <- 0
  }
  return(list(S_mat, I_mat, R_mat))
}


#' Simulate full model 
#' @param n_species
#' @param rmatrix 
#' @param bmatrix - Disease transmission
#' @returns
#' @export
#'
#' @examples
simulate_full_model <- function(n_species,param_type, bmatrix, rmatrix) {
  
  params <- create_parameters(365 *1,param_type)

  infection_time <- params["infection_time"]
  mu <- params["mu"]
  K <- params["K"]
  gamma <- params["gamma"]
  times <- params["time"]
  initial_values <- params["initial_values"]
  delta_T <- params["delta_T"]

  model <- ricker_SIR_model(
    bmatrix, mu, K, gamma, infection_time, rmatrix, n_species, times,
    initial_values, delta_T
  )

  return(model)
}
