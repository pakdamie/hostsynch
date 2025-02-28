#' Retrieve parameter to be used
#'
#' Standard sets of parameter values to be used for modeling;
#' Certain parameter set may be only used once, for supplements.
#'
#' @param type (character) "standard" is the usual, "no inf" is if
#' you want to see the dynamics without any infection
#'
#' @returns
#' @export
#'
#' @examples
create_parameters <- function(timesteps, type = "standard") {
  
  ### Standard parameter 
  param_df <-
    c(mu = 2e-2, # Mortality rate
      K = 100, # Carrying capacity
      gamma = 1 / 7, # Recovery rate
      initial_values = 100, #I nitial number.
      delta_T = 1, #Time step 
      timesteps = timesteps,
      infection_time = 1e30)
  

  ### Infection time starts at 501
  param_inf_df <- param_df
  param_inf_df["infection_time"] <- 501
  
    

  # Returns a parameter set depend on
  switch_param <- switch(type,
    "standard" = param_df,
    "inf" = param_inf_df
  )

  return(switch_param)
}

#' Simulate the different parameters for each species' performance curve
#'
#' To create different performance curves (growth rate versus the environment)
#' for each species, we must change the parameters associated with the gaussian
#' function.
#'
#'
#' @param (numeric) sigma What is the mean breadth 
#' @param (numeric) n_species How many species are there?
#'
#' @returns a matrix object
#' @export
#'
#' @examples
#' 
simulate_variablity_species <- function(breadth_var, n_species) {
  
  #Draw the optimum environment for each species using a uniform distribution
  rnorm_Eopt <- runif(n_species, min = -10, max = 10) 
  
  #Draw the niche breadth for each species using a lognormal distribution
  #with mean = 2 (first pass; but might need to account for skewness)
  rlnorm_sigma <- rlnorm(n_species, mean = 2, sd = breadth_var) # std.

  #Given the optimum environment as well as the breadth, optimize
  #to find maximum r that will give an area of 1: 
  #This is the function to optimize
  area_under_fitcurve <- function(r_max, Eopt, Sigma) {
    
    E <- seq(-30, 30, length = 500) #Environment variable
    y_output <- r_max * exp(-((E - Eopt)^2) / (2 * Sigma^2)) #Gaussian curve
    
    #Use approxfun to create a function which we can integrate.
    
    integrated_area <- integrate(approxfun(E, y_output), -30, 30, subdivisions = 2000)$value

    #Optimize for 1
    return(abs(integrated_area) - 5)
  }
  
  #Gives the r-max for species.
  result <- lapply(1:n_species, function(i) {
    uniroot(
      f = area_under_fitcurve,
      interval = c(0, 5),
      Eopt = rnorm_Eopt[i],
      Sigma = rlnorm_sigma[i]
    )$root
  })

  #For each species 
  return(cbind.data.frame(
    rmax = do.call(rbind, result),
    rnorm_Eopt,
    rlnorm_sigma
  ))
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
  environ <- 10 * sin((2 * pi * seq(1, timestep)) / seasonal) +
    rnorm(timestep, mean = 0, sd = sd_envir)
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
    n_species = 10, breadth_var = 0.25,
    sd_envir = 0.15, 
    timestep = 365 * 10, seasonal = 10) {
  
  environmental_factor <- simulate_env_flucs(sd_envir, timestep, seasonal)
  species_trait <- simulate_variablity_species(breadth_var, n_species)


  # For each time-step (col)
  r_mat <- matrix(0, nrow = timestep, ncol = n_species)

  for (e in seq(1, nrow(environmental_factor))) {
    environ_at_t <- environmental_factor[e, "environ"]

    species_r <- apply(species_trait, 1, function(x) simulate_gaussian_curve(environ_at_t, x))

    r_mat[e, ] <- species_r
  }
  return(r_mat)
}

#' A miscellaneous model just for checking the gamma distribution.
#'
#' @param n #How many betas you want to see
#' @param mean_beta #the baseline beta
#' @param CV_desired #the desired coefficient of variation for the gamma
#'
#' @returns
#' @export
#'
#' @examples
simulate_beta_gamma <- function(n, mean_beta, CV_desired) {
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
simulate_betas <- function(n_species = 10, mean_beta, inter_mult, CV_desired) {
  std_wanted <- CV_desired * mean_beta
  varianced_wanted <- std_wanted^2

  # Calculate scale to achieve the target mean
  shape <- (mean_beta^2) / (varianced_wanted)
  rate <- (mean_beta) / (varianced_wanted)

  beta_intra <- rgamma(n_species, shape = shape, rate = rate)

  Beta_Mat <- matrix(0,
    nrow = n_species,
    ncol = n_species
  )

  for (m in seq(1, n_species)) {
    Beta_Mat[, m] <- beta_intra[m] * inter_mult
    Beta_Mat[m, ] <- beta_intra[m] * inter_mult
  }


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
                             infection_time, rmatrix, 
                             n_species = 100, times,
                             initial_values, delta_T) {
  
  # These are the matrix where are tracking the abundance.
  compartment_label <- c(
    "S_mat", "I_mat", "R_mat"
  )
  # Create the matrix where columns are the number of species
  # and the rows are the number of times.
  for (i in 1:length(compartment_label)) {
    assign(
      compartment_label[i],
      matrix(0, nrow = times, ncol = n_species)
    )
  }

  #We initialize all species to have the same initial values
  S_mat[1, ] <- initial_values 

  for (j in seq(1, times - 1)) {
    
    # Total individuals in each species at time j
    N <- S_mat[j, ] + I_mat[j, ] + R_mat[j, ]

    # Ricker births
    new_births <- (rmatrix[j, ] * N * exp(- 1/K * N))

    # How many individuals get infected by other individuals both within and
    # between species

    if ((j == infection_time) == TRUE){
      I_mat[j, ] <- 1 # Initial seed of infection
    }
    else{
      I_mat[j,] <- I_mat[j,]
    }

    new_infections_sp <- S_mat[j, ] * (bmatrix %*% I_mat[j, ])
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
simulate_full_model <- function(n_species, time, params, bmatrix, rmatrix) {
  
  infection_time <- params["infection_time"]
  mu <- params["mu"]
  K <- params["K"]
  gamma <- params["gamma"]
  times <- params["timesteps"]
  initial_values <- params["initial_values"]
  delta_T <- params["delta_T"]

  model <- ricker_SIR_model(
    bmatrix = bmatrix, 
    mu = mu, K = K,
    gamma = gamma,
    infection_time = infection_time,
    rmatrix = rmatrix,
    n_species = n_species, 
    times = times,
    initial_values = initial_values,
    delta_T = delta_T
  )

  return(model)
}
