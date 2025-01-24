# Function to evaluate how community synchrony changes with sigma_i
evaluate_community_synchrony <- function(n_species = 10, 
                                         sigma_i_seq,
                                         r0 = 0.10,
                                         timestep = 365, 
                                         sd_envir = 10, 
                                         seasonal = 1, 
                                         mean_beta = 1e-3, 
                                         CV = 0.10, 
                                         E_0 = 0,
                                         param = "standard", reps = 100) {
  
  # Initialize the result list
  RE_synchrony_list <- NULL
  
  # Iterate over each sigma_i value
  for (i in seq_along(sigma_i_seq)) {
    
    # Generate r matrices for the current sigma_i
    rmatrix_list <- lapply(1:reps, function(x) simulate_r_matrix(
      n_species = n_species,
      r0 = r0, E_0 = E_0, timestep = timestep, sd_envir = sd_envir, 
      seasonal = seasonal, sigma_i = sigma_i_seq[i]
    ))
    
    # Generate the beta matrix
    beta_matrix <- simulate_betas(n_species = n_species, mean_beta = mean_beta, 
                                  inter_mult = 0.2, CV_desired =  CV)
    
    # Initialize the list for replicates
    RE_synchrony_rep <- NULL
    
    # Run simulations for each replicate
    for (r in seq_len(reps)) {
      # Simulate the full model
      model_sim <- simulate_full_model(n_species = n_species, param = param, 
                                       bmatrix = beta_matrix, rmatrix = rmatrix_list[[r]])
      
      # Wrangle model output
      model_sim_df <- wrangle_model_output(model_sim, "Yes")
      
      # Calculate the effective reproductive number
      RE_time <- calculate_R_effective(model_sim, time = timestep, b_matrix = beta_matrix, 
                                       params = create_parameters(timestep, param))
      
      # Extract and process RE values
      sub_RE <- RE_time[RE_time$time > 49, ]
      average_RE <- mean(sub_RE$RE)
      CV_RE <- sd(sub_RE$RE) / average_RE
      
      # Calculate community synchrony
      synchrony_value <- community.sync(subset(model_sim_df, select = -time))$obs
      
      # Store replicate results
      RE_synchrony_rep[[r]] <- data.frame(
        synch_value = synchrony_value,
        avg_RE = average_RE,
        CV_RE = CV_RE
      )
    }
    
    # Combine replicates for the current sigma_i
    RE_synchrony_list[[i]] <- do.call(rbind, RE_synchrony_rep)
  }
  
  # Combine all results into a single data frame
  RE_synchrony_DF <- do.call(rbind, RE_synchrony_list)
  RE_synchrony_DF$sigma_i <- rep(sigma_i_seq, each = reps)
  
  return(RE_synchrony_DF)
}

RE_synchrony_DF <- evaluate_community_synchrony(sigma_i_seq = c(5,30,50))


plot_RE_synch_value(RE_synchrony_DF) 

ggplot(RE_synchrony_DF, aes(x = as.factor(sigma_i), y = synch_value )) + geom_boxplot()
