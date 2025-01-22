### Does community synchrony change the 

sigma_i <- seq(0,30,length = 10)

RE_synchrony_list<- NULL

for (i in seq_along(sigma_i)){
  
  rmatrix_interest <- lapply(1:100, function(x) simulate_r_matrix(sigma_i = sigma_i[i]))
  beta_matrix <- simulate_betas(n_species, mean_beta = 0.01, sd_beta = 1.5)
  
  RE_synchrony_rep <- NULL
  
  for (r in seq(1,100)){
    
  model_sim <- simulate_full_model(n_species = 10,
                                   param = "standard",
                                   bmatrix = beta_matrix,
                                   rmatrix =  rmatrix_interest[[r]]) 

  model_sim_df <- wrangle_model_output(model_sim, "Yes")

  RE_time <- calculate_R_effective(model_sim, time = 365 * 5, b_matrix = beta_matrix, 
                        params = create_parameters('standard'))
  
  average_RE <- max(RE_time$RE)
  CV_RE <- sd(RE_time$RE)/average_RE
  
  synchrony_value <- community.sync(subset(model_sim_df, select = -time))$obs
  
  
  RE_synchrony_rep[[r]] <- cbind.data.frame(synch_value = synchrony_value,
                   avg_RE = average_RE,
                   CV_RE = CV_RE)
   
  }
  RE_synchrony_list[[i]] <- do.call(rbind, RE_synchrony_rep)
}
RE_synchrony_DF <- do.call(rbind, RE_synchrony_list)
RE_synchrony_DF$sigma_i  <- rep(sigma_i, each = 5)


plot_RE_synch_value(RE_synchrony_DF) 
