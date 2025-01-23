### Does community synchrony change the 

sigma_i <- seq(0,10,length = 10)

RE_synchrony_list<- NULL

for (i in seq_along(sigma_i)){
  
  rmatrix_interest <- lapply(1:100, function(x) simulate_r_matrix(r0 = 0.2, E_0 = 0, 
                                                                  timestep = 365 * 1,
                                                                  sd_envir = 0, seasonal = 30,sigma_i = sigma_i[i]))
  beta_matrix <- simulate_betas(n_species = 10, mean_beta = 9e-4, gamma_shape = 0.8)
  
  RE_synchrony_rep <- NULL
  
  for (r in seq(1,100)){
    
  model_sim <- simulate_full_model(n_species = 10,
                                   param = "standard",
                                   bmatrix = beta_matrix,
                                   rmatrix =  rmatrix_interest[[r]]) 

  model_sim_df <- wrangle_model_output(model_sim, "Yes")
  model_s_df <- wrangle_sus_inf_output(model_sim)[[1]]
  
  
  
  RE_time <- calculate_R_effective(model_sim, time = 365 * 1, b_matrix = beta_matrix, 
                        params = create_parameters(365 *1,'standard'))
  
  sub_RE<- RE_time[RE_time$time > 99,]
  
  average_RE <- mean(sub_RE$RE)
  CV_RE <- sd(sub_RE$RE)/average_RE
  
  synchrony_value <- community.sync(subset(model_sim_df, select = -time))$obs
  
  
  RE_synchrony_rep[[r]] <- cbind.data.frame(synch_value = synchrony_value,
                   avg_RE = average_RE,
                   CV_RE = CV_RE)
   
  }
  RE_synchrony_list[[i]] <- do.call(rbind, RE_synchrony_rep)
}
RE_synchrony_DF <- do.call(rbind, RE_synchrony_list)
RE_synchrony_DF$sigma_i  <- rep(sigma_i, each = 100)


plot_RE_synch_value(RE_synchrony_DF) 

ggplot(RE_synchrony_DF, aes(x = as.factor(sigma_i), y = synch_value )) + geom_boxplot()
