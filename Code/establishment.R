
sigma_i <- seq(0,30,length = 10)

RE_invadable_list<- NULL

for (i in seq_along(sigma_i)){
  
  rmatrix_interest <- lapply(1:5, function(x) simulate_r_matrix(r0 = 0.02, E_0 = 0, 
                                                                  timestep = 365 * 1,
                                                                  sd_envir = 2, seasonal = 30,
                                                                  sigma_i = sigma_i[i]))
  
  beta_matrix <- simulate_betas(n_species = 10, mean_beta = 9e-4, gamma_shape = 0.8)
  
  
  RE_synchrony_rep <- NULL
  
  for (r in seq(1,5)){
    
    model_sim <- simulate_full_model(n_species = 10,
                                     param = "no_inf",
                                     bmatrix = beta_matrix,
                                     rmatrix =  rmatrix_interest[[r]]) 
    
    model_sim_df <- wrangle_model_output(model_sim, "Yes")

    RE_time <- calculate_R_effective(model_sim, time = 365 * 1, b_matrix = beta_matrix, 
                                     params = create_parameters(365 * 1,"no_inf"))
    
    RE_greater_1 <- nrow(RE_time[RE_time$RE>=1,])/365
    
    synchrony_value <- community.sync(subset(model_sim_df, select = -time))$obs
    
    
    RE_synchrony_rep[[r]] <- cbind.data.frame(synch_value = synchrony_value,
                                              prop_invadable = RE_greater_1 )
    
  }
  RE_invadable_list[[i]] <- do.call(rbind, RE_synchrony_rep)
}
RE_synchrony_DF <- do.call(rbind,  RE_invadable_list)
RE_synchrony_DF$sigma_i  <- rep(sigma_i, each = 100)



ggplot(RE_synchrony_DF, aes( x= synch_value, y= prop_invadable)) + geom_point()




