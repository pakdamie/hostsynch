### Calculating establishment


# Function to evaluate how community synchronous changes with sigma_i
evaluate_community_synchrony <- function(n_species = 10, 
                                         niche_breadth,
                                         timestep = 365 * 5, 
                                         sd_envir = 10, 
                                         seasonal = 30,
                                         mean_beta = 2e-4,
                                         CV = 0.25, 
                                         param = "standard", 
                                         reps = 100) {
  
  # Initialize the result list
  RE_synchrony_list <- NULL
  
  # Iterate over each sigma_i value
  for (i in seq_along(niche_breadth)) {
    
    # Generate Rmatrices
    rmatrix_list <- lapply(1:reps, function(x) 
      simulate_r_matrix(
        n_species = n_species,
        breadth_var = niche_breadth,
        timestep = timestep,
        sd_envir = sd_envir, 
        seasonal = seasonal
      )
    )
    
    # Generate Bmatrices
    bmatrix_list <- lapply(1:reps, function(x) 
      simulate_betas(
        n_species = n_species, 
        mean_beta = mean_beta, 
        inter_mult = 0.25, 
        CV_desired =  CV)
      )
    
    
    # Initialize the list for replicates
    RE_synchrony_rep <- NULL
    
    # Run simulations for each replicate
    for (r in seq_len(reps)) {
      
      # Simulate the full model
      model_sim <- simulate_full_model(n_species = n_species, 
                                       time = timestep,
                                       params  = create_parameters(timestep, "standard"),
                                       bmatrix =  bmatrix_list[[r]], 
                                       rmatrix = rmatrix_list[[r]])
      
      # Wrangle model output
      model_sim_df <- wrangle_model_output(model_sim, "Yes")
      model_sim_df<- subset(model_sim_df,model_sim_df$time > 499)

      
      # Calculate the effective reproductive number
      RE_time <- calculate_R_effective(model_sim, time = timestep, 
                                       b_matrix = bmatrix_list[[r]], 
                                       params = create_parameters(timestep, "standard"));


      # Extract and process RE values
      sub_RE <- RE_time[RE_time$time > 500, ]
      proportion_invadable <- nrow(sub_RE[sub_RE$RE >=1,])/((365 * 5) - 500)
      average_RE <- max(sub_RE$RE)
      CV_RE <- sd(sub_RE$RE) / average_RE
      
    
      # Calculate community synchrony
      synchrony_value <- community.sync(subset(model_sim_df, select = -time))$obs
      
      time_points_invadable <-  sub_RE[sub_RE$RE >=1,]
      
      persistence_time_LIST <- NULL
      
      if(nrow(time_points_invadable) != 0){
      
      for (m in 1:nrow(time_points_invadable)){
        
        param_RE= create_parameters(timesteps =  timestep, type = "standard")
        param_RE["infection_time"] <- time_points_invadable$time[m]
        
        model_sim_RE <- simulate_full_model(n_species = n_species, 
                                            time = timestep,
                                            params = param_RE,
                                            bmatrix =  bmatrix_list[[r]], 
                                            rmatrix = rmatrix_list[[r]])
        
        RE_time <- calculate_R_effective(model_sim_RE, time = timestep, 
                                         b_matrix = bmatrix_list[[r]], 
                                         params =  param_RE);
        
        sub_RE <- RE_time[RE_time$time >= time_points_invadable$time[m], ]
        
        time_points<- sub_RE[sub_RE$RE <=1,]
        
        persistence_time <- ifelse(nrow(time_points)==0,  timestep,  time_points[1]$time)
        
        persistence_time_LIST[[m]] <- data.frame(start=time_points_invadable$time[m],
                                                 end = persistence_time)
      }
        
        persistence_time <-   do.call(rbind,  persistence_time_LIST)
        persistence_time_mean<- mean(persistence_time$end - persistence_time$start,
                                     na.rm = TRUE)
        persistence_time_sd<- sd(persistence_time$end - persistence_time$start,
                                 na.rm = TRUE)
        
        }
      else{
        
      }
      persistence_time_mean<- NA
      persistence_time_sd<- NA
    
      # Store replicate results
      RE_synchrony_rep[[r]] <- data.frame(
        proportion_invadable =proportion_invadable,
        synch_value = synchrony_value,
        avg_RE = average_RE,
        CV_RE = CV_RE,
        persistence_time_mean = persistence_time_mean,
        persistence_time_sd = persistence_time_sd
        
      )
    }
  
    # Combine replicates for the current sigma_i
    RE_synchrony_list[[i]] <- do.call(rbind, RE_synchrony_rep)
  }
  
  # Combine all results into a single data frame
  RE_synchrony_DF <- do.call(rbind, RE_synchrony_list)
  RE_synchrony_DF$sigma_i <- rep(niche_breadth, each = reps)
  
  return(RE_synchrony_DF)
}


 
RE_synchrony_DF <- evaluate_community_synchrony(niche_breadth = seq(0,0.5, length = 5))


ggplot(RE_synchrony_DF, aes(x =synch_value, 
                            y = persistence_time_mean)) + 
  geom_line()


ggplot(RE_synchrony_DF, aes( x = synch_value, y= avg_RE)) + 
  geom_point() + 
  xlim(0,1) + 
  xlab("Community synchrony") + 
  ylab("Mean R0") + 
  theme_classic() + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15)) + 
 
ggplot(RE_synchrony_DF, aes( x = synch_value, y= proportion_invadable)) + 
  geom_point()  +
  xlim(0,1) + 
  geom_point() + 
  xlab("Community synchrony") + 
  ylab("Proportion of time invadable") + 
  theme_classic() + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15)) +


ggplot(RE_synchrony_DF, aes( x = synch_value, y= CV_RE)) + 
  geom_point()  + geom_point() + 
  xlab("Community synchrony") + 
  ylab("CV of R0") + 
  theme_classic() + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15)) 

