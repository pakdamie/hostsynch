calculate_estab_persist <- function(n_species = 10,
                                    timestep,
                                    full_time,
                                    bmatrix,
                                    rmatrix) {
  
  
  #Simulate the original model (no infection) to calculate R0
  model_sim <- simulate_full_model(
    n_species = n_species,
    params = create_parameters(timestep, "standard"),
    bmatrix = bmatrix,
    rmatrix = rmatrix)
  
  # Calculate the R0 over time
  RE_time <- data.frame(
    time = 500:timestep,
    RE = Calculate_R(model_sim[[1]],
                     time = timestep,
                     b_matrix = bmatrix,
                     param = create_parameters(timestep, "standard")
    )[500:timestep]
  )
  
  # Calculate mean and CV of R0
  R0_values <- data.frame(
    mean_RE = mean(RE_time$RE, na.rm = TRUE),
    CV_RE = sd(RE_time$RE, na.rm = TRUE) / mean(RE_time$RE, na.rm = TRUE)
  )
  
  # Proportion of time steps where the system is invadable (R0 ≥ 1)
  R0_invadable <- nrow(RE_time[RE_time$RE >= 1, , drop = FALSE]) / (timestep - 500)
  
  # Identify time points where invasion is possible
  time_invadable <- RE_time$time[RE_time$RE >= 1]
  
  # If invasion is possible, calculate the average infection duration
  if (length(time_invadable) > 0) {
    infection_lengths <- NULL
    for (i in seq_along(time_invadable)) {
      
      invaded_time <- time_invadable[i]
      mod_params <- create_parameters(  full_time, "standard")
      mod_params["infection_time"] <- invaded_time
      
      model_sim_inf <- simulate_full_model(
        n_species = n_species,
        params = mod_params,
        bmatrix = bmatrix,
        rmatrix = rmatrix
      )
     
      
      R_calc <- Calculate_R(model_sim_inf[[1]],
                            time = full_time,
                            b_matrix = bmatrix,
                            param = mod_params)
      
      R_calc_df <- cbind.data.frame(time = 1:full_time, RE = R_calc)
      
      Post_R <- R_calc_df[invaded_time:full_time,]
      
      index_inf <- Post_R[which(Post_R $RE < 1),]
      infection_lengths[[i]] <- ifelse(length(index_inf) > 0, 
                                     index_inf$time[1], 
                                     full_time)
    }
    tmp <- do.call(rbind,infection_lengths)
    avg_length <- mean(tmp, na.rm = TRUE)
  } else {
    avg_length <- NA
  }
  
  # Process model output
  model_sim_df <- wrangle_model_output(model_sim, "Yes")
  
  # Compute synchrony value
  synchrony_value <- community.sync(subset(model_sim_df, select = -time))$obs
  
  # Create summary dataframe
  summary_df <- data.frame(
    synch_value = synchrony_value,
    R0_values,
    prop_invadable = R0_invadable,
    inf_length = avg_length
  )
  
  return(summary_df)
}
simulate_estab_persist <- function(full_expand,
                                   rmat_list,
                                   reps = reps, 
                                   timestep = 365 + 500,
                                   full_time = 365 * 2,
                                   n_species = 10,
                                   mean_beta = 2.5e-4,
                                   CV_desired = 0.30) {
  
  collected_info_list <- NULL # Preallocate list

  for (i in seq_len(nrow(full_expand))) {
    rmats <- rmat_list[[i]]

    
    collected_info_df <- NULL # Preallocate list for reps

    for (r in 1:reps) {
      
      rmat_interest <- as.matrix(rmats[[r]])

      
      bmat <- simulate_betas(
        n_species = n_species,
        mean_beta = mean_beta,
        inter_mult = 0.30,
        CV_desired = CV_desired
      )


      collected_info_df[[r]] <-
        calculate_estab_persist(n_species,
          timestep,
          full_time,
          bmatrix = bmat,
          rmatrix = rmat_interest
        )
    }

    # Convert the list of data frames to a single data frame
    collect_df <- do.call(rbind, collected_info_df)

    # Add metadata columns
    collect_df$breadth_var <- full_expand[i, 1]
    collect_df$seasonal <- full_expand[i, 2]

    # Store in main list
    collected_info_list[[i]] <- collect_df
  }

  return(do.call(rbind, collected_info_list))
}



reps <- 100
timestep <- 365 + 500
breadth_var <- c(0.25)
seasonal_var <- c(180)
full_expand <- expand.grid(breadth_var, seasonal_var)

breadth_rep_list <- NULL
for (i in seq(1, nrow(full_expand))) {
  breadth_rep_list[[i]] <- lapply(1:reps, function(x) {
    simulate_r_matrix(
      n_species = 10,
      breadth_var = full_expand[i, 1],
      timestep = timestep,
      sd_envir = 10,
      seasonal = full_expand[i, 2]
    )
  })
}


tmpo4 <- simulate_estab_persist(full_expand = full_expand, 
                              reps = reps,
                              rmat_list = breadth_rep_list)





ggplot(tmpo4, aes(x = synch_value, y = mean_RE)) +
  geom_point() + 
  xlab("Community synchrony") + 
  ylab("Mean RE") + 
  theme_classic()

ggplot(tmpo4, aes(x = synch_value, y = prop_invadable,color = seasonal)) +
  geom_point()
ggplot(tmpo4, aes(x = synch_value, y = CV_RE)) +
  geom_point() 
ggplot(tmpo4, aes(x = synch_value,y =inf_length )) + geom_point()
