n_species <- 10
timestep <- 365 * 5
reps <- 1000
breadth_var <- seq(0, 0.5, 0.1)
seasonal_var <- 180

full_expand<- expand.grid(breadth_var,seasonal_var)


breadth_rep_list <- NULL
for (i in seq(1, nrow(full_expand))) {
  breadth_rep_list[[i]] <- lapply(1:reps, function(x) {
    simulate_r_matrix(
      n_species = n_species,
      breadth_var = full_expand[i,1],
      timestep = timestep,
      sd_envir = 10,
      seasonal = full_expand[i,2]
    )
  })
}

bmatrix_list <- NULL
bmatrix_list <- lapply(1:reps, function(x) {
  simulate_betas(
    n_species = n_species,
    mean_beta = 2e-4,
    inter_mult = 0.25,
    CV_desired = 0.10
  )
})

bmatrix_interest<-  
  simulate_betas(
    n_species = n_species,
    mean_beta = 2.5e-4,
    inter_mult = 0.25,
    CV_desired = 0.10
  )





calculate_establishment_proportion<- function(
    n_species,timestep, bmatrix, rmatrix_list){
  
  
  summary_list = NULL
  for (i in seq(1,length(breadth_rep_list))){
    
    model_sim <- lapply(1:reps, function(x)
      simulate_full_model(n_species = n_species,
                          params = create_parameters(timestep,"standard"),
                          bmatrix = bmatrix,
                          rmatrix = rmatrix_list[[i]][[x]])) 
    
    
    RE_time <- lapply(model_sim, function(x) calculate_R_effective(x,
               time = timestep, 
               b_matrix = bmatrix, 
               params = 
              create_parameters(timestep,"standard"))[500:(timestep),])
    
    
    RE_greater_1 <- lapply(RE_time, function(x) {
      valid_rows <- x[x$RE > 1, ]  # Subset rows where RE > 1
      
      cbind(
        nrow(valid_rows) / (timestep - 499),  # Compute proportion
        mean(valid_rows[,"total_sus"], na.rm = TRUE)  # Compute mean safely
      )
    })
    
    model_sim_df <- lapply(model_sim, function(x)
      wrangle_model_output(x, "Yes"))
    
    synchrony_value <- lapply(model_sim_df,function(x)
      community.sync(subset(x, select = -time))$obs)
    
    
    summary_list[[i]] <- cbind.data.frame(synch_value = do.call(rbind,synchrony_value),
                                          prop_invadable = do.call(rbind,RE_greater_1)[,1],
                                          total_sus = do.call(rbind,RE_greater_1)[,2],
                                          niche = full_expand[i,1],
                                          envir = full_expand[i,2])
    
  }
  
  return(do.call(rbind,summary_list))
}


calculate_establishment<- function(n_species,timestep, 
                                   bmatrix_list, 
                                   rmatrix_list){
  summary_list = NULL
  for (i in seq(1,nrow(full_expand))){
    
  model_sim <- lapply(1:reps, function(x)
                      simulate_full_model(n_species = n_species,
                                   params = create_parameters(timestep,"standard"),
                                   bmatrix = bmatrix_list[[x]],
                                   rmatrix = rmatrix_list[[i]][[x]])) 
  

  RE_time <- lapply(model_sim, function(x) calculate_R_effective(x,
                                   time = timestep, 
                                   b_matrix = bmatrix_list[[x]], 
                                   params = 
                                   create_parameters(timestep,"standard"))[500:
                                                                          (timestep),])
  
  RE_values <- do.call(rbind,lapply(RE_time, function(x) cbind(mean_RE =mean(x$RE),
                                    CV_RE = sd(x$RE) / mean(x$RE))))
  
  
  RE_greater_1 <- lapply(RE_time,function(x)
                       c(nrow(x[x$RE>1,])/(timestep - 500)))
                       

  model_sim_df <- lapply(model_sim, function(x)
                         wrangle_model_output(x, "Yes"))
  
  synchrony_value <- lapply(model_sim_df,function(x)
                            community.sync(subset(x, select = -time))$obs)
  
  
  summary_list[[i]] <- cbind.data.frame(synch_value = do.call(rbind,synchrony_value),
                                            prop_invadable = do.call(rbind,RE_greater_1),
                                            niche = full_expand[i,1],
                                            envir = full_expand[i,2],
                                            RE_values)
  
  }
  
  return(do.call(rbind,summary_list))
  }


tmp <- calculate_establishment(n_species = 10 ,
                         timestep = 365 * 5, 
                         bmatrix_list = bmatrix_list, 
                         rmatrix_list = breadth_rep_list)

tmp2 <- calculate_establishment_proportion(
  n_species = 10 ,
  timestep = 365 * 5, 
  bmatrix = bmatrix_interest, 
  rmatrix_list = breadth_rep_list)


gam.check(gam1)
gam1<-gam(prop_invadable~s(synch_value) ,data =tmp2, 
         family = betar(link = "logit")) 
summary(gam1)
lot(gam1, se = TRUE, rug = TRUE)



ggplot(tmp2, aes(x = synch_value, y = prop_invadable)) + geom_point() 
  
ggplot(tmp, aes( x= synch_value, y=  prop_invadable)) +
  geom_point() + 
  facet_wrap(~envir) +
  xlab("Community synchrony") + 
  ylab("Proportion of time invasible") + 
  theme_classic()


ggplot(tmp, aes( x= synch_value, y= CV_RE)) +
  geom_point() + 
  facet_wrap(~envir) + 
  xlab("Community synchrony") + 
  ylab("Proportion of time invasible") + 
  theme_classic()

  
  ggplot(tmp, aes( x= synch_value, y= mean_RE)) +
    geom_point() + 
    xlab("Community synchrony") + 
    ylab("Proportion of time invasible") + 
    theme_classic()
  
  
  
  
  


