
###Example code 

beta_matrix <- beta_creator(n_species = 3, mean_beta = 0.1, sd_norm = 0.5)

Rmat_test <- r_matrix(n_species = 3, 
                      times = 365, 
                      base_r = 1.2, 
                      sd_shift = 10, 
                      sd_env = 0.5)



###Rewrite in rcpp 
model_sim <- ricker_SIR_model(n_species = 3, times = 100, rmatrix = Rmat_test, 
                              beta = beta_matrix , 
                              mu = 1.1, 
                              K = 100,
                              gamma = 1e-2, 
                              initial_values = 10, 
                              delta_T = 1)


full_SIR_DF <- data.frame(model_sim[[1]] +  model_sim[[2]] + model_sim[[3]]   )
full_SIR_DF$time <- seq(1, nrow(full_SIR_DF))
full_SIR_DF_melt <- melt(full_SIR_DF, id.vars = 'time')


ggplot(full_SIR_DF_melt, aes(x = time, y = (value ), color= variable)) + 
  geom_line(size = 0.5, alpha = 0.5) + 
  scale_color_viridis(discrete= TRUE, option = 'turbo') +
  xlab("Time") + 
  ylab("Total Abundance") + 
  theme_classic() + 
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.position = 'none') 


RE_matrix_1 <- calculate_R_effective(model_sim,time = 100, 
                                     bmat = beta_matrix , 
                                     mu = 1.2,
                                     gamma = 1e-2)


ggplot(full_SIR_DF_melt, 
       aes(x = time, 
           y = (value),
           color= variable)) + 
  geom_line(size = 0.5, alpha = 0.5) + 
  scale_color_viridis(discrete= TRUE, option = 'turbo') +
  xlab("Time") + 
  ylab("Total Abundance") + 
  theme_classic() + 
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.position = 'none') +
  ggplot(RE_matrix_1, aes(x = time, y= RE)) + 
  theme_classic() + 
  geom_line() + xlab("Time") + ylab("RE")   + 
  
  ggplot(RE_matrix_1, aes(x = time, y= `1`)) + 
  geom_line() + 
  geom_line(data =RE_matrix_1, aes(x = time, y= `2`), col= 'red') + 
  geom_line(data =RE_matrix_1, aes(x = time, y= `3`), col= 'green') + 
  
  theme_classic() + 
  geom_line() + xlab("Time") + ylab("RE") 



#Thoughts about using the community synchrony metric


time_series <- subset(full_SIR_DF, select = -c(time))
community.sync(time_series)
