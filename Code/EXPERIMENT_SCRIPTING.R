### Play around with the environmental factor
#sd_envir, timestep, and seasonal period 
environmental_factor <- simulate_env_flucs(1, 365 * 1, 180)
plot(environmental_factor,type = 'l')


growth_rate_df<- simulate_variablity_species(sigma_i = 10, r0 = 0.2, 
                                             E_0 = 0, r_sigma_0 = 4, 10)|>
  apply(1, function(x) simulate_gaussian_curve(x,E = seq(-50,50)))|>
  data.frame()

growth_rate_df$E <-  seq(-50,50)
growth_rate_melt <- melt(growth_rate_df, id.vars ='E')

ggplot(growth_rate_melt, aes( x= E, y= value, group = variable, color = variable)) + 
  geom_line() + 
  scale_color_viridis(discrete = TRUE) + 
  theme_classic() + 
  xlab("Environment") + 
  ylab("Growth Rate") + 
  geom_vline(xintercept = 0) +
  theme(legend.position = 'none')


rmatrix_interest <- simulate_r_matrix(r0 = 0.10, E_0 = 0, sigma_i = 10, 
                                      timestep = 365 * 1,
                                      sd_envir = 10, seasonal = 10)

####Play around with the beta matrix
beta_matrix <- simulate_betas(n_species = 10, mean_beta = 9e-4, CV_desired =  0.01)
hist(beta_matrix)


##Play around with the full model
full_model <-  simulate_full_model(n_species = 10,
                                 param = "standard",
                                 bmatrix = beta_matrix,
                                 rmatrix =  rmatrix_interest) 

wrangle_model_output(full_model,synch = "No") |>
  plot_total_abundance() +

(wrangle_sus_inf_output(full_model)[[1]] |>
  plot_total_abundance()/

wrangle_sus_inf_output(full_model)[[2]] |>
  plot_total_abundance())

model<- wrangle_model_output(full_model, "No") 
agg_model <- aggregate(model$value, list(model$time), 'sum')

ggplot(agg_model, aes(x = Group.1, y = x)) + geom_line()

RE <- calculate_R_effective(full_model,
                            time = 365 * 1, beta_matrix,
                            create_parameters(365, "standard"))

