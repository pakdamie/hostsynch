n_species = 10
#Debugging simulate funciton
infection_time = 1e30
mu = 2e-2 # Mortality rate
K = 200 # Carrying capacity
gamma = 1 / 7 # Recovery rate
initial_values = 100
delta_T = 1
timestep = 1000
seasonal = 180
sd_envir = 30
breadth_var = 0.50


sim_r_matrix<- simulate_r_matrix(
  n_species = n_species,
  breadth_var = breadth_var,
  timestep = timestep,
  sd_envir = sd_envir, 
  seasonal = seasonal
)


beta_matrix <- simulate_betas(n_species = n_species, 
                              mean_beta = 3e-4 , 
                              inter_mult = 0.2, 
                              CV_desired = 0.25)


tmp <- ricker_SIR_model(bmatrix = beta_matrix,mu =  mu, K = K,gamma = gamma,
                       infection_time = infection_time, rmatrix = sim_r_matrix, 
                       n_species = n_species,times = timestep, initial_values, delta_T)
 
ricker_rcpp <- Ricker_Model(bmatrix = beta_matrix,mu =  mu, K = K,gamma = gamma,
                        rmatrix = sim_r_matrix, 
                        n_species = n_species,times = timestep, initial_values, delta_T)


wrangle_modeled_df <- wrangle_model_output(tmp, "No") 
wrangle_modeled_tmp2_df <- wrangle_model_output(ricker_rcpp, "No") 


ggplot(wrangle_modeled_df, aes(x = time, y= value, color = variable, group = variable )) +
  geom_line() + 
  xlab("Time") + 
  ylab("TotaAbundance") + 
  scale_color_viridis(discrete = TRUE, option = 'viridis') + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16))

ggplot(wrangle_modeled_tmp2_df , aes(x = time, y= value, color = variable, group = variable )) +
  geom_line() + 
  xlab("Time") + 
  ylab("TotaAbundance") + 
  scale_color_viridis(discrete = TRUE, option = 'mako') + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16))





ggsave(here("Figures","Timeseries.pdf"), units = "in", height = 3, width = 5)




tmp_inf <- tmp[[2]]

ggplot(data = melt(tmp_inf), 
      aes(x = Var1, y= value, color = as.factor(Var2), group = as.factor(Var2) )) + 
  geom_line() + 
  xlab("Time") + 
  ylab("Infected individuals") + 
  scale_color_viridis(discrete = TRUE, option = 'viridis') + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16))
ggsave(here("Figures","Inf_Timeseries.pdf"), units = "in", height = 3, width = 5)




#Experimental debugging script
#Environmental conditions can range from -30 to 30
E <- seq(-50,50,1)

#Does everyone tend to have similar niche breadth?
r_matrix_param <- simulate_variablity_species(breadth_var = breadth_var, n_species = 10)
gaus_curve <- data.frame(apply(r_matrix_param , 1, function(x) simulate_gaussian_curve(E,x)))
gaus_curve$E <- E
gaus_curve_melt <- melt(gaus_curve, id.vars = "E")

ggplot(gaus_curve_melt, aes(x = E, y= value, color = variable)) + 
  geom_line(size = 1.2) + 
  xlab("Time") + 
  ylab("Intrinsic growth rate") + 
  scale_color_viridis(discrete = "TRUE") + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16))

ggsave(here("Figures","niche.pdf"), units = "in", height = 4.5, width = 5)

