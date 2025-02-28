### What is the community synchrony across niche breadth as well as envrionmental noise
###Vary the environment variability


n_species = 100 
timestep = 365 * 5
seasonal = 30 


sd_envir <-  seq(0,50,10)
#Vary the niche breadth variability
niche_breadth <-  seq(0,0.50, .1)
envir_niche <- expand.grid(sd_envir = sd_envir, niche_breadth = niche_breadth)


r_matrix_list <- NULL

for (i in seq(1, nrow(envir_niche))){

  sd_envir_i <- envir_niche[i,1]
  niche_breadth_i <- envir_niche[i,2]
  
  
  r_matrix_list [[i]] <- lapply(1:100, function(i){
  
  simulate_r_matrix(
    n_species = n_species,
    breadth_var =  niche_breadth_i,
    timestep = timestep,
    sd_envir =  sd_envir_i,
    seasonal = seasonal)
})
}

###After generating.
# Generate the beta matrix
beta_matrix <- simulate_betas(
  n_species = n_species, 
  mean_beta = mean_beta, 
  inter_mult = 0.2, 
  CV_desired =  0.10)



average_list <- NULL
for (i in seq(1, nrow(envir_niche))){
  
  rmat_list <- r_matrix_list[[i]] 
  
  
  reps_list <- NULL
  for (r in seq(1,100)){
  # Simulate the full model
  model_sim <- simulate_full_model(
    n_species = n_species, 
    param = param, 
    time = timestep, 
    bmatrix = beta_matrix, 
    rmatrix = rmat_list[[r]])
  
  model_sim_df <- wrangle_model_output(model_sim, "Yes")
  model_sim_df <- model_sim_df[500:(365*5),]
  
  synchrony_value <- community.sync(subset(model_sim_df, select = -time))$obs
  
  reps_list[[r]] <-  synchrony_value
  
}
 mean_synch<- do.call(rbind, reps_list) |>
   mean()
 
 
 average_list[[i]] <- cbind(mean_synch = mean_synch,envir_niche[i,])
}


average_list_df <- do.call(rbind, average_list )


ggplot(average_list_df, aes(x = niche_breadth, y = sd_envir, fill = mean_synch)) + 
  geom_tile() + 
  scale_fill_viridis(option = 'mako',name = "Synchrony value") +
  xlab("Variability in niche breadth ") + 
  ylab("Variablity in environment") + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme(axis.text = element_text(size = 14, color = 'black'),
        axis.title = element_text(size = 15))

xggsave(here("Figures", "Supp", "Supp_Niche_envir.pdf"), units = "in", width = 5, height =4)
