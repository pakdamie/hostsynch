###Experiment


#Create the beta-matrix of interest
n_species = 10

###Beta matrix - transmission model
beta_matrix <- simulate_betas(n_species, mean_beta = 0.01, sd_beta = 1.5)

###Growth rate
r_matrix <- simulate_r_matrix(n_species, sigma_i = 10, 
                              r0 = 2, E_0 = 0, r_sigma_0 = 5,
                              sd_envir = 1.5, timestep = 365, 
                              seasonal = 1)
###Rewrite in rcpp 
model_sim <- simulate_full_model(n_species,
                                 param = "standard",
                                 bmatrix = beta_matrix,
                                 rmatrix = r_matrix) 

full_SIR_DF<- wrangle_model_output(model_sim, "Yes")



RE_matrix_1 <- calculate_R_effective(model_sim,time = 365, b_matrix = beta_matrix, 
                                     params = create_parameters('standard'))



#Thoughts about using the community synchrony metric


