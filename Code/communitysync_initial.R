## How is community synchrony influenced by seasonal fluctuations 
## as well as the variability in niche breadth?


source("Code/Functions/calculate_functions.R")
source("Code/Functions/simulate_functions.R")
sourceCpp("Code/Functions/Ricker.cpp")
library(parallel)

# Detect number of cores
numCores <- detectCores()
numCoresToUse <- max(1, numCores - 2)  # Ensure there's at least 1 core available

n_species <- 10
timestep <- 365 * 5
reps <- 1000
breadth_var <- seq(0, 0.5, length =10)
seasonal_var <- seq(1, 180, length = 10)

# Create a grid for breadth and seasonal variables
full_expand <- expand.grid(breadth_var, seasonal_var)

# Initialize list to store results
breadth_rep_list <- NULL

# Run parallel simulations
for (i in seq(1, nrow(full_expand))) {
  breadth_rep_list[[i]] <- mclapply(1:reps, function(x) {
    simulate_r_matrix(
      n_species = n_species,
      breadth_var = full_expand[i, 1],
      timestep = timestep,
      sd_envir = 10,
      seasonal = full_expand[i, 2]
    )
  }, mc.cores = numCoresToUse)  # Use the correct number of cores
}

beta_matrix <- simulate_betas(n_species = n_species, 
                              mean_beta = 3e-100, 
                              inter_mult = 0.2, 
                              CV_desired = 0.25)


calculate_synchrony  <- function(
    n_species,timestep,bmatrix,rmatrix_list){
  
  
  summary_list = NULL
  for (i in seq(1,length(breadth_rep_list))){
    
    model_sim <- lapply(1:reps, function(x)
      simulate_full_model(n_species = n_species,
                          params = create_parameters(timestep,"standard"),
                          bmatrix = bmatrix,
                          rmatrix = rmatrix_list[[i]][[x]])) 
    
    
    model_sim_df <- lapply(model_sim, function(x)
      wrangle_model_output(x, "Yes"))
    
    synchrony_value <- lapply(model_sim_df,function(x)
      community.sync(subset(x, select = -time))$obs)
    
    
    summary_list[[i]] <- 
      cbind.data.frame(synch_value = do.call(rbind,synchrony_value),
                                          niche = full_expand[i,1],
                                          envir = full_expand[i,2])
    
  }
  
  return(do.call(rbind,summary_list))
}

full_synchrony<- calculate_synchrony(n_species ,timestep , beta_matrix,breadth_rep_list )

aggregated_dat<- aggregate(synch_value~niche + envir, mean, data = full_synchrony)


ggplot(aggregated_dat, aes(x = niche, y = envir, fill = synch_value))+ 
  geom_tile() + 
  scale_fill_viridis(name = "Community \nsynchrony") + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  xlab("Variability of niche breadth") + 
  ylab("Environmental fluctuation") + 
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'))

ggsave(here("Figures","initial_synchrony.pdf"), width = 6, height =4.5, units = "in")
