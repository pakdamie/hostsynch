## How is community synchrony influenced by seasonal fluctuations 
## as well as the variability in niche breadth?

library(parallel)

# Detect number of cores
numCores <- detectCores()
numCoresToUse <- max(1, numCores - 2)  # Ensure there's at least 1 core available

n_species <- 10
timestep <- 365 * 5
reps <- 5
breadth_var <- seq(0, 0.5, 0.1)
seasonal_var <- seq(1, 180, length = 5)

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
