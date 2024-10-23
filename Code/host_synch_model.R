
ricker_SIR_model <- function(n_sp = 100, times = 100,
   rmat2 =rmat,   beta = 0.01, mu = 2, d = 0.1, K = 5000 ,gamma = 0.005,initial_values = 10,
    delta_T = 0.5){

  compartment_label <- c(
    "S_mat", "I_mat", "R_mat"
  )
  
  for (i in 1:length(compartment_label)) {
    assign(
      compartment_label[i],
      matrix(0, nrow = times, ncol = n_sp))
  }

  S_mat[1,] <- initial_values
  I_mat[1,] <- 1
  
  for (j in seq(1, times - 1)){

  N <- S_mat[j,] + I_mat[j,] + R_mat[j,]
  
  new_births <- sapply(N * exp(rmat2[j,] * (1-d) * (1-N)/K), function(x) rpois(lambda = x, n= 1))

  new_infection_sp <- matrix(0, ncol = 1,nrow = n_sp)
  
  for (sp in 1:(n_sp)){
    
    new_infection_sp[sp] <- sum(sapply(beta * S_mat[j,sp] * I_mat[j,], function(x) rpois(lambda = x, n= 1)))
  }
  
  new_deaths_S <- sapply(mu * S_mat[j,], function(x) rpois(n =1, lambda = x))
  new_deaths_I <-  sapply(mu * I_mat[j,], function(x) rpois(n =1, lambda = x))
  new_deaths_R <- sapply(mu * R_mat[j,], function(x) rpois(n =1, lambda = x))
  
  new_recoveries <- sapply(gamma* I_mat[j,], function(x) rpois(n =1, lambda = x))
  

  
  S_CHANGE <- new_births + new_deaths_S - new_infections_sp
  I_CHANGE <- new_infections_sp -   new_recoveries -   new_deaths_I
  R_CHANGE <- new_recoveries -   new_deaths_R
  
  S_CHANGE[is.na(S_CHANGE) == TRUE] <- 0
  
  S_mat[j +1,] <- S_mat[j,] +   (S_CHANGE * delta_T)
  I_mat[j +1,] <- I_mat[j,] +   (I_CHANGE * delta_T)
  R_mat[j +1,] <- R_mat[j,] +   (R_CHANGE * delta_T)
  
  S_mat [S_mat  < 0] <- 0
  I_mat [I_mat < 0] <- 0
  R_mat[ R_mat < 0] <- 0

  
  
  }
  return(list(S_mat, I_mat, R_mat))  
    
    
}

test<- ricker_SIR_model()
full <- data.frame(test[[1]]+ test[[2]] + test[[3]])
full$time <- seq(1,100)
full_df <- melt(full, id.vars = 'time')

ggplot(full_df, aes(x = time, y =(value), color = variable, group = variable)) + 
  geom_line() + 
  scale_color_viridis(option = 'viridis', discrete = TRUE) + 
  theme_classic() + 
  theme(legend.position= 'none')


r_matrix <- function(nspecies, times, base_r, differences){

rmat <- matrix(0,nrow = times, ncol = nspecies )
for (k in seq(1,nspecies)){
shift = rnorm(1, mean = 0, sd =2)
  
rmat[,k]<-  2.5 + sin(seq(1,1000)/365 + shift)+ rnorm(1000,mean = 0,sd =0.1)
}
}
