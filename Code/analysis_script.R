###How does change in the sd_temp 

damie_lol_whut <- function(){

beta_matrix <- beta_creator(n_species = 20, 
                            mean_beta = 0.1, 
                            sd_norm = 1e-3)

sd_shift <- c(0, 0.01, 0.1, 10,100)

Rmat_test <- r_matrix(n_species = 20, 
                      times = 100, 
                      base_r = 1.2, 
                      sd_shift = 10, 
                      sd_env = 0.5)

Rmat_test_list = NULL
for (i in seq(1,length(sd_shift))){

  Rmat_test_list[[i]]=  r_matrix(n_species =20, 
                                 times = 100, 
                                 base_r = 1.2, 
                                 sd_shift = sd_shift[i], 
                                 sd_env = 0.5)

}


Model_output <- NULL

for (i in seq(1,length(sd_shift))){
  
Model_output[[i]] = ricker_SIR_model(
                   n_species = 20, 
                   times = 100, 
                   rmatrix = Rmat_test_list[[i]], 
                   beta = beta_matrix , 
                   mu = 1.1, 
                   K = 100,
                   gamma = 1e-2, 
                   initial_values = 10, 
                   delta_T = 1)
}

cleaner_df <- function(x){
  full_SIR_DF <- data.frame(x[[1]] +  x[[2]] + x[[3]]   )
return(  full_SIR_DF )
  }
  
  
community_synch <- data.frame(do.call(rbind,
                           lapply(Model_output, 
                                  function(x) community.sync(cleaner_df(x)))))


community_RE <- lapply(Model_output, function(x) 
  calculate_R_effective(x,time = 100, gamma =1e-2, mu = 1.1))

community_synch $CV_RE <- do.call(rbind,lapply(community_RE, function(x) sd(x$RE)/mean(x$RE)))

return(community_synch)
}

rep <- NULL

for (i in seq(1,40)){
  rep[[i]] <- damie_lol_whut()
}

rep_all <- do.call(rbind,rep)
rep_all$sd_phen <- sd_shift 

ggplot(rep_all , 
       aes(x = as.numeric(obs), 
           y = as.numeric(CV_RE),
          color = sd_phen )) +
  scale_color_viridis()+
geom_point(size = 2) + theme_classic() + 
  xlab("Community synchrony") + 
  ylab("Coeffiicent of variation of RE")


