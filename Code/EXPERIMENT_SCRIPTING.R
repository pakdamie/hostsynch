rmatrix <- simulate_variablity_species(0.50, 10)

E <- seq(-50, 50, length = 1500)

rmeet <- data.frame(apply(rmatrix, 1, function(rmat_row) {
  simulate_gaussian_curve(E, rmat_row)
}), E = E) |>
  melt(id.vars = "E") |>
  ggplot(aes(x = E, y = value, color = variable)) +
  geom_line() +
  xlab("Environment") +
  ylab("Growth rate") +
  scale_color_viridis(discrete = TRUE) +
  theme_classic() +
  theme(legend.position = "none")
rmeet




calculate_establishment_proportion <- function(
    n_species, timestep, bmatrix, rmatrix_list) {
  summary_list <- NULL


  for (i in seq(1, length(breadth_rep_list))) {
    model_sim <- lapply(1:reps, function(x) {
      simulate_full_model(
        n_species = n_species,
        params = create_parameters(timestep, "standard"),
        bmatrix = bmatrix,
        rmatrix = rmatrix_list[[i]][[x]]
      )
    })


    RE_time <- lapply(model_sim, function(x) {
      data.frame(
        time = 500:(timestep),
        RE = Calculate_R(x,
          time = timestep,
          b_matrix = bmatrix,
          param = create_parameters(timestep, "standard")
        )[500:(timestep), ]
      )
    })


    RE_values <- do.call(rbind, lapply(RE_time, function(x) {
      cbind(
        mean_RE = mean(x$RE),
        CV_RE = sd(x$RE) / mean(x$RE)
      )
    }))

    RE_greater_1 <- do.call(
      rbind,
      lapply(RE_time, function(x) {
        nrow(x[x$RE > 1, drop = FALSE, ]) / (timestep - 500) # Fixed syntax
      })
    )


    model_sim_df <- lapply(model_sim, function(x) {
      wrangle_model_output(x, "Yes")
    })

    synchrony_value <- do.call(
      rbind,
      lapply(model_sim_df, function(x) {
        community.sync(subset(x, select = -time))$obs
      })
    )





    summary_list[[i]] <- cbind.data.frame(
      synch_value = synchrony_value,
      niche = full_expand[i, 1],
      envir = full_expand[i, 2],
      RE_values,
      RE_greater_1
    )
  }

  return(do.call(rbind, summary_list))
}
