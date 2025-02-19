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
  theme(legend.position = 'none') ; rmeet

