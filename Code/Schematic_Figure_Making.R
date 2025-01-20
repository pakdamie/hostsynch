


E <- seq(-10, 10, length = 100)

rmatrix <- generate_variablity_species(0.5, 2.5, 0, 3, 100)

gaussian_curves <- apply(rmatrix,1, function(x) generate_gaussian_curve(E, x), simplify = FALSE) 

for (i in seq(1, 100)){
  gaussian_curves[[i]] <- data.frame(gaussian_curves[[i]])
  gaussian_curves[[i]]$id <- i
}  

gaussian_curves_df <- do.call(rbind, gaussian_curves)

ggplot(gaussian_curves_df, aes( x= E, y = gaus, color = id, group = id))  + 
  geom_line(size = 0.8) + 
  scale_color_viridis() + 
  xlab("Environment") + 
  ylab("Growth rate (r)") + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 14, color = 'black'))

ggsave(here("Figures", "Schematic","gaussian_curve.pdf"), units = 'in', width = 5, height =5)
