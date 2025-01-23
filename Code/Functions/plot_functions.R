#' Plot total abundance of each species over time
#'
#' @param df 
#'
#' @returns
#' @export
#'
#' @examples
plot_total_abundance <- function(df) {
  
  ab_GG <- ggplot(df, aes(x = time, y = (value), color = variable)) +
    geom_line(size = 0.5, alpha = 0.5) +
    scale_color_viridis(discrete = TRUE, option = "turbo") +
    xlab("Time") +
    ylab("Total Abundance") +
    theme_classic() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 15),
      legend.position = "none"
    )

  return(ab_GG)
}


#' Plot the RE against the community synchrony value 
#'
#' @param df 
#'
#' @returns
#' @export
#'
#' @examples
plot_RE_synch_value <- function(df){
  
  avg_RE_GG <- ggplot(df, aes(x = synch_value, y = avg_RE)) + 
    geom_point(aes(color=sigma_i), size = 2) + 
    
    xlab("Community synchrony") + 
    ylab("Average RE") + 
    scale_color_viridis(name = expression(sigma[i])) + 
    theme_classic() + 
    theme(axis.text = element_text(color = 'black', size = 14),
          axis.title = element_text(color = 'black', size = 15))
  
  CV_RE_GG <- ggplot(df, aes(x = synch_value, y = CV_RE)) + 
    geom_point(aes(color=sigma_i)) + 
    xlab("Community synchrony") + 
    ylab("CV RE") + 
    scale_color_viridis(name = expression(sigma[i])) + 
    theme_classic() + 
    theme(axis.text = element_text(color = 'black', size = 14),
          axis.title = element_text(color = 'black', size = 15))
  
  return(avg_RE_GG + CV_RE_GG + plot_layout(guides = 'collect'))
  
}