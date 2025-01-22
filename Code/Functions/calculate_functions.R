
#' Wrangle model output
#'
#' @param list_output 
#'
#' @returns
#' @export
#'
#' @examples
wrangle_model_output <- function(list_output, synch){
  full_SIR_DF <- data.frame(list_output[[1]] +  list_output[[2]] + list_output[[3]])
  full_SIR_DF$time
  
  
  full_SIR_DF$time <- seq(1, nrow(full_SIR_DF))
  full_SIR_DF_melt <- melt(full_SIR_DF, id.vars = 'time')
  
  df_to_return <- switch(synch, "Yes" = full_SIR_DF, "No" = full_SIR_DF_melt)
  
  
  return(df_to_return)
  
}

#' Calculate the RE
#'
#' @param List 
#' @param time 
#' @param param 
#' @param bmat 
#' @param mu 
#' @param gamma 
#'
#' @returns
#' @export
#'
#' @examples
calculate_R_effective <- function(List, time, b_matrix, params) {
  
  mu <- params["mu"]
  K <- params["K"]
  gamma <- params["gamma"]

  
  
  S <- List[[1]]

  RE_matrix <- NULL

  for (j in (1:time)) {
    F_mat <- matrix(
      rep(
        S[j, ],
        length(S[j, ])
      ),
      ncol = length(S[j, ]), byrow = TRUE
    ) *
      b_matrix

    V_mat <- diag(1 / (gamma + mu), ncol = length(S[j, ]), nrow = length(S[j, ]))

    RE <- max(abs(Re(eigen(F_mat %*% V_mat)$values)))
    
    contribution_species <- (colSums(F_mat %*% V_mat))

    RE_matrix[[j]] <- cbind.data.frame(time = j, RE = RE)
  }

  return(do.call(rbind, RE_matrix))
}

#' Calculate the community synchrony
#'
#' @param df 
#'
#' @returns
#' @export
#'
#' @examples
calculate_comm_synch <- function(df){
}