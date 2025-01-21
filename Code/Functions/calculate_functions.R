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
calculate_R_effective <- function(List, time, param, bmat, mu, gamma) {
  
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
      beta_matrix

    V_mat <- diag(1 / (gamma + mu), ncol = length(S[j, ]), nrow = length(S[j, ]))

    RE <- max(eigen(F_mat %*% V_mat)$values)
    
    contribution_species <- (colSums(F_mat %*% V_mat))

    RE_matrix[[j]] <- cbind.data.frame(time = j, RE = RE)
  }

  return(do.call(rbind, RE_matrix))
}
