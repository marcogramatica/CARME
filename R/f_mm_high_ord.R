#' Helper function for MM matrix simulation
#'
#' @keywords internal
#' @import expm
#' @import stats
f_mm_high_ord <- function(
    W, ord, mat_ord_weight, ind_range
){
  n <- nrow(W)
  # Create base matrix of weights
  weight_ext <- matrix(0, ncol = n, nrow = length(ind_range))

  # Compute necessary multiple order of neighbours matrices
  W_cml <- W
  for(i in 3:(ord+1)){ # Starts from 3 bc: 1=identical order ; 2=first order neigh
    # Name for i-th order neighbour matrix
    # name_ind <- paste("W_", i, sep = "")

    # Compute all areas reachable in EXACTLY i-1 steps
    W_cml_current <- W %^% (i-1)
    ## Transform in adj matrix
    W_cml_current[W_cml_current != 0] <- 1
    ## Remove unwanted diagonal
    diag(W_cml_current) <- rep(0, nrow(W_cml_current))

    # Find exactly i-th order neighbour
    W_now <- W_cml_current - W_cml
    ## Adjust matrix
    W_now[W_now < 0] <- 0
    W_now[W_now != 0] <- 1

    # Simulate weights for current neighbour order
    ## Simulate random weights
    W_now[which(W_now != 0)] <- runif(length(which(W_now != 0)))
    for(j in ind_range){
      W_now[j, ] <- (W_now[j,]/sum(W_now[j,]))*mat_ord_weight[j,i]
    }

    # Save global matrix
    weight_ext <- W_now[ind_range,] + weight_ext

    # Move one order forward
    W_cml <- W_cml_current + W_cml
    W_cml[W_cml != 0] <- 1

    # Control when graph is fully connected
    if(length(which(W %^% (i-1) > 1)) == n^2) break
  }

  return(weight_ext)
}
