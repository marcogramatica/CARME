#' Simulation of MM matrix based
#'
#' @export
#' @description
#' `sim_MM_matrix` returns a multiple membership matrix simulated based on an
#' adjacency matrix according to the method described in
#'
#' @param W Symmetric adjacency matrix of size `n`
#' @param m Integer. Number of membership to simulate
#' @param ord Integer. Maximum order of neighbours to be used to simulate the
#' memberships based on the adjacency matrix `W`
#' @param w_ord A vector of length `ord` that specifies the weights of each
#' order of neighbours
#' @param id_vec Vector of zeros and ones of length `n`. Defaults to
#' a vector of ones. It indicates whether an area is included in the
#' simulation of a membership
#' @param excess_areas if different from FALSE it indicates the indices of the
#' areas to reuse in simulating memberships, whenever `m` > `n`. It defaults to
#' FALSE, and if omitted randomly selects without replacement
#' (if `m` - `n` <= `n`, otherwise with replacement) a subset of areas
#' @param red_areas vector of indices of areas to use if `m` < `n`
#' 
#' @references Marco Gramatica. Silvia Liverani. Peter Congdon. 
#' "Structure Induced by a Multiple Membership Transformation on the Conditional
#'  Autoregressive Model." Bayesian Analysis Advance Publication 1 - 25, 2023.
#'  https://doi.org/10.1214/23-BA1370
#'
#' @return an m x n matrix of weights
#' @import expm
#' @import stats
#'
#' @examples
#' set.seed(455)
#'
#' #---- Load data
#' data(W_sel)
#' ## Number of areas
#' n <- nrow(W_sel)
#' ## Number of memberships
#' m <- 153
#'
#' #---- Simulate MM matrix
#' w_ord <- c(.5, .35, .15) # Weight of each neighbours orders
#' ord <- length(w_ord) - 1 # Order of neighbours to include
#' H_sel_sim <- sim_MM_matrix(
#'   W = W_sel, m = m, ord = ord, w_ord = w_ord, id_vec = rep(1, nrow(W_sel))
#' )
sim_MM_matrix <- function(
    W, m, ord = 3, w_ord, id_vec, excess_areas = FALSE, red_areas
){
  # Number of neighbours
  n <- nrow(W) ; if(n < 4) stop("Too few areas")

  # Create matrix of weights
  # weight_ext <- matrix(0, ncol = n, nrow = n)

  ### Create weights for each new membership
  ### If the area identically included in the membership does not renormalise
  ### otherwise it will have to

  # General checks for weights
  ## First order (is area identically included in the membership?)
  if(exists("id_vec") == FALSE) id_vec <- rep(1, n)
  ## Length of neighbour order weights vector
  if(length(w_ord) < (ord + 1)) stop("Not enough weights")
  ## Force 0 on weight for identical area weight if necessary
  if(sum(id_vec) == 0) w_ord[1] <- 0
  ## Check weights sum to 1
  if(sum(w_ord[1:(ord+1)]) != 1){
    warning(
      "Neighbour weights do not sum to one. Function will normalise and proceed"
    )
    w_ord <- abs(w_ord)/sum(w_ord)
  }

  # Identical area inclusion weight
  weights_W1 <- diag(id_vec)*w_ord[1]

  # Normalise weights for membership without identical area
  ## Create register of weights for higher orders
  mat_ord_weight <- cbind(
    id_vec*w_ord[1], # Weights of identical order
    matrix(rep(w_ord[-1], n), nrow = n, ncol = length(w_ord[-1]), byrow = T)
  )
  ## Normalise register
  mat_ord_weight <- mat_ord_weight/rowSums(mat_ord_weight)

  # Compute necessary multiple order of neighbours matrices
  W_cml <- W
  weight_ext <- f_mm_high_ord(
    W, ord = ord, mat_ord_weight = mat_ord_weight, ind_range = 1:n
  )

  # Add Identical neighbours and first order
  ## First order neigh
  weights_W2 <- W
  weights_W2[which(weights_W2 != 0)] <- runif(length(which(weights_W2 != 0)))
  for(j in 1:n){
    weights_W2[j, ] <- (weights_W2[j,]/sum(weights_W2[j,]))*mat_ord_weight[j,2]
  }
  ## Finalise
  weight_ext <- weight_ext + weights_W1 + weights_W2 ; rowSums(weight_ext)
  ## Check
  if(sum(rowSums(weight_ext)) != n){
    stop("Weights of weight_ext do not sum to 1")
  }

  # More MEMBERSHIPS than AREAS (m > n)
  if(m > n){
    # Check
    if(missing(excess_areas)){
      if(m - n > n){
        warning("m - n > n, so excess areas sampled with replacement")
        excess_areas <- sample(1:n, m - n, replace = T)
      } else{
        excess_areas <- sample(1:n, m - n)
      }
    }

    # Create additional matrix
    ## Higher order
    excess_weight_high <- f_mm_high_ord(
      W, ord = ord, mat_ord_weight = mat_ord_weight,
      ind_range = excess_areas)
    ## First order
    excess_weight_2 <- W
    excess_weight_2[which(excess_weight_2 != 0)] <-
      runif(length(which(excess_weight_2 != 0)))
    for(j in excess_areas){
      excess_weight_2[j, ] <-
        (excess_weight_2[j,]/sum(excess_weight_2[j,]))*mat_ord_weight[j,2]
    }
    ## Subset to excess areas in second order
    excess_weight_2 <- excess_weight_2[excess_areas,]
    ## Identical areas
    excess_weight_1 <- matrix(0, nrow = length(excess_areas), ncol = n)
    excess_weight_1[cbind(1:length(excess_areas), excess_areas)] <- w_ord[1]
    ## Finalise
    excess_weight <- excess_weight_1 + excess_weight_2 + excess_weight_high

    # Final save
    weight_ext <- rbind(weight_ext, excess_weight)

    ## Check
    if(sum(rowSums(weight_ext)) != m){
      stop("Weights of weight_ext do not sum to 1")
    }
  }

  # More AREAS than MEMBERSHIPS (m < n)
  if(m < n){
    weight_ext <- weight_ext[sample(1:n, m), ]
  }

  # Return values
  # output <- list(
  #   weight_ext = weight_ext,
  #   w_ord = w_ord,
  #   id_vec = id_vec,
  #   weights_W1 = weights_W1,
  #   weights_W2 = weights_W2,
  #   W_cml = W_cml
  # )
  return(weight_ext)
}
