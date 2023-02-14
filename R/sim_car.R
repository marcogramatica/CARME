#' Simulation of proper CAR random effects
#'
#' @export
#' @description
#' `sim_car` returns a vector of CAR distributed random effects
#'
#' @param W Symmetric adjacency matrix of size `n`
#' @param alpha properness parameter between 0 and 1. Defaults to 0.5
#' @param tau marginal precision. Defaults to 5
#'
#' @import MASS
#' @return a vector of length `n`
#' @examples
#'data(W_sel)
#'sim_car(W = W_sel, alpha = .9, tau = 5)
sim_car <- function(W, alpha = .5, tau = 5){
  # Number of Neighbours per area
  D <- diag(rowSums(W))

  # Precision matrix for CAR
  Q <- solve(tau*(D - alpha*W))

  # Generate phi's
  phi <- MASS::mvrnorm(n = 1, mu = rep(0, nrow(W)), Sigma = Q)

  return(phi)
}
