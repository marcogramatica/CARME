#' CAR-MM prior model
#'
#' @export
#' @param d_list List of data inputs for the stan model.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#' 
#' @importFrom Rdpack reprompt
#' @references
#' {
#' Marco Gramatica. Silvia Liverani. Peter Congdon. 
#' "Structure Induced by a Multiple Membership Transformation on the Conditional
#'  Autoregressive Model." Bayesian Analysis Advance Publication 1 - 25, 2023.
#'  https://doi.org/10.1214/23-BA1370
#'  
#'  Petrof, O, Neyens, T, Nuyts, V, Nackaerts, K, Nemery, B, Faes, C. On the 
#'  impact of residential history in the spatial analysis of diseases with a 
#'  long latency period: A study of mesothelioma in Belgium. 
#'  Statistics in Medicine. 2020; 39: 3840– 3866.
#'  https://doi.org/10.1002/sim.8697
#' }
#'
#' @examples
#'\donttest{
#' set.seed(455)
#'
#' #---- Load data
#' data(W_sel)
#' ## Number of areas
#' n <- nrow(W_sel)
#' ## Number of memberships
#' m <- 153
#'
#' #---- Simulate covariates
#' X <- cbind(rnorm(nrow(W_sel)), rnorm(nrow(W_sel)))
#' ## Min-max normalisation
#' X_cent <- apply(X, 2, function(x) (x - min(x))/diff(range(x)))
#'
#' #---- Simulate MM matrix
#' w_ord <- c(.5, .35, .15) # Weight of each neighbours orders
#' ord <- length(w_ord) - 1 # Order of neighbours to include
#' H_sel_sim <- sim_MM_matrix(
#'   W = W_sel, m = m, ord = ord, w_ord = w_ord, id_vec = rep(1, nrow(W_sel))
#' )
#'
#' #---- Simulate outcomes
#' ## Linear term parameters
#' gamma <- -.5 # Intercept
#' beta <- c(1, .5) # Covariates coefficients
#' ## CAR random effects
#' phi_car <- sim_car(W = W_sel, alpha = .9, tau = 5)
#' # Areal log relative risks
#' l_RR <- X_cent %*% beta + phi_car
#' ## Membership log relative risks
#' l_RR_mm <- as.numeric(apply(H_sel_sim, 1, function(x) x %*% l_RR))
#' ## Expected rates
#' exp_rates <- rpois(m, lambda = 20)
#' ## Outcomes
#' y <- rpois(m, lambda = exp_rates*exp(l_RR_mm))
#'
#' #---- Create dataset for stan function
#' d_sel <- list(
#'   # Number of areas
#'   n = nrow(W_sel),
#'   # Covariates
#'   k = ncol(X_cent),
#'   X_cov = X_cent,
#'   # Adjacency
#'   W_n = sum(W_sel) / 2,
#'   # Number of neighbour pairs
#'   W = W_sel,
#'   # Memberships
#'   m = nrow(H_sel_sim),
#'   H = H_sel_sim,
#'   # Outcomes
#'   y = y,
#'   log_offset = log(exp_rates),
#'   # Prior parameters
#'   ## Intercept (mean and sd of normal prior)
#'   mu_gamma = 0, sigma_gamma = 1,
#'   ## Covariates (mean and sd of normal prior)
#'   mu_beta = 0, sigma_beta = 1,
#'   ## Marginal precision gamma prior
#'   tau_shape = 2,
#'   tau_rate = 0.2
#' )
#'
#' #---- HMC parameters
#' niter <- 1E4
#' nchains <- 4
#'
#' #---- Stan sampling
#' fit <- car_mm(
#'   d_list = d_sel,
#'   # arguments passed to sampling
#'   iter = niter, chains = nchains, refresh = 500,
#'   control = list(adapt_delta = .99, max_treedepth = 15)
#' )
#'
#'}
car_mm <- function(d_list, ...) {
  out <- rstan::sampling(stanmodels$CARMM_COV_P, data = d_list, ...)
  return(out)
}


