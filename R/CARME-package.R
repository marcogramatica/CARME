#' The 'CARME' package.
#'
#' @description CAR-MM modelling in Stan
#'
#' @docType package
#' @name CARME-package
#' @aliases CARME
#' @useDynLib CARME, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (2023). RStan: the R interface to Stan. R package version 
#' 2.26.11. https://mc-stan.org
#' 
#'  Marco Gramatica. Silvia Liverani. Peter Congdon. 
#'  Structure Induced by a Multiple Membership Transformation on the Conditional
#'  Autoregressive Model. Bayesian Analysis Advance Publication 1 - 25, 2023.
#'  https://doi.org/10.1214/23-BA1370
#'  
#'  Petrof, O, Neyens, T, Nuyts, V, Nackaerts, K, Nemery, B, Faes, C. On the 
#'  impact of residential history in the spatial analysis of diseases with a 
#'  long latency period: A study of mesothelioma in Belgium. 
#'  Statistics in Medicine. 2020; 39: 3840– 3866.
#'  https://doi.org/10.1002/sim.8697
#'  
#'  Marco Gramatica, Peter Congdon, Silvia Liverani, Bayesian Modelling for 
#'  Spatially Misaligned Health Areal Data: A Multiple Membership Approach, 
#'  Journal of the Royal Statistical Society Series C: Applied Statistics, 
#'  Volume 70, Issue 3, June 2021, Pages 645–666,
#'  https://doi.org/10.1111/rssc.12480
#'
NULL
