#' Simulated Information IV Dataset
#'
#' A simulated dataset for demonstrating information IV methods. Contains 
#' experimental data with randomized information treatments and Bayesian updating.
#'
#' @format A data.table with 500 observations and 9 variables:
#' \describe{
#'   \item{tau}{Individual-specific treatment effect parameters (normalized to mean 1)}
#'   \item{alpha}{Learning rates, negatively correlated with tau}
#'   \item{Z}{Binary treatment assignment (0/1)}
#'   \item{signal}{Information signal (1 for high signal, -1 for low signal)}
#'   \item{prior}{Prior beliefs (affected by unobserved confounder V)}
#'   \item{posterior}{Posterior beliefs after updating}
#'   \item{Y}{Outcome variable}
#'   \item{Y0}{Counterfactual outcome (pre-treatment)}
#'   \item{dX}{Change in beliefs (posterior - prior)}
#'   \item{dY}{Change in outcomes (Y - Y0)}
#' }
#' @details 
#' This dataset is generated with endogeneity where unobserved factors (V, U) 
#' affect both beliefs and outcomes. The true average partial effect (APE) is 1.
#' People with higher belief effects (tau) have lower learning rates (alpha).
#' 
#' @source Simulated data 
#' @examples
#' data(info.sim)
#' head(info.sim)
#' 
#' # Check the structure
#' str(info.sim)
#' 
#' # Summary statistics
#' summary(info.sim)
"info.sim"
