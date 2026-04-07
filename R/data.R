#' Simulated Dynamic Panel Dataset
#'
#' A synthetic panel dataset generated from the dynamic common correlated
#' effects DGP of Chudik and Pesaran (2015), equation (1). Generated with
#' N = 30 cross-sectional units and T = 50 time periods. The true mean
#' group parameters are: autoregressive coefficient 0.50, slope on x 1.00.
#' Factor structure uses a single AR(1) common factor with persistence 0.60.
#'
#' The true parameter values are stored in the companion object
#' `dcce_sim_truth` and can be used in tests to verify estimator
#' consistency.
#'
#' @format A data frame with 1500 rows and 4 columns:
#' \describe{
#'   \item{unit}{Cross-sectional unit identifier (integer, 1--30)}
#'   \item{time}{Time period (integer, 1--50)}
#'   \item{y}{Simulated dependent variable}
#'   \item{x}{Simulated regressor (cross-sectionally dependent via common
#'     factor)}
#' }
#' @references
#' Chudik, A. and Pesaran, M. H. (2015). Common correlated effects
#' estimation of heterogeneous dynamic panel data models with weakly
#' exogenous regressors. \emph{Journal of Econometrics}, 188(2), 393--420.
#' @seealso [dcce_sim_truth] for the true parameter values.
"dcce_sim"

#' True Parameters for the Simulated Panel Dataset
#'
#' A named list containing the true mean group parameter values used to
#' generate [dcce_sim]. Use in tests to verify that [dcce()] recovers
#' parameters close to their true values.
#'
#' @format A named list with elements:
#' \describe{
#'   \item{beta1_mg}{True mean group autoregressive coefficient (~0.50)}
#'   \item{beta2_mg}{True mean group slope on x (~1.00)}
#'   \item{N}{Number of cross-sectional units (30)}
#'   \item{T}{Number of time periods (50)}
#'   \item{seed}{Random seed used for generation (20240101)}
#' }
"dcce_sim_truth"

#' Penn World Tables Growth Panel Dataset
#'
#' Panel dataset from Jan Ditzen's xtdcce2 Stata package, originally
#' derived from the Penn World Tables 8. Used in Ditzen (2018, Stata
#' Journal 18:3, 585--617) to illustrate MG, CCE, and DCCE estimation.
#' Contains 93 countries observed from 1960 to 2007 (balanced).
#'
#' @format A data frame with 4464 rows and 8 columns:
#' \describe{
#'   \item{id}{Country numeric identifier}
#'   \item{year}{Year (1960--2007)}
#'   \item{log_rgdpo}{Log real GDP (output-side)}
#'   \item{log_hc}{Log human capital index}
#'   \item{log_ck}{Log physical capital stock}
#'   \item{log_ngd}{Log of population growth plus depreciation}
#'   \item{country}{Country identifier (character)}
#'   \item{d_log_rgdpo}{First difference of log_rgdpo (GDP growth)}
#' }
#' @source
#' Downloaded from
#' \url{https://github.com/JanDitzen/xtdcce2}.
#' @references
#' Ditzen, J. (2018). Estimating dynamic common-correlated effects in
#' Stata. \emph{The Stata Journal}, 18(3), 585--617.
"pwt8"

# dcce_debt dataset to be added in a future release.
# See data-raw/prepare_dcce_debt.R (to be created).
