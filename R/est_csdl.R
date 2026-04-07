#' CS-DL Estimator Internals
#'
#' Internal functions for the Cross-Sectionally augmented Distributed Lag
#' (CS-DL) estimator of Chudik et al. (2016). Long-run coefficient `w2`
#' is estimated directly as the coefficient on the level of x.
#'
#' @name csdl_estimator
#' @keywords internal
NULL

# CS-DL implementation is handled within dcce() by adding first differences
# of x as regressors alongside the level. The long-run coefficient is directly
# estimated as the coefficient on the level of x.
# This file provides helper documentation. The actual estimation logic will
# be integrated into dcce() when model="csdl" is specified.
