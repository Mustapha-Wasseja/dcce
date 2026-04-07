#' marginaleffects Compatibility for dcce_fit Objects
#'
#' Provides S3 methods on the \pkg{marginaleffects} internal generics
#' \code{get_coef}, \code{get_vcov}, \code{get_predict}, and
#' \code{find_predictors} so that
#' \code{marginaleffects::avg_slopes()},
#' \code{marginaleffects::avg_predictions()}, and
#' \code{marginaleffects::hypotheses()} work directly on
#' \code{dcce_fit} objects.
#'
#' The methods are registered dynamically in \code{.onLoad()} if the
#' \pkg{marginaleffects} package is available; this keeps the dependency
#' in \code{Suggests} rather than \code{Imports}.
#'
#' @name marginaleffects_compat
#' @keywords internal
NULL


#' @keywords internal
get_coef.dcce_fit <- function(model, ...) {
  stats::coef(model)
}

#' @keywords internal
get_vcov.dcce_fit <- function(model, ...) {
  stats::vcov(model)
}

#' @keywords internal
get_predict.dcce_fit <- function(model, newdata = NULL, ...) {
  if (is.null(newdata)) {
    preds <- as.numeric(stats::fitted(model))
    if (length(preds) == 0L) {
      preds <- as.numeric(stats::predict(model))
    }
    data.frame(rowid = seq_along(preds), estimate = preds)
  } else {
    preds <- as.numeric(stats::predict(model, newdata = newdata))
    data.frame(rowid = seq_len(nrow(newdata)), estimate = preds)
  }
}

#' @keywords internal
find_predictors.dcce_fit <- function(model, ...) {
  list(conditional = model$x_names)
}

#' @keywords internal
find_response.dcce_fit <- function(model, ...) {
  model$y_name
}


#' Register marginaleffects S3 methods at package load
#'
#' Called from \code{.onLoad()} in \code{R/zzz.R}. Dynamically registers
#' methods on the \pkg{marginaleffects} generics if that package is
#' installed. Safe no-op otherwise.
#' @keywords internal
.register_marginaleffects_methods <- function() {
  if (!requireNamespace("marginaleffects", quietly = TRUE)) return(invisible())

  # Only register if the generic is actually exported (defensive)
  me_ns <- asNamespace("marginaleffects")

  maybe_register <- function(generic, class, method) {
    if (exists(generic, envir = me_ns, inherits = FALSE)) {
      registerS3method(generic, class, method, envir = me_ns)
    }
  }

  maybe_register("get_coef",        "dcce_fit", get_coef.dcce_fit)
  maybe_register("get_vcov",        "dcce_fit", get_vcov.dcce_fit)
  maybe_register("get_predict",     "dcce_fit", get_predict.dcce_fit)
  maybe_register("find_predictors", "dcce_fit", find_predictors.dcce_fit)
  maybe_register("find_response",   "dcce_fit", find_response.dcce_fit)

  invisible()
}
