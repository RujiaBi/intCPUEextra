#' Get observation-level predicted values and residuals
#'
#' Extracts observation-level fitted values and residual diagnostics from a
#' fitted `intCPUE` model.
#'
#' `intCPUE` is a delta model. Response-scale fitted values are unconditional
#' means, while log-scale diagnostics are defined only for positive observations
#' and use the conditional positive-CPUE mean.
#'
#' @param object An object of class `intCPUE` returned by [intCPUE()].
#' @param data Optional original data.frame to bind to the returned fitted
#'   values. If supplied, it must have the same number of rows and order as the
#'   fitted observations.
#'
#' @return A data.frame with fitted values and residuals.
#' @export
get_predicted <- function(object, data = NULL) {
  if (!inherits(object, "intCPUE")) {
    stop("`object` must be an `intCPUE` fit from intCPUE().", call. = FALSE)
  }

  par_best <- object$obj$env$last.par.best
  if (is.null(par_best)) {
    stop("Could not find `last.par.best` in fitted object.", call. = FALSE)
  }

  rep_obj <- object$obj$report(par_best)
  req <- c(
    "eta_hat_encounter_i",
    "eta_hat_positive_i",
    "encounter_prob_i",
    "log_positive_mean_i",
    "mu_hat_i"
  )
  miss <- setdiff(req, names(rep_obj))
  if (length(miss)) {
    stop(
      "Missing fitted-value outputs in `obj$report()`: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  observed <- object$data_tmb$b_i
  observed_log <- ifelse(observed > 0, log(observed), NA_real_)
  predicted_log <- as.numeric(rep_obj$log_positive_mean_i)

  out <- data.frame(
    obs = seq_along(rep_obj$mu_hat_i),
    observed = observed,
    encounter = object$data_tmb$e_i,
    tid = object$data_tmb$t_i,
    vesid = object$data_tmb$v_i,
    flagid = object$data_tmb$f_i,
    eta_hat_encounter = rep_obj$eta_hat_encounter_i,
    eta_hat_positive = rep_obj$eta_hat_positive_i,
    encounter_prob = rep_obj$encounter_prob_i,
    positive_mean = exp(rep_obj$log_positive_mean_i),
    fitted = rep_obj$mu_hat_i,
    observed_log = observed_log,
    predicted_log = predicted_log,
    residual_raw = observed - rep_obj$mu_hat_i,
    residual_log = observed_log - predicted_log
  )

  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      stop("`data` must be a data.frame if supplied.", call. = FALSE)
    }
    if (nrow(data) != nrow(out)) {
      stop(
        "`data` must have the same number of rows as the fitted model. Got ",
        nrow(data), " rows but expected ", nrow(out), ".",
        call. = FALSE
      )
    }
    out <- cbind(data, out)
  }

  out
}

#' Get observation-level fitted values and residuals
#'
#' Backward-compatible wrapper for [get_predicted()].
#'
#' @inheritParams get_predicted
#' @return A data.frame with predicted values and residuals.
#' @export
get_fitted <- function(object, data = NULL) {
  get_predicted(object = object, data = data)
}
