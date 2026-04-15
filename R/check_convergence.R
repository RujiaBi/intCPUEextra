#' Check Model Convergence
#'
#' Summarizes basic optimizer diagnostics for a fitted `intCPUE` model.
#'
#' @param object A fitted `intCPUE` object.
#'
#' @return A one-row data.frame containing the optimizer convergence code,
#'   optimizer message, and maximum absolute gradient component at the reported
#'   optimum.
#' @export
check_convergence <- function(object) {
  if (!inherits(object, "intCPUE")) {
    stop("`object` must be an `intCPUE` fit from intCPUE().", call. = FALSE)
  }

  data.frame(
    convergence = object$opt$convergence,
    message = object$opt$message,
    max_grad = if (!is.null(object$opt$max_grad)) {
      object$opt$max_grad
    } else {
      max(abs(object$obj$gr(object$opt$par)))
    }
  )
}
