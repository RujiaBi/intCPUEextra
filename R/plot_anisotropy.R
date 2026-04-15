#' Plot anisotropy ellipses for the fitted spatial fields
#'
#' Visualizes the anisotropy implied by `ln_H_input` together with the fitted
#' range parameters for the two linear predictors.
#'
#' @param object An object of class `intCPUE` returned by [intCPUE()], or a TMB
#'   object list containing `obj$env$last.par.best`.
#' @param colors Character vector of length 2 giving ellipse colors for the
#'   first and second linear predictors.
#' @param labels Character vector of length 2 used in the legend.
#' @param n_points Number of points used to draw each ellipse.
#'
#' @return A ggplot object.
#' @export
plot_anisotropy <- function(
    object,
    colors = c("#2E7D32", "#111111"),
    labels = c("1st linear predictor", "2nd linear predictor"),
    n_points = 361L
) {
  par_best <- NULL
  rep_obj <- NULL

  if (inherits(object, "intCPUE")) {
    par_best <- object$obj$env$last.par.best
  } else if (is.list(object) &&
             !is.null(object$obj) &&
             !is.null(object$obj$env) &&
             !is.null(object$obj$env$last.par.best)) {
    par_best <- object$obj$env$last.par.best
  }

  if (is.null(par_best)) {
    stop(
      "`object` must be an `intCPUE` fit or contain `obj$env$last.par.best`.",
      call. = FALSE
    )
  }

  if ((inherits(object, "intCPUE") || (is.list(object) && !is.null(object$obj)))) {
    rep_obj <- try(object$obj$report(par_best), silent = TRUE)
  }

  ln_H_input <- as.numeric(par_best[grep("^ln_H_input", names(par_best))])
  if (length(ln_H_input) != 2L) {
    stop("Could not find the two anisotropy parameters `ln_H_input`.", call. = FALSE)
  }

  H <- matrix(
    c(
      exp(ln_H_input[1]), ln_H_input[2],
      ln_H_input[2], (1 + ln_H_input[2]^2) / exp(ln_H_input[1])
    ),
    nrow = 2,
    byrow = TRUE
  )

  eig <- eigen(H)

  extract_range <- function(idx) {
    nm_range <- paste0("ln_range_", idx)
    nm_kappa <- paste0("ln_kappa_", idx)

    if (nm_range %in% names(par_best)) {
      return(exp(unname(par_best[[nm_range]])))
    }
    if (nm_kappa %in% names(par_best)) {
      return(sqrt(8.0) / exp(unname(par_best[[nm_kappa]])))
    }

    stop(
      "Could not find `", nm_range, "` or `", nm_kappa, "` in `last.par.best`.",
      call. = FALSE
    )
  }

  ranges <- c(extract_range(1), extract_range(2))

  if (length(colors) != 2L) {
    stop("`colors` must have length 2.", call. = FALSE)
  }
  if (length(labels) != 2L) {
    stop("`labels` must have length 2.", call. = FALSE)
  }
  n_points <- as.integer(n_points)
  if (is.na(n_points) || n_points < 20L) {
    stop("`n_points` must be an integer >= 20.", call. = FALSE)
  }

  make_axes <- function(range_val, label) {
    major <- eig$vectors[, 1] * eig$values[1] * range_val
    minor <- eig$vectors[, 2] * eig$values[2] * range_val
    list(major = major, minor = minor, label = label)
  }

  axis_list <- list(
    make_axes(ranges[1], labels[1]),
    make_axes(ranges[2], labels[2])
  )

  theta <- seq(0, 2 * pi, length.out = n_points)
  ellipse_df <- do.call(
    rbind,
    lapply(seq_along(axis_list), function(i) {
      ax <- axis_list[[i]]
      xy <- vapply(
        theta,
        function(th) ax$major * cos(th) + ax$minor * sin(th),
        numeric(2)
      )
      data.frame(
        x = xy[1, ],
        y = xy[2, ],
        component = ax$label
      )
    })
  )

  axis_df <- do.call(
    rbind,
    lapply(seq_along(axis_list), function(i) {
      ax <- axis_list[[i]]
      data.frame(
        x = c(-ax$major[1], ax$major[1], -ax$minor[1], ax$minor[1]),
        y = c(-ax$major[2], ax$major[2], -ax$minor[2], ax$minor[2]),
        xend = c(ax$major[1], -ax$major[1], ax$minor[1], -ax$minor[1]),
        yend = c(ax$major[2], -ax$major[2], ax$minor[2], -ax$minor[2]),
        component = ax$label,
        axis_type = rep(c("major", "major", "minor", "minor"), each = 1L)
      )
    })
  )

  lim <- max(abs(c(ellipse_df$x, ellipse_df$y)), na.rm = TRUE) * 1.1

  ggplot2::ggplot() +
    ggplot2::geom_path(
      data = ellipse_df,
      ggplot2::aes(x = x, y = y, colour = component),
      linewidth = 1
    ) +
    ggplot2::geom_segment(
      data = axis_df,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, colour = component),
      linewidth = 0.4,
      alpha = 0.7
    ) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.4) +
    ggplot2::coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    ggplot2::scale_colour_manual(values = stats::setNames(colors, labels)) +
    ggplot2::labs(
      x = "Scaled Easting",
      y = "Scaled Northing",
      colour = NULL,
      title = "Distance at 10% correlation"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top"
    )
}
