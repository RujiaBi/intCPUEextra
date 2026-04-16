#' Plot standardized CPUE index and CV
#'
#' Creates a two-panel plot showing the standardized CPUE index and its
#' coefficient of variation (CV) over time.
#'
#' @param index_df A data.frame returned by [get_index()] containing at least
#'   `time`, `index`, and `cv`. If `areaid` is present, each area is drawn as a
#'   separate colored series.
#' @param time_values Optional vector of labels for the x-axis. It may have the
#'   same length as `nrow(index_df)` or the number of unique values in
#'   `index_df$time`.
#' @param time_positions Optional numeric vector of x-axis positions. It may
#'   have the same length as `nrow(index_df)` or the number of unique values in
#'   `index_df$time`.
#' @param x_text_angle Optional numeric angle for x-axis text. If `NULL`,
#'   `plot_index()` uses `90` for year-month-like labels and `45` otherwise.
#' @param x_text_hjust Optional horizontal justification for x-axis text. If
#'   `NULL`, a value matching `x_text_angle` is chosen automatically.
#' @param x_text_vjust Optional vertical justification for x-axis text. If
#'   `NULL`, a value matching `x_text_angle` is chosen automatically.
#'
#' @return A ggplot object.
#' @export
plot_index <- function(
    index_df,
    time_values = NULL,
    time_positions = NULL,
    x_text_angle = NULL,
    x_text_hjust = NULL,
    x_text_vjust = NULL
) {
  req <- c("time", "index", "cv")
  miss <- setdiff(req, names(index_df))
  if (length(miss)) {
    stop(
      "`index_df` is missing required columns: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  has_area <- "areaid" %in% names(index_df)
  area_vals <- if (has_area) as.character(index_df$areaid) else rep("All data", nrow(index_df))
  time_levels <- sort(unique(index_df$time))

  if (!is.null(time_positions)) {
    if (!is.numeric(time_positions)) {
      stop("`time_positions` must be numeric.", call. = FALSE)
    }
    if (length(time_positions) == length(time_levels)) {
      time_seq <- unname(time_positions[match(index_df$time, time_levels)])
    } else if (length(time_positions) == nrow(index_df)) {
      time_seq <- time_positions
    } else {
      stop(
        "`time_positions` must have length `nrow(index_df)` or the number of unique time values.",
        call. = FALSE
      )
    }
  } else {
    time_seq <- index_df$time
  }

  if (!is.null(time_values)) {
    if (length(time_values) == length(time_levels)) {
      x_breaks <- if (is.null(time_positions)) time_levels else time_positions
      x_labels <- time_values
    } else if (length(time_values) == nrow(index_df)) {
      x_breaks <- time_seq
      x_labels <- time_values
    } else {
      stop(
        "`time_values` must have length `nrow(index_df)` or the number of unique time values.",
        call. = FALSE
      )
    }
  } else {
    x_breaks <- sort(unique(time_seq))
    x_labels <- x_breaks
  }

  plot_dat <- rbind(
    data.frame(
      time = time_seq,
      value = index_df$index,
      panel = "Standardized CPUE",
      areaid = area_vals
    ),
    data.frame(
      time = time_seq,
      value = index_df$cv,
      panel = "CV",
      areaid = area_vals
    )
  )

  plot_dat$panel <- factor(
    plot_dat$panel,
    levels = c("Standardized CPUE", "CV")
  )

  if (has_area) {
    p <- ggplot2::ggplot(
      plot_dat,
      ggplot2::aes(x = time, y = value, colour = areaid, group = areaid)
    ) +
      ggplot2::geom_line(linewidth = 1.05, lineend = "round") +
      ggplot2::geom_point(size = 2.1)
  } else {
    p <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = time, y = value)) +
      ggplot2::geom_line(
        data = subset(plot_dat, panel == "Standardized CPUE"),
        colour = "#0F766E",
        linewidth = 1.1,
        lineend = "round"
      ) +
      ggplot2::geom_point(
        data = subset(plot_dat, panel == "Standardized CPUE"),
        colour = "#0F766E",
        fill = "white",
        shape = 21,
        stroke = 0.7,
        size = 2.2
      ) +
      ggplot2::geom_line(
        data = subset(plot_dat, panel == "CV"),
        colour = "#B45309",
        linewidth = 1,
        lineend = "round"
      ) +
      ggplot2::geom_point(
        data = subset(plot_dat, panel == "CV"),
        colour = "#B45309",
        fill = "white",
        shape = 21,
        stroke = 0.7,
        size = 1.9
      )
  }

  p <- p +
    ggplot2::facet_wrap(~ panel, ncol = 1, scales = "free_y", strip.position = "top") +
    ggplot2::labs(x = "Time", y = NULL, colour = if (has_area) "Area" else NULL) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey96", colour = "grey80"),
      strip.text.x = ggplot2::element_text(face = "bold"),
      strip.placement = "outside",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "grey92"),
      panel.border = ggplot2::element_rect(colour = "grey40"),
      axis.title.x = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(colour = "grey20"),
      axis.text.y = ggplot2::element_text(colour = "grey20"),
      legend.position = if (has_area) "top" else "none"
    )

  if (!is.null(time_values) || !is.null(time_positions)) {
    x_labels_chr <- as.character(x_labels)
    has_month_labels <- any(grepl("-", x_labels_chr, fixed = TRUE))

    if (is.null(x_text_angle)) {
      x_text_angle <- if (has_month_labels || any(nchar(x_labels_chr) > 4L)) 90 else 45
    }

    angle_norm <- ((x_text_angle %% 360) + 360) %% 360
    if (is.null(x_text_hjust)) {
      x_text_hjust <- 1
    }
    if (is.null(x_text_vjust)) {
      x_text_vjust <- if (angle_norm %in% c(90, 270)) 0.5 else 1
    }

    p <- p +
      ggplot2::scale_x_continuous(
        breaks = x_breaks,
        labels = x_labels
      ) +
      ggplot2::labs(x = if (has_month_labels) "Year-Month" else "Year") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = x_text_angle,
          hjust = x_text_hjust,
          vjust = x_text_vjust,
          colour = "grey20"
        )
      )
  }

  p
}
