# Internal: build TMB map to switch model components on/off
#
# The map controls which parameters are estimated vs fixed.
# In TMB, map entries set to NA are fixed at their initial value.
#
# Important statistical principle:
# We automatically prevent estimation of variance parameters that are not identifiable
# (e.g., flag temporal random effects with < 2 time points).
#
.build_flag_t_constraint <- function(has_tf, n_f = ncol(has_tf) + 1L) {
  if (is.null(has_tf)) {
    stop("`has_tf` is required to build temporal flag constraints.", call. = FALSE)
  }
  if (!is.matrix(has_tf)) {
    stop("`has_tf` must be a matrix.", call. = FALSE)
  }
  if (length(n_f) != 1L || is.na(n_f) || n_f < 1L) {
    stop("`n_f` must be a positive integer.", call. = FALSE)
  }

  keep_tf <- has_tf > 0
  n_flag_cols <- max(0L, n_f - 1L)
  if (!identical(ncol(keep_tf), n_flag_cols)) {
    stop("`has_tf` has incompatible column count for the supplied `n_f`.", call. = FALSE)
  }

  estimable_flag <- colSums(keep_tf) >= 2L
  flag_t_index <- matrix(-1L, nrow = nrow(keep_tf), ncol = ncol(keep_tf))
  next_idx <- 0L

  if (ncol(keep_tf) > 0L) {
    for (j in seq_len(ncol(keep_tf))) {
      if (!estimable_flag[j]) {
        keep_tf[, j] <- FALSE
        next
      }
      rows_j <- which(keep_tf[, j])
      # Use the first observed time as the reference level and fix it to zero.
      anchor_row <- rows_j[1L]
      keep_tf[anchor_row, j] <- FALSE

      free_rows <- rows_j[-1L]
      if (length(free_rows) > 0L) {
        idx_j <- seq.int(from = next_idx, length.out = length(free_rows))
        flag_t_index[free_rows, j] <- as.integer(idx_j)
        next_idx <- next_idx + length(free_rows)
      }
    }
  }

  list(
    estimable_flag = estimable_flag,
    any_estimable = any(estimable_flag),
    flag_t_index = flag_t_index,
    n_free = next_idx
  )
}

.make_map_intCPUE <- function(parameters, n_f,
                              obs_sd = c("shared","flag"),
                              pop_spatiotemporal_type = c("rw", "ar1"),
                              q_diffs_time = c("on","off"),
                              has_tf = NULL,
                              estimable_flag = NULL) {
  q_diffs_time <- match.arg(q_diffs_time)
  obs_sd <- match.arg(obs_sd)
  pop_spatiotemporal_type <- match.arg(pop_spatiotemporal_type)
  
  map <- list()
  has_param <- function(name) !is.null(parameters[[name]])

  # ---------------------------------------------------------
  # Observation SD
  # ---------------------------------------------------------
  if (obs_sd == "shared") {
    if (has_param("ln_sd_flag")) {
      map$ln_sd_flag <- factor(rep(NA, length(parameters$ln_sd_flag)))
    }
  } else {
    if (has_param("ln_sd")) {
      map$ln_sd <- factor(rep(NA, length(parameters$ln_sd)))
    }
  }

  # Fixed core architecture:
  # spatial + spatiotemporal + vessel + q_diffs_system are always on.
  # q_diffs_time / q_diffs_spatial are handled by separate TMB templates, not map.
  # When the spatiotemporal process is RW, AR1 correlation parameters are unused
  # and should be fixed at their initial value.
  if (pop_spatiotemporal_type == "rw") {
    if (has_param("transf_rho_1")) {
      map$transf_rho_1 <- factor(rep(NA, length(parameters$transf_rho_1)))
    }
    if (has_param("transf_rho_2")) {
      map$transf_rho_2 <- factor(rep(NA, length(parameters$transf_rho_2)))
    }
  }

  if (has_param("flag_t_ln_std_dev_1")) {
    if (q_diffs_time == "off" || n_f <= 1L) {
      map$flag_t_ln_std_dev_1 <- factor(rep(NA, length(parameters$flag_t_ln_std_dev_1)))
      if (has_param("flag_t_ln_std_dev_2")) {
        map$flag_t_ln_std_dev_2 <- factor(rep(NA, length(parameters$flag_t_ln_std_dev_2)))
      }
    } else {
      if (is.null(has_tf)) {
        stop("q_diffs_time='on' requires `has_tf` to construct temporal flag effects.", call. = FALSE)
      }
      if (is.null(estimable_flag)) {
        tf_constraint <- .build_flag_t_constraint(has_tf, n_f = n_f)
        estimable_flag <- tf_constraint$estimable_flag
      }
      if (!any(estimable_flag)) {
        map$flag_t_ln_std_dev_1 <- factor(rep(NA, length(parameters$flag_t_ln_std_dev_1)))
        if (has_param("flag_t_ln_std_dev_2")) {
          map$flag_t_ln_std_dev_2 <- factor(rep(NA, length(parameters$flag_t_ln_std_dev_2)))
        }
      }
    }
  }

  # If there is only one flag in the data, the flag-system terms are not identifiable.
  if (n_f <= 1L) {
    if (has_param("flag_f_1")) map$flag_f_1 <- factor(rep(NA, length(parameters$flag_f_1)))
    if (has_param("flag_f_2")) map$flag_f_2 <- factor(rep(NA, length(parameters$flag_f_2)))
    if (has_param("flag_ln_std_dev_1")) map$flag_ln_std_dev_1 <- factor(rep(NA, length(parameters$flag_ln_std_dev_1)))
    if (has_param("flag_ln_std_dev_2")) map$flag_ln_std_dev_2 <- factor(rep(NA, length(parameters$flag_ln_std_dev_2)))
  }

  map
}
