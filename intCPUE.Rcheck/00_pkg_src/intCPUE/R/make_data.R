.resolve_intCPUE_formulas <- function(formula = NULL,
                                      formula_catchability = NULL,
                                      formula_population = NULL,
                                      caller = "make_data") {
  if (!is.null(formula) && !is.null(formula_catchability)) {
    stop(
      sprintf(
        "In %s(), supply either `formula` or `formula_catchability`, not both.",
        caller
      ),
      call. = FALSE
    )
  }

  if (!is.null(formula) && is.null(formula_catchability)) {
    formula_catchability <- formula
  }

  list(
    formula_catchability = formula_catchability,
    formula_population = formula_population
  )
}

.normalize_smoother_output <- function(sm, n_rows) {
  if (!isTRUE(sm$has_smooths)) {
    return(list(
      has_smooths = FALSE,
      Xs = matrix(0, nrow = n_rows, ncol = 0L),
      Zs = list(),
      basis_out = list(),
      labels = list(),
      classes = list(),
      sm_dims = integer(0),
      b_smooth_start = integer(0)
    ))
  }

  Xs <- sm$Xs
  Zs <- sm$Zs
  sm_dims <- as.integer(sm$sm_dims)
  b_smooth_start <- as.integer(sm$b_smooth_start)

  if (!identical(nrow(Xs), n_rows)) {
    stop("parse_smoothers() returned Xs with nrow != expected nrow.", call. = FALSE)
  }

  if (length(Zs)) {
    for (k in seq_along(Zs)) {
      if (!identical(nrow(Zs[[k]]), n_rows)) {
        stop("parse_smoothers() returned Zs[[k]] with nrow != expected nrow.", call. = FALSE)
      }
      if (!inherits(Zs[[k]], "sparseMatrix")) {
        Zs[[k]] <- Matrix::Matrix(Zs[[k]], sparse = TRUE)
      }
    }
  }

  list(
    has_smooths = TRUE,
    Xs = Xs,
    Zs = Zs,
    basis_out = sm$basis_out,
    labels = sm$labels,
    classes = sm$classes,
    sm_dims = sm_dims,
    b_smooth_start = b_smooth_start
  )
}

.smooth_vars_from_basis <- function(basis_out) {
  .keep_var_names <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x) & x != "NA"]
    x
  }

  needed <- character(0)

  if (!length(basis_out)) {
    return(needed)
  }

  for (basis_i in basis_out) {
    if (!length(basis_i)) next
    for (sm in basis_i) {
      term_i <- sm$term
      if (length(term_i)) {
        term_i <- .keep_var_names(term_i)
        needed <- c(needed, term_i)
      }
      if (!is.null(sm$by) && is.character(sm$by) && length(sm$by) == 1L &&
          length(.keep_var_names(sm$by)) == 1L) {
        needed <- c(needed, sm$by)
      }
    }
  }

  unique(.keep_var_names(needed))
}

.expand_projection_over_time <- function(projection_data, tid_values) {
  projection_data <- as.data.frame(projection_data)
  tid_values <- as.integer(tid_values)

  out <- projection_data[rep(seq_len(nrow(projection_data)), times = length(tid_values)), , drop = FALSE]
  out$tid <- rep(tid_values, each = nrow(projection_data))
  rownames(out) <- NULL
  out
}

.warn_projection_na <- function(projection_data, needed_vars) {
  if (!length(needed_vars)) {
    return(invisible(NULL))
  }

  na_counts <- vapply(
    needed_vars,
    function(nm) sum(is.na(projection_data[[nm]])),
    integer(1)
  )
  na_counts <- na_counts[na_counts > 0L]

  if (!length(na_counts)) {
    return(invisible(NULL))
  }

  warning(
    paste0(
      "`projection_data` contains NA in population covariates: ",
      paste(sprintf("%s (%d)", names(na_counts), as.integer(na_counts)), collapse = ", "),
      ". For projection rows with NA, the affected population smooth term(s) will contribute 0."
    ),
    call. = FALSE
  )

  invisible(NULL)
}

.align_projection_data_to_key <- function(projection_data, key, needed_vars, tid_values) {
  projection_data <- as.data.frame(projection_data)
  tid_values <- as.integer(tid_values)

  coord_cols <- c("utm_x_scale", "utm_y_scale")
  has_tid <- "tid" %in% names(projection_data)
  req <- unique(c(coord_cols, if (has_tid) "tid", needed_vars))
  miss <- setdiff(req, names(projection_data))
  if (length(miss)) {
    stop(
      "`projection_data` is missing required columns: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  if (has_tid) {
    if (anyNA(projection_data$tid)) {
      stop("`projection_data$tid` must not contain NA.", call. = FALSE)
    }

    tid_num <- suppressWarnings(as.numeric(as.character(projection_data$tid)))
    if (any(!is.finite(tid_num)) || any(abs(tid_num - round(tid_num)) > 1e-9)) {
      stop("`projection_data$tid` must be integer-valued and use the same coding as `data_utm$tid`.", call. = FALSE)
    }
    projection_data$tid <- as.integer(round(tid_num))

    bad_tid <- setdiff(unique(projection_data$tid), tid_values)
    if (length(bad_tid)) {
      stop(
        "`projection_data$tid` contains values not used by the model: ",
        paste(sort(bad_tid), collapse = ", "),
        call. = FALSE
      )
    }

    proj_id <- paste(projection_data$utm_x_scale, projection_data$utm_y_scale, projection_data$tid, sep = "\r")
    key_gt <- .expand_projection_over_time(key[, coord_cols, drop = FALSE], tid_values)
    key_id <- paste(key_gt$utm_x_scale, key_gt$utm_y_scale, key_gt$tid, sep = "\r")

    if (anyDuplicated(proj_id)) {
      stop("`projection_data` must contain at most one row per extrapolation cell-time combination.", call. = FALSE)
    }

    idx <- match(key_id, proj_id)
    if (anyNA(idx)) {
      stop(
        "`projection_data` must contain one row for every extrapolation grid cell-time combination ",
        "used by the model. Missing combinations detected.",
        call. = FALSE
      )
    }

    out <- projection_data[idx, req, drop = FALSE]
  } else {
    proj_id <- paste(projection_data$utm_x_scale, projection_data$utm_y_scale, sep = "\r")
    key_id <- paste(key$utm_x_scale, key$utm_y_scale, sep = "\r")

    if (anyDuplicated(proj_id)) {
      stop("`projection_data` must contain at most one row per extrapolation cell.", call. = FALSE)
    }

    idx <- match(key_id, proj_id)
    if (anyNA(idx)) {
      stop(
        "`projection_data` must contain one row for every extrapolation grid cell ",
        "used by the model. Missing cells detected.",
        call. = FALSE
      )
    }

    out <- .expand_projection_over_time(
      projection_data[idx, unique(c(coord_cols, needed_vars)), drop = FALSE],
      tid_values = tid_values
    )
  }

  rownames(out) <- NULL
  out
}

.default_population_projection_data <- function(data_utm, key, basis_out, tid_values) {
  needed_vars <- .smooth_vars_from_basis(basis_out)
  coord_cols <- c("utm_x_scale", "utm_y_scale")

  if (!length(needed_vars)) {
    out <- .expand_projection_over_time(key[, coord_cols, drop = FALSE], tid_values = tid_values)
    rownames(out) <- NULL
    return(out)
  }

  miss <- setdiff(needed_vars, names(data_utm))
  if (length(miss)) {
    stop(
      "Population smooth terms require covariates not found in `data_utm`: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  cell_id <- paste(data_utm$utm_x_scale, data_utm$utm_y_scale, sep = "\r")
  key_id <- paste(key$utm_x_scale, key$utm_y_scale, sep = "\r")

  out_g <- key[, coord_cols, drop = FALSE]
  for (nm in needed_vars) {
    proto <- data_utm[[nm]]
    out_g[[nm]] <- proto[rep(NA_integer_, nrow(key))]
  }

  for (g in seq_len(nrow(key))) {
    rows <- which(cell_id == key_id[g])
    if (!length(rows)) {
      stop("Internal error: extrapolation cell not found in `data_utm`.", call. = FALSE)
    }

    for (nm in needed_vars) {
      vals <- data_utm[[nm]][rows]
      vals_non_na <- vals[!is.na(vals)]
      uniq <- unique(vals_non_na)
      if (length(uniq) > 1L) {
        stop(
          "Population smooth covariate `", nm, "` varies within an extrapolation grid cell. ",
          "Please supply `projection_data` explicitly. Use a `tid` column there ",
          "if the covariate varies over time.",
          call. = FALSE
        )
      }
      if (length(uniq) == 1L) {
        out_g[[nm]][g] <- uniq
      }
    }
  }

  out <- .expand_projection_over_time(out_g, tid_values = tid_values)
  rownames(out) <- NULL
  out
}

#' Prepare data objects and mesh for intCPUE workflows
#'
#' Main data-prep function for intCPUE.
#' - mesh/SPDE/A matrices
#' - extrapolation key grid + areas
#' - parse catchability and population smooths (mgcv `s()`) into Xs/Zs
#'
#' @param formula Legacy model formula. If supplied, it is treated as
#'   `formula_catchability`.
#' @param data_utm A data.frame containing required columns (with utm_x/y_scale)
#' @param mesh intCPUEmesh built from make_mesh(), or a custom mesh
#' @param formula_catchability Optional one-sided or two-sided formula defining
#'   smooth terms that affect catchability only.
#' @param formula_population Optional one-sided or two-sided formula defining
#'   smooth terms that affect both the observation model and the projected
#'   population density surface.
#' @param projection_data Optional data.frame giving projection-grid covariates
#'   for `formula_population`. It must contain `utm_x_scale` and `utm_y_scale`
#'   matching the extrapolation grid. If it also contains `tid`, it is treated
#'   as time-varying projection data with one row per grid cell-time
#'   combination. If omitted, `make_data()` attempts to derive a static
#'   projection table from `data_utm` and replicate it across time; this only
#'   works when each population-smooth covariate is uniquely defined within
#'   each grid cell. In summary:
#'   (1) static population covariates: one row per grid cell;
#'   (2) time-varying population covariates: one row per grid cell-time
#'   combination, with `tid` coded the same way as `data_utm$tid`;
#'   (3) if both static and time-varying covariates are used together, provide
#'   one row per grid cell-time combination and repeat the static covariate
#'   values across `tid`.
#' @param area_scale Numeric or "auto". Scaling factor for area_km2.
#'
#' @return A list with elements mesh, data, key, scales, smooth_basis, smooth_info.
#' @author Rujia Bi \email{rbi@@iattc.org}
#' @export
make_data <- function(
    formula = NULL,
    data_utm,
    mesh,
    formula_catchability = NULL,
    formula_population = NULL,
    projection_data = NULL,
    area_scale = "auto"
) {
  formulas <- .resolve_intCPUE_formulas(
    formula = formula,
    formula_catchability = formula_catchability,
    formula_population = formula_population,
    caller = "make_data"
  )
  formula_catchability <- formulas$formula_catchability
  formula_population <- formulas$formula_population

  data_utm <- as.data.frame(data_utm)

  .check_required_cols(data_utm, c("cpue", "encounter", "lon", "lat", "vesid", "tid", "flagid", "utm_x_scale", "utm_y_scale"))
  .check_numeric(data_utm, c("cpue", "encounter", "lon", "lat", "vesid", "tid", "flagid", "utm_x_scale", "utm_y_scale"))

  if (anyNA(data_utm$lon) || anyNA(data_utm$lat)) {
    stop("`lon`/`lat` must not contain NA.", call. = FALSE)
  }

  if (anyNA(data_utm$utm_x_scale) || anyNA(data_utm$utm_y_scale)) {
    stop("`utm_x_scale`/`utm_y_scale` must not contain NA.", call. = FALSE)
  }

  # ---- SPDE + A matrix (handle intCPUEmesh or bare mesh) ----
  loc_xy <- as.matrix(data_utm[, c("utm_x_scale", "utm_y_scale"), drop = FALSE])

  mesh_in <- mesh
  mesh_obj <- .as_intCPUEmesh(
    mesh = mesh_in,
    loc_xy = loc_xy,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    recompute_A = "auto"
  )

  if (!is.null(mesh_obj$loc_xy)) {
    r1 <- range(mesh_obj$loc_xy[, 1])
    r2 <- range(loc_xy[, 1])
    if (is.finite(r1[1]) && is.finite(r2[1])) {
      if (abs(diff(r1) - diff(r2)) / max(1e-12, diff(r2)) > 0.5) {
        warning("Mesh coordinate scale may not match `utm_x_scale/utm_y_scale`. Check scaling.", call. = FALSE)
      }
    }
  }

  mesh <- mesh_obj$mesh
  spde <- mesh_obj$spde

  # ---- A matrices ----
  A_is <- mesh_obj$A
  A_isT <- methods::as(A_is, "TsparseMatrix")
  Ais_ij <- cbind(A_isT@i, A_isT@j)
  Ais_x <- A_is@x

  # ---- key/extrapolation grid ----
  key_out <- .prep_key_area(data_utm, mesh, area_scale = area_scale)
  key <- key_out$key
  A_gs <- key_out$A_gs

  n_i <- nrow(data_utm)
  n_g <- nrow(key)

  # ---- user-supplied 0-based indices: validate only ----
  t_chk <- .check_0based_contiguous(data_utm$tid, "tid")
  v_chk <- .check_0based_contiguous(data_utm$vesid, "vesid")
  f_chk <- .check_0based_contiguous(data_utm$flagid, "flagid")

  t_i <- t_chk$x
  v_i <- v_chk$x
  f_i <- f_chk$x

  n_t <- t_chk$n
  n_v <- v_chk$n
  n_f <- f_chk$n
  tid_values <- seq.int(0L, n_t - 1L)

  # ---- smooth parsing: catchability layer ----
  sm_catch <- .normalize_smoother_output(
    parse_smoothers(
      formula = formula_catchability,
      data = data_utm,
      knots = NULL,
      newdata = NULL,
      basis_prev = NULL
    ),
    n_rows = n_i
  )

  # ---- smooth parsing: population layer (obs + projection) ----
  sm_pop_obs <- .normalize_smoother_output(
    parse_smoothers(
      formula = formula_population,
      data = data_utm,
      knots = NULL,
      newdata = NULL,
      basis_prev = NULL
    ),
    n_rows = n_i
  )

  if (isTRUE(sm_pop_obs$has_smooths)) {
    needed_vars_pop <- .smooth_vars_from_basis(sm_pop_obs$basis_out)
    projection_data_use <- if (is.null(projection_data)) {
      .default_population_projection_data(
        data_utm = data_utm,
        key = key,
        basis_out = sm_pop_obs$basis_out,
        tid_values = tid_values
      )
    } else {
      .align_projection_data_to_key(
        projection_data = projection_data,
        key = key,
        needed_vars = needed_vars_pop,
        tid_values = tid_values
      )
    }

    .warn_projection_na(
      projection_data = projection_data_use,
      needed_vars = needed_vars_pop
    )

    sm_pop_proj <- .normalize_smoother_output(
      parse_smoothers(
        formula = formula_population,
        data = data_utm,
        knots = NULL,
        newdata = projection_data_use,
        basis_prev = sm_pop_obs$basis_out
      ),
      n_rows = n_g * n_t
    )
  } else {
    projection_data_use <- .expand_projection_over_time(
      key[, c("utm_x_scale", "utm_y_scale"), drop = FALSE],
      tid_values = tid_values
    )
    sm_pop_proj <- .normalize_smoother_output(
      parse_smoothers(
        formula = NULL,
        data = data_utm,
        knots = NULL,
        newdata = projection_data_use,
        basis_prev = NULL
      ),
      n_rows = n_g * n_t
    )
  }

  # ---- build has_tf: n_t x (n_f-1), for flag-specific time effects ----
  has_tf <- matrix(FALSE, nrow = n_t, ncol = max(0L, n_f - 1L))
  if (n_f > 1L) {
    ii <- which(f_i > 0L)
    if (length(ii) > 0L) {
      tt <- t_i[ii] + 1L
      ff <- f_i[ii]
      has_tf[cbind(tt, ff)] <- TRUE
    }
  }

  data <- list(
    n_i = n_i,
    n_t = n_t,
    n_v = n_v,
    n_f = n_f,
    n_g = n_g,

    b_i = data_utm$cpue,
    e_i = data_utm$encounter,
    t_i = t_i,
    v_i = v_i,
    f_i = f_i,

    has_tf = has_tf * 1L,

    area_g = key$area_km2_scaled,

    A_is = A_is,
    A_gs = A_gs,
    Ais_ij = Ais_ij,
    Ais_x = Ais_x,

    matern_range = diff(range(mesh$loc[, 1])) / 5,
    range_prob = 0.5,
    matern_sigma_0 = 1,
    matern_sigma_t = 1,
    matern_sigma_flag = 1,
    sigma_prob = 0.05,

    has_smooths_catch = as.integer(isTRUE(sm_catch$has_smooths)),
    Xs_catch = sm_catch$Xs,
    Zs_catch = sm_catch$Zs,
    b_smooth_start_catch = as.integer(sm_catch$b_smooth_start),

    has_smooths_pop = as.integer(isTRUE(sm_pop_obs$has_smooths)),
    Xs_pop_i = sm_pop_obs$Xs,
    Zs_pop_i = sm_pop_obs$Zs,
    Xs_pop_g = sm_pop_proj$Xs,
    Zs_pop_g = sm_pop_proj$Zs,
    b_smooth_start_pop = as.integer(sm_pop_obs$b_smooth_start)
  )

  data$spde <- .prep_anisotropy(mesh = mesh, spde = spde)

  list(
    data = data,
    key = key,
    scales = list(area_scale = key_out$area_scale_val),
    projection_data = projection_data_use,
    smooth_basis = list(
      catchability = sm_catch$basis_out,
      population = sm_pop_obs$basis_out
    ),
    smooth_info = list(
      catchability = list(
        labels = sm_catch$labels,
        classes = sm_catch$classes,
        sm_dims = sm_catch$sm_dims,
        b_smooth_start = sm_catch$b_smooth_start,
        K_smooth = ncol(sm_catch$Xs),
        n_smooth = length(sm_catch$Zs)
      ),
      population = list(
        labels = sm_pop_obs$labels,
        classes = sm_pop_obs$classes,
        sm_dims = sm_pop_obs$sm_dims,
        b_smooth_start = sm_pop_obs$b_smooth_start,
        K_smooth = ncol(sm_pop_obs$Xs),
        n_smooth = length(sm_pop_obs$Zs)
      )
    )
  )
}
