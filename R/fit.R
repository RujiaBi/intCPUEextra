#' @useDynLib intCPUEextra, .registration = TRUE
NULL

#' Fit a multi-area CPUE standardization model based on intCPUE
#'
#' Fits a multi-area extension of `intCPUE` for CPUE standardization. This
#' version is designed for settings where vessel effects and catchability
#' smooths are shared across areas, fishery-level mean differences and
#' flag-specific temporal deviations are shared across areas, and temporal
#' effects and latent population structure are estimated separately by area.
#' Flag-specific spatial deviations are modeled on a dedicated full-domain
#' mesh via `flag_mesh`, with their own anisotropy, range, and SD parameters.
#'
#' The model is implemented using Template Model Builder (TMB). Spatial and
#' spatiotemporal random fields are represented using the SPDE approach for
#' computational efficiency. See the vignette for details:
#' https://github.com/RujiaBi/intCPUEextra/blob/main/vignettes/intCPUEextra-intro.Rmd
#'
#' @param formula Legacy one-sided or two-sided model formula with optional
#'   mgcv::s() smooth terms. If supplied, it is treated as
#'   `formula_catchability`.
#' @param data_utm A data.frame with required columns.
#' @param mesh An `intCPUEmesh` or a bare fmesher mesh.
#' @param formula_catchability Optional one-sided or two-sided formula defining
#'   smooth terms that affect catchability only.
#' @param formula_population Optional one-sided or two-sided formula defining
#'   smooth terms that affect the latent population surface and therefore also
#'   enter projection.
#' @param projection_data Optional projection-grid data for population smooth
#'   covariates. Passed to [make_data()]. For static population covariates,
#'   supply one row per extrapolation grid cell with columns `utm_x_scale`,
#'   `utm_y_scale`, and the covariate names used in `formula_population`. For
#'   time-varying population covariates, additionally include `tid` and supply
#'   one row per grid cell-time combination. If `formula_population` mixes
#'   static and time-varying covariates, use the grid cell-time format and
#'   repeat static covariate values across `tid`.
#' @param area_col Optional area column in `data_utm`. If supplied and it
#'   identifies multiple areas, areas share vessel effects and catchability
#'   smooths and flag-specific non-spatial effects, while fishery-level mean
#'   differences and temporal/spatial population structure remain area-specific.
#'   Use a named list of meshes with one entry per area when fitting multiple
#'   areas.
#' @param flag_mesh Optional dedicated mesh for the shared flag-specific
#'   spatial field (`flag_s`). When fitting multiple areas with
#'   `q_diffs_spatial = "on"`, supply a single full-domain mesh here.
#' @param q_diffs_time "on" or "off". Controls whether flag-specific temporal
#'   deviations are included. Implemented via dedicated TMB templates rather
#'   than `map`.
#' @param q_diffs_spatial "on" or "off". Controls whether flag-specific spatial
#'   deviations are included. Implemented via dedicated TMB templates rather
#'   than `map`.
#' @param pop_spatiotemporal_type `"rw"` or `"ar1"`. Controls whether the
#'   always-on spatiotemporal population field follows a random-walk or AR1
#'   evolution over time.
#' @param obs_sd_flag `"shared"` or `"flag"`. Controls whether the
#'   positive-catch lognormal observation SD is shared across flags or
#'   estimated separately for each flag.
#' @param obs_sd_area `"area"` or `"shared"`. Controls whether the
#'   positive-catch lognormal observation SD is estimated separately by area
#'   or shared across areas. The default `"area"` preserves the historical
#'   package behavior.
#' @param t_sd Optional positive scalar. Only used when `formula_population`
#'   contains `f(tid)`. The default `NULL` estimates the random time-effect
#'   standard deviation; a numeric value (for example `10`) fixes that standard
#'   deviation at the supplied value.
#' @param control Control list passed to [stats::nlminb()].
#' @param ncores Optional integer. If provided, sets the number of OpenMP threads. Passed to [TMB::openmp()].
#' @param ... Passed to [intCPUEextra::make_data()] (for example, `area_scale`).
#' @param silent Logical. Passed to [TMB::MakeADFun()].
#' @param restart_max Non-negative integer. Maximum number of automatic
#'   `nlminb()` restarts attempted from the current best parameter vector when
#'   convergence or gradient checks remain unsatisfactory.
#' @param newton_max Non-negative integer. Maximum number of
#'   [TMB::newton()] refinement attempts after `nlminb()`/restart stages.
#' @param coord_max Non-negative integer. Maximum number of coordinate-polish
#'   iterations attempted after the restart/Newton stages. Set
#'   `restart_max = 0`, `newton_max = 0`, and `coord_max = 0` to run only a
#'   single `nlminb()` pass.
#'
#' @return An object of class `intCPUE` with elements `obj`, `opt`, `rep`, `prep`, etc.
#' @author Rujia Bi \email{rbi@@iattc.org}
#' @export
intCPUE <- function(
    formula = NULL,
    data_utm,
    mesh,
    formula_catchability = NULL,
    formula_population = NULL,
    projection_data = NULL,
    area_col = NULL,
    flag_mesh = NULL,
    q_diffs_time = c("on", "off"),
    q_diffs_spatial = c("on", "off"),
    pop_spatiotemporal_type = c("rw", "ar1"),
    obs_sd_flag = c("shared", "flag"),
    obs_sd_area = c("area", "shared"),
    t_sd = NULL,
    control = list(eval.max = 1e5, iter.max = 1e5),
    ncores = NULL,
    ...,
    silent = FALSE,
    restart_max = 1L,
    newton_max = 2L,
    coord_max = 5L
) {
  q_diffs_time <- match.arg(q_diffs_time)
  q_diffs_spatial <- match.arg(q_diffs_spatial)
  pop_spatiotemporal_type <- match.arg(pop_spatiotemporal_type)
  obs_sd_flag <- match.arg(obs_sd_flag)
  obs_sd_area <- match.arg(obs_sd_area)
  if (is.character(t_sd) && length(t_sd) == 1L && identical(tolower(t_sd), "null")) {
    t_sd <- NULL
  }
  if (!is.null(t_sd)) {
    if (!is.numeric(t_sd) || length(t_sd) != 1L || !is.finite(t_sd) || t_sd <= 0) {
      stop("`t_sd` must be NULL or a positive finite numeric scalar.", call. = FALSE)
    }
    t_sd <- as.numeric(t_sd)
  }
  restart_max <- .validate_nonneg_count(restart_max, "restart_max")
  newton_max <- .validate_nonneg_count(newton_max, "newton_max")
  coord_max <- .validate_nonneg_count(coord_max, "coord_max")
  
  data_utm <- as.data.frame(data_utm)

  if (!is.null(area_col) && q_diffs_spatial == "on" && area_col %in% names(data_utm)) {
    n_areas_input <- length(unique(as.character(data_utm[[area_col]])))
    if (n_areas_input > 1L && is.null(flag_mesh)) {
      stop(
        "Multi-area fits with `q_diffs_spatial = 'on'` require `flag_mesh` ",
        "so the shared flag spatial field can use a dedicated full-domain mesh.",
        call. = FALSE
      )
    }
  }
  
  # ---- 1) Data prep (single source of truth) ----
  prep <- make_data(
    formula = formula,
    data_utm = data_utm,
    mesh = mesh,
    formula_catchability = formula_catchability,
    formula_population = formula_population,
    projection_data = projection_data,
    area_col = area_col,
    flag_mesh = flag_mesh,
    q_diffs_spatial = q_diffs_spatial,
    ...
  )
  
  data_tmb <- prep$data
  use_random_tid_effect <- isTRUE(as.integer(data_tmb$use_random_tid_effect) == 1L)
  if (!use_random_tid_effect && !is.null(t_sd)) {
    stop("`t_sd` can only be used when `formula_population` contains `f(tid)`.", call. = FALSE)
  }
  data_tmb$fix_t_sd <- as.integer(use_random_tid_effect && !is.null(t_sd))
  data_tmb$t_sd <- if (is.null(t_sd)) 1.0 else t_sd
  
  # ---- 2) Defensive checks (catch mismatches early) ----
  n_a <- data_tmb$n_a
  n_t <- data_tmb$n_t
  n_v <- data_tmb$n_v
  n_f <- data_tmb$n_f
  n_s <- sum(data_tmb$n_s_area)
  n_s_flag <- data_tmb$n_s_flag

  data_tmb$use_q_diffs_time <- as.integer(q_diffs_time == "on" && n_f > 1L)
  data_tmb$use_q_diffs_spatial <- as.integer(q_diffs_spatial == "on" && n_f > 1L)
  data_tmb$use_pop_spatiotemporal_rw <- as.integer(pop_spatiotemporal_type == "rw")
  data_tmb$use_pop_spatiotemporal_ar1 <- as.integer(pop_spatiotemporal_type == "ar1")
  data_tmb$use_flag_sd <- as.integer(obs_sd_flag == "flag" && n_f > 0L)

  has_tf <- NULL
  flag_t_constraint <- NULL
  if (q_diffs_time == "on" && n_f > 1L) {
    has_tf <- data_tmb$has_tf > 0L
    flag_t_constraint <- .build_flag_t_constraint(has_tf, n_f = n_f)
    data_tmb$flag_t_index <- matrix(
      as.integer(flag_t_constraint$flag_t_index),
      nrow = nrow(flag_t_constraint$flag_t_index),
      ncol = ncol(flag_t_constraint$flag_t_index)
    )
  } else {
    data_tmb$flag_t_index <- matrix(
      integer(n_t * max(0L, n_f - 1L)),
      nrow = n_t,
      ncol = max(0L, n_f - 1L)
    )
  }
  
  # Smooth dims
  has_smooths_catch <- isTRUE(data_tmb$has_smooths_catch == 1L)
  K_smooth_catch <- if (has_smooths_catch) ncol(data_tmb$Xs_catch) else 0L
  n_smooth_catch <- if (has_smooths_catch) length(data_tmb$Zs_catch) else 0L
  sum_k_catch <- if (has_smooths_catch && n_smooth_catch > 0L) sum(vapply(data_tmb$Zs_catch, ncol, 0L)) else 0L

  has_smooths_pop <- isTRUE(data_tmb$has_smooths_pop == 1L)
  K_smooth_pop <- if (has_smooths_pop) ncol(data_tmb$Xs_pop_i) else 0L
  n_smooth_pop <- if (has_smooths_pop) length(data_tmb$Zs_pop_i) else 0L
  sum_k_pop <- if (has_smooths_pop && n_smooth_pop > 0L) sum(vapply(data_tmb$Zs_pop_i, ncol, 0L)) else 0L
  
  # ---- 3) Initial parameters (must match cpp) ----
  parameters <- .make_parameters_intCPUE(
    n_a = n_a,
    n_t = n_t, n_v = n_v, n_f = n_f, n_s = n_s,
    n_s_flag = n_s_flag,
    n_flag_t_free = if (is.null(flag_t_constraint)) 0L else flag_t_constraint$n_free,
    use_q_diffs_time = (q_diffs_time == "on" && n_f > 1L),
    use_q_diffs_spatial = (q_diffs_spatial == "on" && n_f > 1L),
    K_smooth_catch = K_smooth_catch,
    n_smooth_catch = n_smooth_catch,
    sum_k_catch = sum_k_catch,
    K_smooth_pop = K_smooth_pop,
    n_smooth_pop = n_smooth_pop,
    sum_k_pop = sum_k_pop
  )

  # Match the active TMB template exactly: q-time/q-spatial parameters do not
  # exist at all in the corresponding off-templates.
  if (q_diffs_time == "off" || n_f <= 1L) {
    parameters$flag_t_1 <- NULL
    parameters$flag_t_2 <- NULL
    parameters$flag_t_ln_std_dev_1 <- NULL
    parameters$flag_t_ln_std_dev_2 <- NULL
  }
  if (q_diffs_spatial == "off" || n_f <= 1L) {
    parameters$ln_H_flag_input <- NULL
    parameters$flag_s_1 <- NULL
    parameters$flag_s_2 <- NULL
    parameters$ln_range_flag_1 <- NULL
    parameters$ln_range_flag_2 <- NULL
    parameters$ln_sigma_flag_1 <- NULL
    parameters$ln_sigma_flag_2 <- NULL
  }
  .check_parameter_shapes_intCPUE(
    parameters = parameters,
    data_tmb = data_tmb,
    q_diffs_time = q_diffs_time,
    q_diffs_spatial = q_diffs_spatial
  )
  
  # ---- 4) MAP (turn on/off components without touching cpp) ----
  map <- .make_map_intCPUE(
    parameters = parameters,
    n_f = n_f,
    obs_sd_flag = obs_sd_flag,
    obs_sd_area = obs_sd_area,
    use_random_tid_effect = use_random_tid_effect,
    fix_t_sd = !is.null(t_sd),
    pop_spatiotemporal_type = pop_spatiotemporal_type,
    q_diffs_time = q_diffs_time,
    has_tf = has_tf,
    estimable_flag = if (is.null(flag_t_constraint)) NULL else flag_t_constraint$estimable_flag
  )
  
  # ---- 5) Random effects list ----
  random <- character(0)

  random <- c(random, "omega_s_1", "omega_s_2")
  random <- c(random, "epsilon_st_1", "epsilon_st_2")

  if (n_v > 0L) {
    random <- c(random, "ves_v_1", "ves_v_2")
  }
  
  if (n_f > 1L) {
    random <- c(random, "flag_f_1", "flag_f_2")
  }
  if (n_f > 1L && q_diffs_time == "on") {
    random <- c(random, "flag_t_1", "flag_t_2")
  }
  if (n_f > 1L && q_diffs_spatial == "on") {
    random <- c(random, "flag_s_1", "flag_s_2")
  }
  
  # Smooth random coeffs
  if (has_smooths_catch && sum_k_catch > 0L) random <- c(random, "b_smooth_catch")
  if (has_smooths_pop && sum_k_pop > 0L) random <- c(random, "b_smooth_pop")
  if (use_random_tid_effect) random <- c(random, "yq_t_1", "yq_t_2")
  
  random <- unique(random)
  
  # ---- 6) Build & optimize ----
  .check_fit_inputs_intCPUE(data_tmb)
  
  DLL <- .ensure_template_dll_intCPUE(
    q_diffs_time = q_diffs_time,
    q_diffs_spatial = q_diffs_spatial
  )
  
  if (!is.null(ncores)) {
    ncores <- as.integer(ncores)
    if (is.na(ncores) || ncores < 1L) {
      stop("`ncores` must be a positive integer.", call. = FALSE)
    }
    if (ncores > 1L) {
      TMB::openmp(ncores, autopar = TRUE, DLL = DLL)
    }
  }
  
  obj <- TMB::MakeADFun(
    data = data_tmb,
    parameters = parameters,
    map = map,
    random = random,
    DLL = DLL,
    silent = silent
  )
  
  opt <- .safe_optimize(
    obj = obj,
    control = control,
    restart_max = restart_max,
    newton_max = newton_max,
    coord_max = coord_max
  )
  par_structured <- try(obj$env$parList(opt$par), silent = TRUE)
  if (!inherits(par_structured, "try-error")) {
    opt$par_list <- par_structured
  }
  
  rep_info <- .safe_sdreport(obj, opt)
  rep <- rep_info$rep
  
  out <- list(
    obj = obj,
    opt = opt,
    rep = rep,
    prep = prep,
    data_tmb = data_tmb,
    map = map,
    random = random,
    control = control,
    settings = list(
      formula = formula,
      formula_catchability = if (is.null(formula_catchability)) formula else formula_catchability,
      formula_population = formula_population,
      area_col = area_col,
      flag_mesh = flag_mesh,
      pop_spatiotemporal_type = pop_spatiotemporal_type,
      q_diffs_time = q_diffs_time,
      q_diffs_spatial = q_diffs_spatial,
      obs_sd_flag = obs_sd_flag,
      obs_sd_area = obs_sd_area,
      t_sd = t_sd,
      use_random_tid_effect = use_random_tid_effect,
      DLL = DLL,
      ncores = ncores,
      restart_max = restart_max,
      newton_max = newton_max,
      coord_max = coord_max
    ),
    diagnostics = list(
      convergence = opt$convergence,
      message = opt$message,
      max_grad = if (!is.null(opt$max_grad)) opt$max_grad else max(abs(obj$gr(opt$par))),
      newton_used = if (!is.null(opt$newton_used)) opt$newton_used else NA_integer_,
      coord_used = if (!is.null(opt$coord_used)) opt$coord_used else NA_integer_,
      sdreport_method = rep_info$method,
      sdreport_error = rep_info$error
    )
  )
  class(out) <- "intCPUE"
  out
}
