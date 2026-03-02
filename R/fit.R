#' @useDynLib intCPUE, .registration = TRUE
NULL

#' Fit an integrated spatiotemporal CPUE standardization model
#'
#' Fits an integrated spatiotemporal model for CPUE standardization that jointly
#' models observation and sampling processes across one or more fisheries or surveys.
#' The framework is designed to reduce bias caused by preferential sampling,
#' targeting behavior, and heterogeneous effort distributions, and to estimate
#' a coherent latent abundance index.
#'
#' The model is implemented using Template Model Builder (TMB). Spatial and
#' spatiotemporal random fields are represented using the SPDE (stochastic partial
#' differential equation) approach for computational efficiency. See the model
#' description vignette for details:
#' https://github.com/RujiaBi/intCPUE/blob/main/vignettes/intCPUE-intro.Rmd
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
#' @param pop_spatial "on" or "off". If "off", the population spatial field is
#'   fixed at 0.
#' @param pop_spatiotemporal "on" or "off". If "off", the population
#'   spatiotemporal field is fixed at 0.
#' @param pop_spatiotemporal_type One of `"rw"`, `"iid"`, or `"ar1"`.
#'   Controls the temporal dependence structure of the population
#'   spatiotemporal field when `pop_spatiotemporal = "on"`. Defaults to `"rw"`
#'   to preserve the previous model behavior.
#' @param vessel_effect "on" or "off". If "off", vessel RE is fixed at 0 via map.
#' @param q_diffs_system "on" or "off". If "off", flag systematic differences are fixed.
#' @param q_diffs_time "on" or "off". If "off", flag-specific temporal effects are fixed.
#' @param q_diffs_spatial "on" or "off". If "off", flag-specific spatial fields are fixed.
#' @param obs_sd `"shared"` or `"flag"`. If `"flag"`, the positive-catch
#'   lognormal observation SD is estimated separately for each flag.
#' @param control Control list passed to [stats::nlminb()].
#' @param ncores Optional integer. If provided, sets the number of OpenMP threads. Passed to [TMB::openmp()].
#' @param silent Logical. Passed to [TMB::MakeADFun()].
#' @param ... Passed to [intCPUE::make_data()] (for example, `area_scale`).
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
    pop_spatial = c("on", "off"),
    pop_spatiotemporal = c("on", "off"),
    pop_spatiotemporal_type = c("rw", "iid", "ar1"),
    vessel_effect = c("on", "off"),
    q_diffs_system = c("on", "off"),
    q_diffs_time = c("on", "off"),
    q_diffs_spatial = c("on", "off"),
    obs_sd = c("shared", "flag"),
    control = list(eval.max = 1e5, iter.max = 1e5),
    ncores = NULL,
    ...,
    silent = FALSE
) {
  pop_spatial <- match.arg(pop_spatial)
  pop_spatiotemporal <- match.arg(pop_spatiotemporal)
  pop_spatiotemporal_type <- match.arg(pop_spatiotemporal_type)
  vessel_effect   <- match.arg(vessel_effect)
  q_diffs_system  <- match.arg(q_diffs_system)
  q_diffs_time    <- match.arg(q_diffs_time)
  q_diffs_spatial <- match.arg(q_diffs_spatial)
  obs_sd <- match.arg(obs_sd)
  
  data_utm <- as.data.frame(data_utm)
  
  # ---- 1) Data prep (single source of truth) ----
  prep <- make_data(
    formula = formula,
    data_utm = data_utm,
    mesh = mesh,
    formula_catchability = formula_catchability,
    formula_population = formula_population,
    projection_data = projection_data,
    ...
  )
  
  data_tmb <- prep$data
  
  # ---- 2) Defensive checks (catch mismatches early) ----
  n_t <- data_tmb$n_t
  n_v <- data_tmb$n_v
  n_f <- data_tmb$n_f
  n_s <- data_tmb$spde$n_s

  data_tmb$use_vessel_effect <- as.integer(vessel_effect == "on" && n_v > 0L)
  data_tmb$use_pop_spatial <- as.integer(pop_spatial == "on")
  data_tmb$use_pop_spatiotemporal <- as.integer(pop_spatiotemporal == "on" && n_t > 0L)
  data_tmb$use_pop_spatiotemporal_rw <- as.integer(
    pop_spatiotemporal == "on" && pop_spatiotemporal_type == "rw" && n_t > 0L
  )
  data_tmb$use_pop_spatiotemporal_ar1 <- as.integer(
    pop_spatiotemporal == "on" && pop_spatiotemporal_type == "ar1" && n_t > 0L
  )
  data_tmb$use_q_diffs_system <- as.integer(q_diffs_system == "on" && n_f > 1L)
  data_tmb$use_q_diffs_time <- as.integer(q_diffs_time == "on" && n_f > 1L)
  data_tmb$use_q_diffs_spatial <- as.integer(q_diffs_spatial == "on" && n_f > 1L)
  data_tmb$use_flag_sd <- as.integer(obs_sd == "flag" && n_f > 0L)
  
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
    n_t = n_t, n_v = n_v, n_f = n_f, n_s = n_s,
    K_smooth_catch = K_smooth_catch,
    n_smooth_catch = n_smooth_catch,
    sum_k_catch = sum_k_catch,
    K_smooth_pop = K_smooth_pop,
    n_smooth_pop = n_smooth_pop,
    sum_k_pop = sum_k_pop
  )
  
  # ---- 4) MAP (turn on/off components without touching cpp) ----
  has_tf <- NULL
  if (!is.null(data_tmb$has_tf)) {
    # has_tf is integer matrix 0/1 in your make_data()
    has_tf <- (data_tmb$has_tf > 0L)
  }
  
  map <- .make_map_intCPUE(
    parameters = parameters,
    n_f = n_f,
    pop_spatial = pop_spatial,
    pop_spatiotemporal = pop_spatiotemporal,
    pop_spatiotemporal_type = pop_spatiotemporal_type,
    vessel_effect = vessel_effect,
    q_diffs_system = q_diffs_system,
    q_diffs_time = q_diffs_time,
    q_diffs_spatial = q_diffs_spatial,
    obs_sd = obs_sd,
    has_tf = has_tf
  )
  
  # ---- 5) Random effects list ----
  random <- character(0)

  if (pop_spatial == "on") {
    random <- c(random, "omega_s_1", "omega_s_2")
  }

  if (pop_spatiotemporal == "on") {
    random <- c(random, "epsilon_st_1", "epsilon_st_2")
  }
  
  if (vessel_effect == "on") {
    random <- c(random, "ves_v_1", "ves_v_2")
  }
  
  # Only include these if flags exist:
  if (n_f > 1L && q_diffs_system == "on") {
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
  
  random <- unique(random)
  
  # ---- 6) Build & optimize ----
  .check_fit_inputs_intCPUE(data_tmb)
  
  DLL <- "intCPUE"
  
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
  
  opt <- .safe_optimize(obj, control)
  
  rep <- TMB::sdreport(obj)
  
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
      pop_spatial = pop_spatial,
      pop_spatiotemporal = pop_spatiotemporal,
      pop_spatiotemporal_type = pop_spatiotemporal_type,
      vessel_effect = vessel_effect,
      q_diffs_system = q_diffs_system,
      q_diffs_time = q_diffs_time,
      q_diffs_spatial = q_diffs_spatial,
      obs_sd = obs_sd,
      DLL = DLL,
      ncores = ncores
    ),
    diagnostics = list(
      convergence = opt$convergence,
      message = opt$message,
      max_grad = max(abs(obj$gr(opt$par)))
    )
  )
  class(out) <- "intCPUE"
  out
}
