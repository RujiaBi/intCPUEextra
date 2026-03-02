# Internal: build TMB map to switch model components on/off
#
# The map controls which parameters are estimated vs fixed.
# In TMB, map entries set to NA are fixed at their initial value.
#
# Important statistical principle:
# We automatically prevent estimation of variance parameters that are not identifiable
# (e.g., flag temporal random effects with < 2 time points).
#
.make_map_intCPUE <- function(parameters, n_f,
                              pop_spatial = c("on","off"),
                              pop_spatiotemporal = c("on","off"),
                              pop_spatiotemporal_type = c("rw","iid","ar1"),
                              vessel_effect = c("on","off"),
                              q_diffs_system = c("on","off"),
                              q_diffs_time = c("on","off"),
                              q_diffs_spatial = c("on","off"),
                              obs_sd = c("shared","flag"),
                              has_tf = NULL) {

  pop_spatial <- match.arg(pop_spatial)
  pop_spatiotemporal <- match.arg(pop_spatiotemporal)
  pop_spatiotemporal_type <- match.arg(pop_spatiotemporal_type)
  vessel_effect  <- match.arg(vessel_effect)
  q_diffs_system <- match.arg(q_diffs_system)
  q_diffs_time   <- match.arg(q_diffs_time)
  q_diffs_spatial<- match.arg(q_diffs_spatial)
  obs_sd <- match.arg(obs_sd)
  
  map <- list()

  # ---------------------------------------------------------
  # Observation SD
  # ---------------------------------------------------------
  if (obs_sd == "shared") {
    map$ln_sd_flag <- factor(rep(NA, length(parameters$ln_sd_flag)))
  } else {
    map$ln_sd <- factor(NA)
  }

  # ---------------------------------------------------------
  # Population spatial fields
  # ---------------------------------------------------------
  if (pop_spatial == "off") {
    map$omega_s_1 <- factor(rep(NA, length(parameters$omega_s_1)))
    map$omega_s_2 <- factor(rep(NA, length(parameters$omega_s_2)))
    map$ln_sigma_0_1 <- factor(NA)
    map$ln_sigma_0_2 <- factor(NA)
  }

  # ---------------------------------------------------------
  # Population spatiotemporal fields
  # ---------------------------------------------------------
  if (pop_spatiotemporal == "off") {
    map$epsilon_st_1 <- .map_matrix_NA(parameters$epsilon_st_1)
    map$epsilon_st_2 <- .map_matrix_NA(parameters$epsilon_st_2)
    map$ln_sigma_t_1 <- factor(NA)
    map$ln_sigma_t_2 <- factor(NA)
    map$transf_rho_1 <- factor(NA)
    map$transf_rho_2 <- factor(NA)
  } else if (pop_spatiotemporal_type != "ar1") {
    map$transf_rho_1 <- factor(NA)
    map$transf_rho_2 <- factor(NA)
  }
  
  # ---------------------------------------------------------
  # Vessel random effects
  # ---------------------------------------------------------
  if (vessel_effect == "off") {
    map$ves_v_1 <- factor(rep(NA, length(parameters$ves_v_1)))
    map$ves_v_2 <- factor(rep(NA, length(parameters$ves_v_2)))
    map$ves_ln_std_dev_1 <- factor(NA)
    map$ves_ln_std_dev_2 <- factor(NA)
  }
  
  # ---------------------------------------------------------
  # Systematic flag differences (constant offsets among flags)
  # ---------------------------------------------------------
  if (q_diffs_system == "off" || n_f <= 1L) {
    map$flag_f_1 <- factor(rep(NA, length(parameters$flag_f_1)))
    map$flag_f_2 <- factor(rep(NA, length(parameters$flag_f_2)))
    map$flag_ln_std_dev_1 <- factor(NA)
    map$flag_ln_std_dev_2 <- factor(NA)
  }
  
  # ==========================================================
  # Temporal flag differences (flag_t)
  # ==========================================================
  # flag_t is a per-(time, flag) deviation relative to flag_f.
  #
  # IMPORTANT:
  # - Missing (t,flag) cells should be fixed to 0 (no extrapolation).
  # - To avoid confounding with flag_f, we recommend centering flag_t
  #   in the C++ (sum-to-zero over OBSERVED times per flag).
  #   The map alone cannot enforce sum-to-zero constraints.
  #
  # NOTE on identifiability:
  # - `flag_t_ln_std_dev_*` is a SINGLE global SD (shared across flags)
  #   in your current CPP, so you cannot "fix variance for one flag only".
  # - We only fix the global SD if NO flags have >=2 observed time points.
  #
  if (q_diffs_time == "off" || n_f <= 1L) {
    
    map$flag_t_1 <- .map_matrix_NA(parameters$flag_t_1)
    map$flag_t_2 <- .map_matrix_NA(parameters$flag_t_2)
    map$flag_t_ln_std_dev_1 <- factor(NA)
    map$flag_t_ln_std_dev_2 <- factor(NA)
    
  } else {
    
    if (is.null(has_tf)) {
      stop("q_diffs_time='on' requires `has_tf` (n_t x (n_f-1)) to fix missing (t,flag) cells to 0.", call. = FALSE)
    }
    
    if (!is.matrix(has_tf)) stop("`has_tf` must be a matrix.", call. = FALSE)
    if (!all(dim(has_tf) == dim(parameters$flag_t_1))) {
      stop("`has_tf` dim must match `flag_t_1` dim: ",
           paste(dim(has_tf), collapse = "x"), " vs ",
           paste(dim(parameters$flag_t_1), collapse = "x"),
           call. = FALSE)
    }
    
    # How many observed time points per flag column?
    n_time_per_flag <- colSums(has_tf > 0)
    estimable_flag  <- n_time_per_flag >= 2

    keep_tf <- has_tf
    if (any(!estimable_flag)) {
      keep_tf[, !estimable_flag] <- FALSE
    }

    # Fix missing (t,flag) cells to 0 via map (NA => fixed at initial value 0).
    # If a flag has <2 observed time points, its whole temporal column is fixed to 0.
    map$flag_t_1 <- .map_matrix_partial_fix(parameters$flag_t_1, keep_tf)
    map$flag_t_2 <- .map_matrix_partial_fix(parameters$flag_t_2, keep_tf)
    
    if (any(!estimable_flag)) {
      warning(
        sprintf(
          paste0(
            "Some flags have <2 time points (flags: %s). ",
            "Their entire temporal deviation columns are fixed to 0. ",
            "The global temporal SD is fixed only if NO flags have >=2 time points."
          ),
          paste(which(!estimable_flag), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    
    # If no flag has >=2 observed time points, the global SD is not identifiable
    if (!any(estimable_flag)) {
      map$flag_t_ln_std_dev_1 <- factor(NA)
      map$flag_t_ln_std_dev_2 <- factor(NA)
    }
  }
  
  # ---------------------------------------------------------
  # Spatial flag differences
  # ---------------------------------------------------------
  if (q_diffs_spatial == "off" || n_f <= 1L) {
    map$flag_s_1 <- .map_matrix_NA(parameters$flag_s_1)
    map$flag_s_2 <- .map_matrix_NA(parameters$flag_s_2)
    map$ln_sigma_flag_1 <- factor(NA)
    map$ln_sigma_flag_2 <- factor(NA)
  }

  use_any_spde <- (
    pop_spatial == "on" ||
    pop_spatiotemporal == "on" ||
    (q_diffs_spatial == "on" && n_f > 1L)
  )

  if (!use_any_spde) {
    map$ln_range_1 <- factor(NA)
    map$ln_range_2 <- factor(NA)
    map$ln_H_input <- factor(rep(NA, length(parameters$ln_H_input)))
  }
  
  map
}
