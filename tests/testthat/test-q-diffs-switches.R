test_that("q-diff and vessel switches are honored in TMB when turned off", {
  skip_if_not_installed("sf")
  skip_if_not_installed("fmesher")
  skip_if_not_installed("TMB")

  data_input <- data.frame(
    cpue = c(1.2, 1e-6, 0.8, 1.1, 1e-6, 0.7, 0.9, 1e-6, 1.3, 0.6, 1e-6, 1.0),
    encounter = c(1L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L),
    lon = c(150.0, 150.1, 150.2, 150.3, 150.4, 150.5, 150.6, 150.7, 150.8, 150.9, 151.0, 151.1),
    lat = c(40.0, 40.0, 40.1, 40.1, 40.2, 40.2, 40.3, 40.3, 40.4, 40.4, 40.5, 40.5),
    vesid = rep(0:1, 6),
    tid = rep(0:2, each = 4),
    flagid = rep(c(0L, 1L), 6)
  )

  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  data_utm <- utm$data_utm

  mesh <- make_mesh(
    data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.2
  )

  fit <- suppressWarnings(
    intCPUE(
      formula = cpue ~ 1,
      data_utm = data_utm,
      mesh = mesh,
      vessel_effect = "off",
      q_diffs_system = "off",
      q_diffs_time = "off",
      q_diffs_spatial = "off",
      control = list(eval.max = 200, iter.max = 200),
      silent = TRUE
    )
  )

  rep <- fit$obj$report()

  expect_equal(as.integer(rep$use_vessel_effect), 0L)
  expect_equal(as.integer(rep$use_pop_spatial), 1L)
  expect_equal(as.integer(rep$use_pop_spatiotemporal), 1L)
  expect_equal(as.integer(rep$use_q_diffs_system), 0L)
  expect_equal(as.integer(rep$use_q_diffs_time), 0L)
  expect_equal(as.integer(rep$use_q_diffs_spatial), 0L)
  expect_s3_class(fit, "intCPUE")
})

test_that("flags with fewer than two time points are fully fixed in flag_t map", {
  parameters <- intCPUE:::.make_parameters_intCPUE(
    n_t = 4L,
    n_v = 1L,
    n_f = 3L,
    n_s = 2L,
    K_smooth = 0L,
    n_smooth = 0L,
    sum_k = 0L
  )

  has_tf <- matrix(
    c(
      TRUE,  FALSE,
      FALSE, TRUE,
      FALSE, TRUE,
      FALSE, FALSE
    ),
    nrow = 4,
    byrow = TRUE
  )

  expect_warning(
    map <- intCPUE:::.make_map_intCPUE(
      parameters = parameters,
      n_f = 3L,
      q_diffs_time = "on",
      has_tf = has_tf
    ),
    "<2 time points"
  )

  expect_true(all(is.na(map$flag_t_1[, 1])))
  expect_true(all(is.na(map$flag_t_2[, 1])))
  expect_true(is.na(map$flag_t_1[1, 2]))
  expect_false(any(is.na(map$flag_t_1[2:3, 2])))
  expect_null(map$flag_t_ln_std_dev_1)
  expect_null(map$flag_t_ln_std_dev_2)
})

test_that("population switches and observation-SD mapping are honored", {
  parameters <- intCPUE:::.make_parameters_intCPUE(
    n_t = 3L,
    n_v = 1L,
    n_f = 2L,
    n_s = 4L,
    K_smooth_catch = 0L,
    n_smooth_catch = 0L,
    sum_k_catch = 0L,
    K_smooth_pop = 0L,
    n_smooth_pop = 0L,
    sum_k_pop = 0L
  )

  map <- intCPUE:::.make_map_intCPUE(
    parameters = parameters,
    n_f = 2L,
    pop_spatial = "off",
    pop_spatiotemporal = "off",
    pop_spatiotemporal_type = "rw",
    q_diffs_time = "off",
    q_diffs_spatial = "off",
    obs_sd = "flag"
  )

  expect_true(all(is.na(map$omega_s_1)))
  expect_true(all(is.na(map$omega_s_2)))
  expect_true(all(is.na(map$epsilon_st_1)))
  expect_true(all(is.na(map$epsilon_st_2)))
  expect_true(is.na(map$ln_sigma_0_1))
  expect_true(is.na(map$ln_sigma_t_1))
  expect_true(is.na(map$ln_sigma_0_2))
  expect_true(is.na(map$ln_sigma_t_2))
  expect_true(is.na(map$ln_range_1))
  expect_true(is.na(map$ln_range_2))
  expect_true(all(is.na(map$ln_H_input)))
  expect_true(is.na(map$ln_sd))
  expect_null(map$ln_sd_flag)
})

test_that("population switches and AR1 spatiotemporal setting are reported by TMB", {
  skip_if_not_installed("sf")
  skip_if_not_installed("fmesher")
  skip_if_not_installed("TMB")

  data_input <- data.frame(
    cpue = c(1.2, 1e-6, 0.8, 1.1, 1e-6, 0.7, 0.9, 1e-6, 1.3, 0.6, 1e-6, 1.0),
    encounter = c(1L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L),
    lon = c(150.0, 150.1, 150.2, 150.3, 150.4, 150.5, 150.6, 150.7, 150.8, 150.9, 151.0, 151.1),
    lat = c(40.0, 40.0, 40.1, 40.1, 40.2, 40.2, 40.3, 40.3, 40.4, 40.4, 40.5, 40.5),
    vesid = rep(0:1, 6),
    tid = rep(0:2, each = 4),
    flagid = rep(c(0L, 1L), 6)
  )

  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  data_utm <- utm$data_utm

  mesh <- make_mesh(
    data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.2
  )

  fit_flag_sd <- suppressWarnings(
    intCPUE(
      formula = cpue ~ 1,
      data_utm = data_utm,
      mesh = mesh,
      pop_spatial = "off",
      pop_spatiotemporal = "off",
      obs_sd = "flag",
      vessel_effect = "off",
      q_diffs_system = "off",
      q_diffs_time = "off",
      q_diffs_spatial = "off",
      control = list(eval.max = 200, iter.max = 200),
      silent = TRUE
    )
  )

  rep_flag_sd <- fit_flag_sd$obj$report()
  expect_equal(as.integer(rep_flag_sd$use_pop_spatial), 0L)
  expect_equal(as.integer(rep_flag_sd$use_pop_spatiotemporal), 0L)
  expect_equal(as.integer(rep_flag_sd$use_pop_spatiotemporal_rw), 0L)
  expect_equal(as.integer(rep_flag_sd$use_pop_spatiotemporal_ar1), 0L)
  expect_equal(as.integer(rep_flag_sd$use_flag_sd), 1L)
  expect_length(rep_flag_sd$sd_flag, fit_flag_sd$data_tmb$n_f)

  fit_ar1 <- suppressWarnings(
    intCPUE(
      formula = cpue ~ 1,
      data_utm = data_utm,
      mesh = mesh,
      pop_spatial = "off",
      pop_spatiotemporal = "on",
      pop_spatiotemporal_type = "ar1",
      obs_sd = "shared",
      vessel_effect = "off",
      q_diffs_system = "off",
      q_diffs_time = "off",
      q_diffs_spatial = "off",
      control = list(eval.max = 200, iter.max = 200),
      silent = TRUE
    )
  )

  rep_ar1 <- fit_ar1$obj$report()
  expect_equal(as.integer(rep_ar1$use_pop_spatial), 0L)
  expect_equal(as.integer(rep_ar1$use_pop_spatiotemporal), 1L)
  expect_equal(as.integer(rep_ar1$use_pop_spatiotemporal_rw), 0L)
  expect_equal(as.integer(rep_ar1$use_pop_spatiotemporal_ar1), 1L)
  expect_equal(as.integer(rep_ar1$use_flag_sd), 0L)
})
