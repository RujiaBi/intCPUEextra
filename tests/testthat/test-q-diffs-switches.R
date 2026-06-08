test_that("q-diff switches are honored in TMB when turned off", {
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
      q_diffs_time = "off",
      q_diffs_spatial = "off",
      control = list(eval.max = 200, iter.max = 200),
      silent = TRUE
    )
  )

  rep <- fit$obj$report()

  expect_equal(as.integer(rep$use_q_diffs_time), 0L)
  expect_equal(as.integer(rep$use_q_diffs_spatial), 0L)
  expect_equal(as.integer(rep$use_pop_spatiotemporal_rw), 1L)
  expect_equal(as.integer(rep$use_pop_spatiotemporal_ar1), 0L)
  expect_s3_class(fit, "intCPUE")
})

test_that("flags with fewer than two time points are excluded from flag_t constraints", {
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

  constraint <- intCPUEextra:::.build_flag_t_constraint(has_tf)

  expect_equal(constraint$estimable_flag, c(FALSE, TRUE))
  expect_true(constraint$any_estimable)
  expect_equal(constraint$n_free, 1L)
  expect_equal(
    constraint$flag_t_index,
    matrix(
      c(
        -1L, -1L,
        -1L, -1L,
        -1L,  0L,
        -1L, -1L
      ),
      nrow = 4,
      byrow = TRUE
    )
  )
})

test_that("multi-area q_diffs_spatial requires a dedicated flag mesh", {
  skip_if_not_installed("sf")
  skip_if_not_installed("fmesher")

  data_input <- data.frame(
    cpue = c(1.2, 0.0, 0.8, 1.1, 0.9, 0.0, 1.3, 0.7),
    encounter = c(1L, 0L, 1L, 1L, 1L, 0L, 1L, 1L),
    lon = c(150.0, 150.1, 150.2, 150.3, 151.0, 151.1, 151.2, 151.3),
    lat = c(40.0, 40.0, 40.1, 40.1, 39.8, 39.8, 39.9, 39.9),
    vesid = rep(0:1, 4),
    tid = rep(0:1, each = 4),
    flagid = rep(c(0L, 1L), 4),
    area = rep(c("west", "east"), each = 4)
  )

  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  data_utm <- utm$data_utm
  mesh_list <- lapply(split(data_utm, data_utm$area), function(dd) {
    make_mesh(
      dd,
      xy_cols = c("utm_x_scale", "utm_y_scale"),
      type = "cutoff",
      cutoff = 0.2
    )
  })

  expect_error(
    intCPUE(
      data_utm = data_utm,
      mesh = mesh_list,
      area_col = "area",
      q_diffs_spatial = "on",
      q_diffs_time = "off"
    ),
    "require `flag_mesh`"
  )
})

test_that("make_data skips dedicated flag mesh objects when q_diffs_spatial is off", {
  skip_if_not_installed("sf")
  skip_if_not_installed("fmesher")

  data_input <- data.frame(
    cpue = c(1.2, 0, 0.8, 1.1, 0, 0.7, 0.9, 0),
    encounter = c(1L, 0L, 1L, 1L, 0L, 1L, 1L, 0L),
    lon = c(150.0, 150.1, 150.2, 150.3, 160.0, 160.1, 160.2, 160.3),
    lat = c(40.0, 40.0, 40.1, 40.1, 35.0, 35.0, 35.1, 35.1),
    vesid = rep(0:1, 4),
    tid = rep(0:1, each = 4),
    flagid = rep(c(0L, 1L), 4),
    region = rep(c("A", "B"), 4)
  )

  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  data_utm <- utm$data_utm
  area_levels <- unique(as.character(data_utm$region))
  mesh <- setNames(
    lapply(area_levels, function(a) {
      make_mesh(
        data_utm[data_utm$region == a, ],
        xy_cols = c("utm_x_scale", "utm_y_scale"),
        type = "cutoff",
        cutoff = 0.2
      )
    }),
    area_levels
  )

  prep <- make_data(
    data_utm = data_utm,
    mesh = mesh,
    area_col = "region",
    q_diffs_spatial = "off"
  )

  expect_equal(prep$data$n_s_flag, 0L)
  expect_true(inherits(prep$data$A_flag_is, "sparseMatrix"))
  expect_equal(dim(prep$data$A_flag_is), c(nrow(data_utm), 0L))
  expect_length(prep$data$flag_spde, 0L)
  expect_length(prep$data$matern_range_flag, 0L)
  expect_length(prep$data$matern_sigma_flag, 0L)
})

test_that("observation-SD mapping is honored", {
  parameters <- intCPUEextra:::.make_parameters_intCPUE(
    n_a = 2L,
    n_t = 3L,
    n_v = 1L,
    n_f = 2L,
    n_s = 4L,
    n_s_flag = 4L,
    K_smooth_catch = 0L,
    n_smooth_catch = 0L,
    sum_k_catch = 0L,
    K_smooth_pop = 0L,
    n_smooth_pop = 0L,
    sum_k_pop = 0L
  )

  map <- intCPUEextra:::.make_map_intCPUE(
    parameters = parameters,
    n_f = 2L,
    pop_spatiotemporal_type = "rw",
    q_diffs_time = "off",
    obs_sd_flag = "flag",
    obs_sd_area = "area"
  )

  expect_true(all(is.na(map$ln_sd)))
  expect_null(map$ln_sd_flag)

  shared_flag_area_map <- intCPUEextra:::.make_map_intCPUE(
    parameters = parameters,
    n_f = 2L,
    pop_spatiotemporal_type = "rw",
    q_diffs_time = "off",
    obs_sd_flag = "shared",
    obs_sd_area = "area"
  )

  expect_true(all(is.na(shared_flag_area_map$ln_sd_flag)))
  expect_null(shared_flag_area_map[["ln_sd", exact = TRUE]])

  shared_flag_shared_area_map <- intCPUEextra:::.make_map_intCPUE(
    parameters = parameters,
    n_f = 2L,
    pop_spatiotemporal_type = "rw",
    q_diffs_time = "off",
    obs_sd_flag = "shared",
    obs_sd_area = "shared"
  )

  expect_true(all(is.na(shared_flag_shared_area_map$ln_sd_flag)))
  expect_equal(nlevels(shared_flag_shared_area_map$ln_sd), 1L)

  flag_shared_area_map <- intCPUEextra:::.make_map_intCPUE(
    parameters = parameters,
    n_f = 2L,
    pop_spatiotemporal_type = "rw",
    q_diffs_time = "off",
    obs_sd_flag = "flag",
    obs_sd_area = "shared"
  )

  expect_true(all(is.na(flag_shared_area_map$ln_sd)))
  expect_equal(as.integer(flag_shared_area_map$ln_sd_flag), c(1L, 2L, 1L, 2L))
})

test_that("RW/AR1 spatiotemporal type maps rho parameters appropriately", {
  parameters <- intCPUEextra:::.make_parameters_intCPUE(
    n_a = 2L,
    n_t = 3L,
    n_v = 1L,
    n_f = 2L,
    n_s = 4L,
    n_s_flag = 4L,
    K_smooth_catch = 0L,
    n_smooth_catch = 0L,
    sum_k_catch = 0L,
    K_smooth_pop = 0L,
    n_smooth_pop = 0L,
    sum_k_pop = 0L
  )

  rw_map <- intCPUEextra:::.make_map_intCPUE(
    parameters = parameters,
    n_f = 2L,
    pop_spatiotemporal_type = "rw",
    q_diffs_time = "off",
    obs_sd_flag = "shared"
  )
  expect_true(all(is.na(rw_map$transf_rho_1)))
  expect_true(all(is.na(rw_map$transf_rho_2)))

  ar1_map <- intCPUEextra:::.make_map_intCPUE(
    parameters = parameters,
    n_f = 2L,
    pop_spatiotemporal_type = "ar1",
    q_diffs_time = "off",
    obs_sd_flag = "shared"
  )
  expect_null(ar1_map[["transf_rho_1", exact = TRUE]])
  expect_null(ar1_map[["transf_rho_2", exact = TRUE]])
})

test_that("observation-SD and AR1 spatiotemporal settings are reported by TMB", {
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
      obs_sd_flag = "flag",
      q_diffs_time = "off",
      q_diffs_spatial = "off",
      control = list(eval.max = 200, iter.max = 200),
      silent = TRUE
    )
  )

  rep_flag_sd <- fit_flag_sd$obj$report()
  expect_equal(fit_flag_sd$settings$pop_spatiotemporal_type, "rw")
  expect_equal(as.integer(rep_flag_sd$use_pop_spatiotemporal_rw), 1L)
  expect_equal(as.integer(rep_flag_sd$use_pop_spatiotemporal_ar1), 0L)
  expect_equal(as.integer(rep_flag_sd$use_flag_sd), 1L)
  expect_length(rep_flag_sd$sd_flag, fit_flag_sd$data_tmb$n_f)

  fit_ar1 <- suppressWarnings(
    intCPUE(
      formula = cpue ~ 1,
      data_utm = data_utm,
      mesh = mesh,
      pop_spatiotemporal_type = "ar1",
      obs_sd_flag = "shared",
      q_diffs_time = "off",
      q_diffs_spatial = "off",
      control = list(eval.max = 200, iter.max = 200),
      silent = TRUE
    )
  )

  rep_ar1 <- fit_ar1$obj$report()
  expect_equal(fit_ar1$settings$pop_spatiotemporal_type, "ar1")
  expect_equal(as.integer(rep_ar1$use_pop_spatiotemporal_rw), 0L)
  expect_equal(as.integer(rep_ar1$use_pop_spatiotemporal_ar1), 1L)
  expect_equal(as.integer(rep_ar1$use_flag_sd), 0L)
})
