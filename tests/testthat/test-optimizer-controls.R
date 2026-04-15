test_that("intCPUE exposes optimizer refinement controls", {
  fmls <- formals(intCPUE)

  expect_equal(fmls$restart_max, 1L)
  expect_equal(fmls$newton_max, 2L)
  expect_equal(fmls$coord_max, 5L)

  expect_error(
    intCPUE(data_utm = data.frame(), mesh = NULL, restart_max = -1L),
    "`restart_max` must be a single non-negative integer."
  )
  expect_error(
    intCPUE(data_utm = data.frame(), mesh = NULL, newton_max = -1L),
    "`newton_max` must be a single non-negative integer."
  )
  expect_error(
    intCPUE(data_utm = data.frame(), mesh = NULL, coord_max = -1L),
    "`coord_max` must be a single non-negative integer."
  )
})

test_that(".safe_optimize can disable all post-nlminb refinement", {
  pkg_loaded <- "intCPUE" %in% loadedNamespaces()
  ns <- if (pkg_loaded) asNamespace("intCPUE") else environment(.safe_optimize)
  nlminb_calls <- 0L
  newton_calls <- 0L
  coord_calls <- 0L

  mocked_run <- function(obj, start, control) {
    nlminb_calls <<- nlminb_calls + 1L
    list(par = start, objective = 10, convergence = 0L, message = "ok")
  }
  mocked_grad <- function(obj, par) 10
  mocked_newton <- function(obj, start, grad_tol, trace = 0L) {
    newton_calls <<- newton_calls + 1L
    list(ok = TRUE, result = list(par = start, value = 9, gradient = 0), error = NULL)
  }
  mocked_coord <- function(obj, start, grad_tol, maxit = 5L, top_k = 5L) {
    coord_calls <<- coord_calls + 1L
    list(par = start, value = 8, gradient = 0, max_grad = 0, used = 1L)
  }

  if (pkg_loaded) {
    testthat::local_mocked_bindings(
      .run_nlminb_once = mocked_run,
      .max_grad_opt = mocked_grad,
      .try_newton_refine = mocked_newton,
      .try_coordinate_refine = mocked_coord,
      .package = "intCPUE",
      .env = ns
    )
  } else {
    orig_run <- get(".run_nlminb_once", envir = ns)
    orig_grad <- get(".max_grad_opt", envir = ns)
    orig_newton <- get(".try_newton_refine", envir = ns)
    orig_coord <- get(".try_coordinate_refine", envir = ns)
    on.exit(assign(".run_nlminb_once", orig_run, envir = ns), add = TRUE)
    on.exit(assign(".max_grad_opt", orig_grad, envir = ns), add = TRUE)
    on.exit(assign(".try_newton_refine", orig_newton, envir = ns), add = TRUE)
    on.exit(assign(".try_coordinate_refine", orig_coord, envir = ns), add = TRUE)
    assign(".run_nlminb_once", mocked_run, envir = ns)
    assign(".max_grad_opt", mocked_grad, envir = ns)
    assign(".try_newton_refine", mocked_newton, envir = ns)
    assign(".try_coordinate_refine", mocked_coord, envir = ns)
  }

  suppressWarnings(
    get(".safe_optimize", envir = ns)(
      obj = list(par = 0),
      control = list(),
      restart_max = 0L,
      newton_max = 0L,
      coord_max = 0L
    )
  )

  expect_equal(nlminb_calls, 1L)
  expect_equal(newton_calls, 0L)
  expect_equal(coord_calls, 0L)
})
