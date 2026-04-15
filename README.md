
<!-- README.md is generated from README.Rmd. Please edit that file -->

# intCPUE <img src="man/figures/logo.png" align="right" height="136" alt="intCPUE logo" />

> Integrated CPUE standardization with spatiotemporal models in TMB

<!-- badges: start -->

[![R-CMD-check](https://github.com/RujiaBi/intCPUE/workflows/R-CMD-check/badge.svg)](https://github.com/RujiaBi/intCPUE/actions)
<!-- badges: end -->

**intCPUE** is a TMB-based framework for integrated CPUE standardization
across multiple fisheries or surveys, with optional preferential
sampling correction (under development).

## Table of Contents

- [Overview](#overview)
- [Contact](#contact)
- [Citation](#citation)
- [Installation](#installation)
- [Data structure](#data-structure)
- [Coordinate projection](#coordinate-projection)
- [Build spatial mesh](#build-spatial-mesh)
- [Fit the model](#fit-the-model)
- [Getting index with bias
  correction](#getting-index-with-bias-correction)
- [Next steps](#next-steps)

## Overview

The model supports:

- spatiotemporal random fields via an SPDE mesh (using `fmesher`)

- multiple fleets / surveys through catchability components

- mgcv-style `s()` smooth terms parsed into mixed-effects from inside
  TMB

This is currently a Poisson-link delta model. The baseline population
surface is always driven by temporal effects together with spatial and
spatiotemporal random fields. Vessel effects and systematic fishery
differences are part of the core model. The user-facing switches
currently control the spatiotemporal population evolution type
(`pop_spatiotemporal_type`), the time-varying and spatial
fishery-specific catchability deviations (`q_diffs_time`,
`q_diffs_spatial`), and the observation-SD structure (`obs_sd`). Smooth
terms can be assigned either to catchability only or to the population
layer. Population-layer smooths are included both in the observation
model and in projection, so they affect the standardized index.

Here provides a minimal workflow: install → data preparation →
coordinate projection → mesh → model fitting → index extraction.

## Contact

For questions, suggestions, or collaboration, please contact:

- **Rujia Bi** — <rbi@iattc.org>

## Citation

If you use **intCPUE** in your work, please cite it as software:

> Bi, R. (2026). *intCPUE: Integrated CPUE standardization with TMB*. R
> package (v0.1.0). <https://github.com/RujiaBi/intCPUE>.

Once a paper or DOI is available, this section will be updated.

## Installation (development version)

``` r
# install.packages("remotes")
remotes::install_github("RujiaBi/intCPUE")
```

``` r
library(intCPUE)
```

## Data structure

An intCPUE model requires a data frame containing the following columns:

- `cpue` — positive catch rate

- `encounter` — encounter indicator (0 = zero catch, 1 = positive catch)

- `lon`, `lat` — geographic coordinates

- `vesid` — vessel ID (0-based)

- `tid` — time index (0-based)

- `flagid` — fishery / survey ID (0-based; 0 = reference fishery)

These columns are required. Additional columns may be included as
needed.

``` r
# Use `pcod` from `sdmTMB` as an example
data_input <- data.frame(
  "cpue" = pcod$density,
  "encounter" = pcod$present,
  "lon" = pcod$lon,
  "lat" = pcod$lat,
  "vesid" = 0,
  "tid" = as.numeric(as.factor(pcod$year))-1,
  "flagid" = 0,
  "depth" = pcod$depth
)
```

### Important

- `vesid`, `tid` and `flagid` must be 0-based contiguous integers

- `flagid` must use 0 as the reference fishery / survey

## Coordinate projection (lon/lat → UTM)

Longitude may be in either -180..180 or 0..360.

`make_utm()` automatically:

- detects longitude convention

- selects an appropriate UTM zone

- scales coordinates for numerical stability

``` r
utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
data_utm <- utm$data_utm
```

## Build spatial mesh

The mesh must be constructed using the scaled projected coordinates:
`utm_x_scale` and `utm_y_scale`

### K-means mesh

``` r
mesh <- make_mesh(data_utm, xy_cols = c("utm_x_scale", "utm_y_scale"), type = "kmeans", n_knots = 50)
plot(mesh)
```

### Cutoff mesh

``` r
mesh <- make_mesh(data_utm, xy_cols = c("utm_x_scale", "utm_y_scale"), type = "cutoff", cutoff = 0.1)
plot(mesh)
```

### Tailor mesh

``` r
mesh <- make_mesh(data_utm, xy_cols = c("utm_x_scale", "utm_y_scale"), type = "tailored",
    convex = -0.1,         # for a finer boundary
    max.edge = c(0.5, 2),   # max triangle edge length; inner and outer meshes
    offset = c(0.1, 0.5),  # inner and outer border widths
    cutoff = 0.05)
plot(mesh)
```

### Custom mesh

``` r
bnd <- INLA::inla.nonconvex.hull(cbind(data_utm$utm_x_scale, data_utm$utm_y_scale), convex = -0.1)
mesh_inla <- INLA::inla.mesh.2d(
  boundary = bnd,
  max.edge = c(0.5, 2)
)
mesh <- make_mesh(data_utm, xy_cols = c("utm_x_scale", "utm_y_scale"), mesh = mesh_inla)
plot(mesh)
```

## Fit the model

``` r
formula_catchability <- ~ s(depth)
formula_population <- ~ s(temp) + s(chl)
```

Interpretation of the baseline model:

Both encounter and positive components include:

- fixed temporal effects (`tid`, as intercepts)

- spatial random field (`omega`)

- spatiotemporal random field (`epsilon`)

`formula_catchability` adds smooth covariates that affect catchability
only.

`formula_population` adds smooth covariates that affect the latent
population surface and therefore also enter projection. When population
covariates are not uniquely defined within each extrapolation cell,
provide them explicitly through `projection_data`. If they also vary
over time, include a `tid` column in `projection_data` so the projection
covariates are matched by grid cell and time.

The legacy `formula = ...` interface is still supported and is treated
as `formula_catchability = ...`.

### `projection_data` format

If `formula_population` is used, population covariates must be available
on the extrapolation grid:

- Static population covariates: `projection_data` should contain one row
  per extrapolation grid cell, with columns `utm_x_scale`,
  `utm_y_scale`, and the covariates used in `formula_population`.
- Time-varying population covariates: `projection_data` should contain
  one row per grid cell-time combination, with columns `utm_x_scale`,
  `utm_y_scale`, `tid`, and the covariates used in `formula_population`.
- If `formula_population` mixes static and time-varying covariates, use
  the grid cell-time format and repeat the static covariate values
  across `tid`.

`tid` must use the same 0-based coding as `data_utm$tid`.

``` r
ncores <- 4
mesh <- make_mesh(data_utm, xy_cols = c("utm_x_scale", "utm_y_scale"), type = "cutoff", cutoff = 0.1)
fit <- intCPUE(
  formula_catchability = formula_catchability,
  data_utm = data_utm,
  mesh = mesh,
  q_diffs_time = "off",  
  q_diffs_spatial = "off",
  obs_sd = "shared",
  ncores = ncores
)
```

### Core model components

- `pop_spatiotemporal_type` — temporal dependence for the population
  spatiotemporal field (`"rw"` or `"ar1"`)

The current package version fits a fixed core model with:

- population spatial and spatiotemporal random fields turned on

- vessel effects turned on

- systematic fishery catchability differences turned on

### User-facing switches

- `q_diffs_time` — time-varying catchability difference

- `q_diffs_spatial` — spatial catchability difference

For the reference fishery (`flagid = 0`), the `q_diffs_*` terms are
constrained to 0.

### Observation error

- `obs_sd = "shared"` uses one lognormal observation SD across all flags

- `obs_sd = "flag"` estimates one lognormal observation SD for each flag

## Getting index with bias correction

``` r
index <- get_index(fit)
plot_index(index)
```

## Diagnostics

``` r
check_convergence(fit)
calc_marginal_aic(fit)

pred <- get_predicted(fit, data = data_utm)
plots <- plot_residuals(pred, observed_col = "cpue")

plots$observed_predicted
plots$spatial_residual

plot_anisotropy(fit)
```

## Next steps

- preferential sampling correction

- length frequency part
