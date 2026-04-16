.template_id_intCPUE <- function(q_diffs_time = c("on", "off"),
                                 q_diffs_spatial = c("on", "off")) {
  q_diffs_time <- match.arg(q_diffs_time)
  q_diffs_spatial <- match.arg(q_diffs_spatial)

  paste0(
    "intCPUE_q",
    if (q_diffs_time == "on") "1" else "0",
    if (q_diffs_spatial == "on") "1" else "0"
  )
}

.ensure_template_dll_intCPUE <- function(q_diffs_time = c("on", "off"),
                                         q_diffs_spatial = c("on", "off")) {
  template_id <- .template_id_intCPUE(q_diffs_time, q_diffs_spatial)

  dll_dir <- system.file("tmb", package = "intCPUEextra")
  if (!nzchar(dll_dir)) {
    stop("Could not locate installed intCPUEextra template DLLs in inst/tmb.", call. = FALSE)
  }

  dll_path <- TMB::dynlib(file.path(dll_dir, template_id))
  if (!file.exists(dll_path)) {
    stop(
      "Missing precompiled template DLL for ", template_id,
      ". Reinstall intCPUEextra so installation can build the template set.",
      call. = FALSE
    )
  }

  loaded_names <- names(getLoadedDLLs())
  if (!(template_id %in% loaded_names)) {
    dyn.load(dll_path)
  }

  template_id
}
