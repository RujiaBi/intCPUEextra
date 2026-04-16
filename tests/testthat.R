# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

Sys.setenv(TESTTHAT_CPUS = "1")

library(testthat)
library(intCPUEextra)

test_path <- if (dir.exists("tests/testthat")) "tests/testthat" else "testthat"
test_files <- list.files(test_path, pattern = "\\.[Rr]$", full.names = TRUE)

for (path in sort(test_files)) {
  source(path, chdir = FALSE, local = new.env(parent = globalenv()))
}
