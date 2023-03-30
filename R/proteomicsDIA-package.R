#' @description
#' A package for the differential proteomics analysis with DIA data. To learn more,
#' start with the tutorial vignette: `browseVignettes(package = "proteoDA")`
#' @keywords internal
#' @importFrom rlang .data
"_PACKAGE" # nocov

NULL # nocov


#' Dummy function for statmod import
#'
#' The limma functions we use need to have statmod installed, but limma
#' only suggests statmod, it doesn't depend on it. So, users can get halfway
#' through our pipeline, then have to install statmod to run a core function.
#' This also causes issues with running vignettes. So, we'll have a dummy
#' function that includes a statmod import.
#'
#' @return NULL
dummy_statmod <- function() { # nocov start
  x <- statmod::vecmat()

  NULL
} # nocov end
