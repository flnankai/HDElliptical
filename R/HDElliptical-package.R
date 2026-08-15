#' HDElliptical: high-dimensional elliptical methods
#'
#' `HDElliptical` provides R and Rcpp implementations accompanying
#' *High-Dimensional Data Analysis for Elliptically Symmetric Distributions*.
#' Functions use observations in rows and variables in columns throughout.
#'
#' @keywords internal
#' @useDynLib HDElliptical, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"

# Internal string concatenation helper used only in error messages.
`%+%` <- function(lhs, rhs) paste0(lhs, rhs)
