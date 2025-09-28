#'
#' Generic derivative of log-lik w.r.t. regression coefficients
#'
#' @param o object we want to get the derivatives of.
#' @param ... arguments to be passed to methods.
#' @rdname DllkDbeta
#' @noMd
#' @noRd
#' @export 
DllkDbeta <- function(o, ...)
    UseMethod("DllkDbeta")


#'
#' Derivatives of Hessian w.r.t. smoothing parameters
#'
#' @param o object we want to get the derivatives of.
#' @param ... arguments to be passed to methods.
#' @rdname DHessDrho
#' @noMd
#' @noRd
#' @export 
DHessDrho <- function(o, ...)
    UseMethod("DHessDrho")
